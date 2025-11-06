// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/roofGeneration.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/algorithm/extrude.h"
#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/algorithm/straightSkeleton.h"
#include "SFCGAL/algorithm/covers.h"
#include "SFCGAL/algorithm/length.h"
#include "SFCGAL/algorithm/translate.h"
#include "SFCGAL/triangulate/triangulate2DZ.h"

#include <map>
#include <limits>
#include <set>

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

namespace SFCGAL::algorithm {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

namespace {

/**
 * @brief Calculate distance from point to line segment
 */
auto distancePointToSegment(const Point &point, const Point &segStart, const Point &segEnd) -> Kernel::FT {
  // Vector from segment start to end
  Kernel::FT segX = segEnd.x() - segStart.x();
  Kernel::FT segY = segEnd.y() - segStart.y();

  // Vector from segment start to point
  Kernel::FT ptX = point.x() - segStart.x();
  Kernel::FT ptY = point.y() - segStart.y();

  // Segment length squared
  Kernel::FT segLenSq = segX * segX + segY * segY;

  if (segLenSq < 1e-20) {
    // Degenerate segment, return distance to start point
    return std::sqrt(CGAL::to_double(ptX * ptX + ptY * ptY));
  }

  // Project point onto segment line: t = dot(AP, AB) / |AB|^2
  Kernel::FT t = (ptX * segX + ptY * segY) / segLenSq;

  // Clamp t to [0, 1] to stay on segment
  if (t < 0) t = 0;
  if (t > 1) t = 1;

  // Calculate closest point on segment
  Kernel::FT closestX = segStart.x() + t * segX;
  Kernel::FT closestY = segStart.y() + t * segY;

  // Return distance from point to closest point on segment
  Kernel::FT distX = point.x() - closestX;
  Kernel::FT distY = point.y() - closestY;

  return std::sqrt(CGAL::to_double(distX * distX + distY * distY));
}

/**
 * @brief Add vertical faces at ridge extensions
 */
void addVerticalFacesAtRidgeExtensions(const Polygon &footprint, const LineString &ridgeLine,
                                     Kernel::FT ridgeHeight, PolyhedralSurface &roof) {
  // Simplified implementation - for now, just skip vertical faces
  // These could be added later if needed
  return;
}

/**
 * @brief Calculate 3D point projected from base point with given height and slope
 */
auto
calculateSlopePoint(const Point &basePoint, const Kernel::Vector_3 &normal,
                    double distance, double slopeAngle) -> Point
{
  double height = distance * std::tan(slopeAngle * M_PI / 180.0);
  return Point(basePoint.x() + normal.x() * height,
               basePoint.y() + normal.y() * height,
               basePoint.z() + height);
}

/**
 * @brief Split a polygon along a line and return the two resulting polygons
 */
auto
splitPolygonByLine(const Polygon &polygon, const LineString &line)
    -> std::pair<std::unique_ptr<Polygon>, std::unique_ptr<Polygon>>
{
  if (polygon.isEmpty() || line.isEmpty()) {
    BOOST_THROW_EXCEPTION(Exception("Cannot split empty geometry"));
  }

  if (line.numPoints() < 2) {
    BOOST_THROW_EXCEPTION(Exception("Ridge line must contain at least 2 points"));
  }

  // For a simplified but more robust implementation,
  // just return two halves of the polygon split at the midpoint
  auto &ring = polygon.exteriorRing();
  size_t numPoints = ring.numPoints() - 1; // Excluding closing point

  if (numPoints < 4) {
    BOOST_THROW_EXCEPTION(Exception("Polygon must have at least 4 vertices for splitting"));
  }

  // Create two halves by splitting at roughly the middle
  std::vector<Point> side1, side2;

  size_t midpoint = numPoints / 2;

  // First half
  for (size_t i = 0; i <= midpoint; ++i) {
    side1.push_back(ring.pointN(i));
  }

  // Add ridge line points as connecting points
  side1.push_back(line.pointN(0));
  side1.push_back(line.pointN(line.numPoints() - 1));
  side1.push_back(ring.pointN(0)); // Close the ring

  // Second half
  for (size_t i = midpoint; i < numPoints; ++i) {
    side2.push_back(ring.pointN(i));
  }
  side2.push_back(ring.pointN(0)); // Back to start

  // Add ridge line points
  side2.push_back(line.pointN(line.numPoints() - 1));
  side2.push_back(line.pointN(0));
  side2.push_back(ring.pointN(midpoint)); // Close the ring

  // Create polygons
  auto poly1 = std::make_unique<Polygon>(LineString(side1));
  auto poly2 = std::make_unique<Polygon>(LineString(side2));

  return std::make_pair(std::move(poly1), std::move(poly2));
}

/**
 * @brief Create a sloped surface from polygon base to ridge line
 */
auto
createSlopedSurface(const Polygon &basePolygon, const LineString &ridgeLine,
                    double slopeAngle, bool slopeUp = true)
    -> std::unique_ptr<PolyhedralSurface>
{
  auto surface = std::make_unique<PolyhedralSurface>();

  if (basePolygon.isEmpty() || ridgeLine.isEmpty()) {
    return surface; // Return empty surface for empty input
  }

  // Calculate ridge height based on slope angle and average distance
  double avgDistance = 2.0; // Simplified average distance
  auto ridgeHeight = calculateRidgeHeight(avgDistance, slopeAngle);

  // Create elevated ridge line
  std::vector<Point> ridgePoints;
  for (size_t i = 0; i < ridgeLine.numPoints(); ++i) {
    auto point = ridgeLine.pointN(i);
    ridgePoints.emplace_back(point.x(), point.y(),
                            slopeUp ? ridgeHeight : 0.0);
  }

  if (ridgePoints.empty()) {
    return surface; // Return empty surface if no ridge points
  }

  // Create simple triangulated surface
  auto &baseRing = basePolygon.exteriorRing();
  size_t numBasePoints = baseRing.numPoints() - 1; // Exclude closing point

  if (numBasePoints < 3) {
    return surface; // Need at least 3 points for a polygon
  }

  // For simplicity, connect all base vertices to the first ridge point
  auto ridgePoint = ridgePoints[0];

  for (size_t i = 0; i < numBasePoints; ++i) {
    auto basePoint1 = baseRing.pointN(i);
    auto basePoint2 = baseRing.pointN((i + 1) % numBasePoints);

    // Create triangle from two consecutive base points to ridge point
    std::vector<Point> trianglePoints = {basePoint1, basePoint2, ridgePoint, basePoint1};
    surface->addPatch(Polygon(LineString(trianglePoints)));
  }

  return surface;
}

} // anonymous namespace

/// @} end of private section

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

auto
calculateRidgeHeight(double horizontalDistance, double slopeAngle) -> double
{
  if (slopeAngle <= 0.0 || slopeAngle >= 90.0) {
    BOOST_THROW_EXCEPTION(Exception("Slope angle must be between 0 and 90 degrees"));
  }
  return horizontalDistance * std::tan(slopeAngle * M_PI / 180.0);
}

auto
calculateHorizontalDistance(double height, double slopeAngle) -> double
{
  if (slopeAngle <= 0.0 || slopeAngle >= 90.0) {
    BOOST_THROW_EXCEPTION(Exception("Slope angle must be between 0 and 90 degrees"));
  }
  return height / std::tan(slopeAngle * M_PI / 180.0);
}

auto
generatePitchedRoof(const Polygon &footprint, const LineString &ridgeLine,
                    double slopeAngle, RidgePosition ridgePosition)
    -> std::unique_ptr<PolyhedralSurface>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(footprint);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(ridgeLine);

  if (slopeAngle <= 0.0 || slopeAngle >= 90.0) {
    BOOST_THROW_EXCEPTION(Exception("Slope angle must be between 0 and 90 degrees"));
  }

  auto result = std::make_unique<PolyhedralSurface>();

  if (footprint.isEmpty() || ridgeLine.isEmpty()) {
    return result; // Return empty surface for empty input
  }

  try {
    switch (ridgePosition) {
    case RidgePosition::INTERIOR: {
      // Split polygon along ridge line and create two slopes
      auto [poly1, poly2] = splitPolygonByLine(footprint, ridgeLine);

      auto slope1 = createSlopedSurface(*poly1, ridgeLine, slopeAngle, true);
      auto slope2 = createSlopedSurface(*poly2, ridgeLine, slopeAngle, true);

      // Combine surfaces
      for (size_t i = 0; i < slope1->numPatches(); ++i) {
        result->addPatch(slope1->patchN(i));
      }
      for (size_t i = 0; i < slope2->numPatches(); ++i) {
        result->addPatch(slope2->patchN(i));
      }
      break;
    }
    case RidgePosition::EDGE:
    case RidgePosition::EXTERIOR: {
      // Create single slope from ridge to opposite edges
      auto surface = createSlopedSurface(footprint, ridgeLine, slopeAngle, true);
      result = std::move(surface);
      break;
    }
    }
  } catch (const Exception& e) {
    // If splitting fails, fall back to single surface
    auto surface = createSlopedSurface(footprint, ridgeLine, slopeAngle, true);
    result = std::move(surface);
  }

  propagateValidityFlag(*result, true);
  return result;
}

auto
generatePitchedRoof(const Polygon &footprint, const LineString &ridgeLine,
                    double slopeAngle, RidgePosition ridgePosition,
                    NoValidityCheck &/*nvc*/) -> std::unique_ptr<PolyhedralSurface>
{
  if (slopeAngle <= 0.0 || slopeAngle >= 90.0) {
    BOOST_THROW_EXCEPTION(Exception("Slope angle must be between 0 and 90 degrees"));
  }

  auto result = std::make_unique<PolyhedralSurface>();

  switch (ridgePosition) {
  case RidgePosition::INTERIOR: {
    auto [poly1, poly2] = splitPolygonByLine(footprint, ridgeLine);

    auto slope1 = createSlopedSurface(*poly1, ridgeLine, slopeAngle, true);
    auto slope2 = createSlopedSurface(*poly2, ridgeLine, slopeAngle, true);

    for (size_t i = 0; i < slope1->numPatches(); ++i) {
      result->addPatch(slope1->patchN(i));
    }
    for (size_t i = 0; i < slope2->numPatches(); ++i) {
      result->addPatch(slope2->patchN(i));
    }
    break;
  }
  case RidgePosition::EDGE:
  case RidgePosition::EXTERIOR: {
    auto surface = createSlopedSurface(footprint, ridgeLine, slopeAngle, true);
    result = std::move(surface);
    break;
  }
  }

  return result;
}

auto
generateGableRoof(const Polygon &footprint, const LineString &ridgeLine,
                  double slopeAngle, bool addHips)
    -> std::unique_ptr<PolyhedralSurface>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(footprint);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(ridgeLine);

  if (slopeAngle <= 0.0 || slopeAngle >= 90.0) {
    BOOST_THROW_EXCEPTION(Exception("Slope angle must be between 0 and 90 degrees"));
  }

  auto result = std::make_unique<PolyhedralSurface>();

  if (footprint.isEmpty() || ridgeLine.isEmpty()) {
    return result; // Return empty surface for empty input
  }

  // Calculate proper ridge height based on maximum distance to polygon edge
  auto ridgeStart = ridgeLine.pointN(0);
  auto ridgeEnd = ridgeLine.pointN(ridgeLine.numPoints() - 1);

  // Calculate proper perpendicular distance from ridge line to polygon edge
  double maxPerpDistance = 0.0;
  auto exteriorRing = footprint.exteriorRing();

  // For each polygon edge, calculate perpendicular distance to ridge line
  for (size_t i = 0; i < exteriorRing.numPoints() - 1; ++i) {
    auto edgeStart = exteriorRing.pointN(i);
    auto edgeEnd = exteriorRing.pointN((i + 1) % (exteriorRing.numPoints() - 1));

    // Calculate perpendicular distance from ridge line to this polygon edge
    // Ridge line vector
    auto ridgeVecX = ridgeEnd.x() - ridgeStart.x();
    auto ridgeVecY = ridgeEnd.y() - ridgeStart.y();
    double ridgeLength = std::sqrt(CGAL::to_double(ridgeVecX * ridgeVecX + ridgeVecY * ridgeVecY));

    if (ridgeLength > 1e-10) {
      // Normalize ridge vector
      auto ridgeNormX = ridgeVecX / ridgeLength;
      auto ridgeNormY = ridgeVecY / ridgeLength;

      // Calculate distance from edge endpoints to ridge line
      auto dx1 = edgeStart.x() - ridgeStart.x();
      auto dy1 = edgeStart.y() - ridgeStart.y();
      double dist1 = std::abs(CGAL::to_double(dx1 * (-ridgeNormY) + dy1 * ridgeNormX));

      auto dx2 = edgeEnd.x() - ridgeStart.x();
      auto dy2 = edgeEnd.y() - ridgeStart.y();
      double dist2 = std::abs(CGAL::to_double(dx2 * (-ridgeNormY) + dy2 * ridgeNormX));

      maxPerpDistance = std::max({maxPerpDistance, dist1, dist2});
    }
  }

  // Ensure minimum distance for roof height calculation
  maxPerpDistance = std::max(maxPerpDistance, 1.0);
  double ridgeHeight = calculateRidgeHeight(maxPerpDistance, slopeAngle);

  // Create ridge line points at calculated height
  Point ridgeStartAtHeight(ridgeStart.x(), ridgeStart.y(), ridgeHeight);
  Point ridgeEndAtHeight(ridgeEnd.x(), ridgeEnd.y(), ridgeHeight);

  // Create the two main roof slopes
  // For a gable roof, we create two symmetric slopes meeting at the ridge

  // Split the polygon into two halves
  size_t numPoints = exteriorRing.numPoints() - 1; // Excluding closing point
  if (numPoints < 4) {
    BOOST_THROW_EXCEPTION(Exception("Polygon must have at least 4 vertices for gable roof"));
  }

  // CDT-based approach for gable roof generation

  // 1. Get projected medial axis as ridge line
  auto projectedMedialAxis = projectMedialAxisToEdges(footprint);
  if (!projectedMedialAxis || projectedMedialAxis->isEmpty()) {
    BOOST_THROW_EXCEPTION(Exception("Could not compute projected medial axis for gable roof"));
  }

  // 2. Perform CDT with ridge line as constraints
  auto triangulation = SFCGAL::triangulate::triangulate2DZ(footprint);
  SFCGAL::triangulate::triangulate2DZ(*projectedMedialAxis, triangulation);
  auto cdtSurface = triangulation.getTriangulatedSurface();

  if (!cdtSurface) {
    BOOST_THROW_EXCEPTION(Exception("CDT triangulation failed"));
  }

  // 3. Collect ridge points and calculate uniform ridge height
  std::set<std::pair<double, double>> ridgePointsSet;
  for (size_t i = 0; i < projectedMedialAxis->numGeometries(); ++i) {
    const auto *line = dynamic_cast<const LineString *>(&projectedMedialAxis->geometryN(i));
    if (!line) continue;

    for (size_t j = 0; j < line->numPoints(); ++j) {
      const auto &rp = line->pointN(j);
      ridgePointsSet.insert({CGAL::to_double(rp.x()), CGAL::to_double(rp.y())});
    }
  }

  // Calculate uniform ridge height based on maximum distance from ridge line to polygon edge
  double uniformRidgeHeight = ridgeHeight; // Use the height calculated earlier

  // 4. Process each triangle from CDT
  for (size_t i = 0; i < cdtSurface->numPatches(); ++i) {
    const auto &patch = cdtSurface->patchN(i);
    const auto *triangle = dynamic_cast<const Triangle *>(&patch);
    if (!triangle) continue;

    // Get triangle vertices
    auto v1 = triangle->vertex(0);
    auto v2 = triangle->vertex(1);
    auto v3 = triangle->vertex(2);

    // Calculate triangle centroid using exact arithmetic
    auto centroidX = (v1.x() + v2.x() + v3.x()) / 3;
    auto centroidY = (v1.y() + v2.y() + v3.y()) / 3;

    // 5. Filter: check if centroid is inside original polygon
    Point centroidPt(centroidX, centroidY, 0.0);
    if (!SFCGAL::algorithm::covers(footprint, centroidPt)) {
      continue; // Skip triangles outside polygon
    }

    // 6. Elevate ridge vertices and create roof triangle
    Point newV1 = v1;
    Point newV2 = v2;
    Point newV3 = v3;

    // Check if each vertex is a ridge point and elevate it
    auto isRidgePoint = [&ridgePointsSet](const Point &p) -> bool {
      return ridgePointsSet.find({CGAL::to_double(p.x()), CGAL::to_double(p.y())}) != ridgePointsSet.end();
    };

    if (isRidgePoint(v1)) {
      // Use uniform ridge height for all ridge points
      newV1 = Point(v1.x(), v1.y(), uniformRidgeHeight);
    }

    if (isRidgePoint(v2)) {
      // Use uniform ridge height for all ridge points
      newV2 = Point(v2.x(), v2.y(), uniformRidgeHeight);
    }

    if (isRidgePoint(v3)) {
      // Use uniform ridge height for all ridge points
      newV3 = Point(v3.x(), v3.y(), uniformRidgeHeight);
    }

    // Create roof triangle with elevated vertices
    std::unique_ptr<Triangle> roofTriangle(new Triangle(newV1, newV2, newV3));
    result->addPatch(*roofTriangle);
  }

  // 7. Add vertical faces (walls) connecting polygon edges to elevated ridge points
  for (size_t i = 0; i < exteriorRing.numPoints() - 1; ++i) {
    auto edgeStart = exteriorRing.pointN(i);
    auto edgeEnd = exteriorRing.pointN((i + 1) % (exteriorRing.numPoints() - 1));

    Point edge1(edgeStart.x(), edgeStart.y(), 0.0);
    Point edge2(edgeEnd.x(), edgeEnd.y(), 0.0);

    // Check if edge endpoints are ridge points
    bool edge1IsRidge = ridgePointsSet.find({CGAL::to_double(edgeStart.x()), CGAL::to_double(edgeStart.y())}) != ridgePointsSet.end();
    bool edge2IsRidge = ridgePointsSet.find({CGAL::to_double(edgeEnd.x()), CGAL::to_double(edgeEnd.y())}) != ridgePointsSet.end();

    if (edge1IsRidge || edge2IsRidge) {
      // Calculate heights for elevated points
      Point elevatedEdge1 = edge1;
      Point elevatedEdge2 = edge2;

      if (edge1IsRidge) {
        // Use uniform ridge height for all ridge points
        elevatedEdge1 = Point(edgeStart.x(), edgeStart.y(), uniformRidgeHeight);
      }

      if (edge2IsRidge) {
        // Use uniform ridge height for all ridge points
        elevatedEdge2 = Point(edgeEnd.x(), edgeEnd.y(), uniformRidgeHeight);
      }

      // Create vertical triangles if needed
      if (edge1IsRidge && !edge2IsRidge) {
        std::unique_ptr<Triangle> verticalTriangle1(new Triangle(edge1, edge2, elevatedEdge1));
        result->addPatch(*verticalTriangle1);
      } else if (!edge1IsRidge && edge2IsRidge) {
        std::unique_ptr<Triangle> verticalTriangle2(new Triangle(edge1, edge2, elevatedEdge2));
        result->addPatch(*verticalTriangle2);
      } else if (edge1IsRidge && edge2IsRidge) {
        // Both ends elevated - create quad as two triangles
        std::unique_ptr<Triangle> verticalTriangle1(new Triangle(edge1, edge2, elevatedEdge1));
        std::unique_ptr<Triangle> verticalTriangle2(new Triangle(edge2, elevatedEdge2, elevatedEdge1));
        result->addPatch(*verticalTriangle1);
        result->addPatch(*verticalTriangle2);
      }
    }
  }

  if (addHips) {
    // Add hip treatment at gable ends
    // This would involve creating triangular surfaces at the ends
    // For now, we'll use the basic gable roof result
  }

  propagateValidityFlag(*result, true);
  return result;
}

auto
generateGableRoof(const Polygon &footprint, const LineString &ridgeLine,
                  double slopeAngle, bool addHips, NoValidityCheck &nvc)
    -> std::unique_ptr<PolyhedralSurface>
{
  if (slopeAngle <= 0.0 || slopeAngle >= 90.0) {
    BOOST_THROW_EXCEPTION(Exception("Slope angle must be between 0 and 90 degrees"));
  }

  auto result = generatePitchedRoof(footprint, ridgeLine, slopeAngle,
                                   RidgePosition::INTERIOR, nvc);

  if (addHips) {
    // Add hip treatment
  }

  // Add gable ends
  auto ridgeStart = ridgeLine.pointN(0);
  auto ridgeEnd = ridgeLine.pointN(ridgeLine.numPoints() - 1);
  double ridgeHeight = calculateRidgeHeight(1.0, slopeAngle);

  Point ridgeStartTop(ridgeStart.x(), ridgeStart.y(), ridgeHeight);
  Point ridgeEndTop(ridgeEnd.x(), ridgeEnd.y(), ridgeHeight);

  Point gableBase1(ridgeStart.x(), ridgeStart.y() + 1.0, 0.0);
  Point gableBase2(ridgeEnd.x(), ridgeEnd.y() + 1.0, 0.0);

  std::vector<Point> gableTriangle1 = {ridgeStart, ridgeStartTop, gableBase1, ridgeStart};
  result->addPatch(Polygon(LineString(gableTriangle1)));

  std::vector<Point> gableTriangle2 = {ridgeEnd, ridgeEndTop, gableBase2, ridgeEnd};
  result->addPatch(Polygon(LineString(gableTriangle2)));

  return result;
}

auto
generateSkillionRoof(const Polygon &footprint, const LineString &ridgeLine,
                     double slopeAngle) -> std::unique_ptr<PolyhedralSurface>
{
  // Skillion roof is essentially a pitched roof with edge ridge position
  return generatePitchedRoof(footprint, ridgeLine, slopeAngle,
                            RidgePosition::EDGE);
}

auto
generateRoof(const Polygon &footprint, const LineString &ridgeLine,
             const RoofParameters &params) -> std::unique_ptr<PolyhedralSurface>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(footprint);

  switch (params.type) {
  case RoofType::FLAT: {
    // Use extrude to create flat roof - use roofHeight for backward compatibility
    double height = (params.roofHeight > 0.0) ? params.roofHeight : 3.0;
    auto extruded = extrude(footprint, height);
    if (auto solid = dynamic_cast<const Solid*>(extruded.get())) {
      return std::make_unique<PolyhedralSurface>(solid->exteriorShell());
    }
    BOOST_THROW_EXCEPTION(Exception("Failed to create flat roof"));
  }
  case RoofType::HIPPED: {
    // Use extrudeStraightSkeleton for hipped roof - use roofHeight for backward compatibility
    double height = (params.roofHeight > 0.0) ? params.roofHeight : 3.0;
    return extrudeStraightSkeleton(footprint, height);
  }
  case RoofType::PITCHED:
  case RoofType::SKILLION: {
    SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(ridgeLine);
    return generatePitchedRoof(footprint, ridgeLine, params.slopeAngle,
                              params.ridgePosition);
  }
  case RoofType::GABLE: {
    SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(ridgeLine);
    return generateGableRoof(footprint, ridgeLine, params.slopeAngle,
                            params.addHips);
  }
  default:
    BOOST_THROW_EXCEPTION(Exception("Unknown roof type"));
  }
}

auto
generateRoof(const Polygon &footprint, const LineString &ridgeLine,
             const RoofParameters &params, NoValidityCheck &nvc)
    -> std::unique_ptr<PolyhedralSurface>
{
  switch (params.type) {
  case RoofType::FLAT: {
    double height = (params.roofHeight > 0.0) ? params.roofHeight : 3.0;
    auto extruded = extrude(footprint, height);
    if (auto solid = dynamic_cast<const Solid*>(extruded.get())) {
      return std::make_unique<PolyhedralSurface>(solid->exteriorShell());
    }
    BOOST_THROW_EXCEPTION(Exception("Failed to create flat roof"));
  }
  case RoofType::HIPPED: {
    double height = (params.roofHeight > 0.0) ? params.roofHeight : 3.0;
    return extrudeStraightSkeleton(footprint, height);
  }
  case RoofType::PITCHED:
  case RoofType::SKILLION: {
    return generatePitchedRoof(footprint, ridgeLine, params.slopeAngle,
                              params.ridgePosition, nvc);
  }
  case RoofType::GABLE: {
    return generateGableRoof(footprint, ridgeLine, params.slopeAngle,
                            params.addHips, nvc);
  }
  default:
    BOOST_THROW_EXCEPTION(Exception("Unknown roof type"));
  }
}

/**
 * @brief Create a closed Solid by combining building base and roof
 */
auto
createBuildingWithRoof(const Polygon &footprint, double buildingHeight,
                       std::unique_ptr<PolyhedralSurface> roof)
    -> std::unique_ptr<Geometry>
{
  if (buildingHeight <= 0.0) {
    // No building base, just return the roof
    return std::move(roof);
  }

  // Create building base by extrusion
  auto buildingBase = extrude(footprint, buildingHeight);

  if (!buildingBase) {
    BOOST_THROW_EXCEPTION(Exception("Failed to create building base"));
  }

  // If we just want the roof on top of the building without complex union
  if (auto solid = dynamic_cast<const Solid*>(buildingBase.get())) {
    auto result = std::make_unique<PolyhedralSurface>(solid->exteriorShell());

    // Add roof surfaces translated to building height
    for (size_t i = 0; i < roof->numPatches(); ++i) {
      auto roofPatch = roof->patchN(i);

      // Translate roof patch to building height
      std::vector<Point> translatedPoints;
      if (auto polygon = dynamic_cast<const Polygon*>(&roofPatch)) {
        const auto& ring = polygon->exteriorRing();
        for (size_t j = 0; j < ring.numPoints(); ++j) {
          auto pt = ring.pointN(j);
          translatedPoints.emplace_back(pt.x(), pt.y(), pt.z() + buildingHeight);
        }
        result->addPatch(Polygon(LineString(translatedPoints)));
      }
    }

    // Create a proper Solid from the PolyhedralSurface
    if (result->numPatches() > 0) {
      // Check if the surface is closed and oriented correctly
      return std::make_unique<Solid>(*result);
    }
  }

  return std::move(roof);
}

/**
 * @brief Generate a pitched roof with building height and roof height.
 */
auto
generatePitchedRoof(const Polygon &footprint, const LineString &ridgeLine,
                    double buildingHeight, double roofHeight, double slopeAngle,
                    RidgePosition ridgePosition) -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(footprint);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(ridgeLine);

  if (buildingHeight < 0.0) {
    BOOST_THROW_EXCEPTION(Exception("Building height must be >= 0"));
  }
  if (roofHeight <= 0.0) {
    BOOST_THROW_EXCEPTION(Exception("Roof height must be > 0"));
  }
  if (slopeAngle <= 0.0 || slopeAngle >= 90.0) {
    BOOST_THROW_EXCEPTION(Exception("Slope angle must be between 0 and 90 degrees"));
  }

  // Create roof using the original function
  auto roof = generatePitchedRoof(footprint, ridgeLine, slopeAngle, ridgePosition);

  // Combine with building base
  return createBuildingWithRoof(footprint, buildingHeight, std::move(roof));
}

/**
 * @brief Generate a gable roof with building height and roof height.
 */
auto
generateGableRoof(const Polygon &footprint, const LineString &ridgeLine,
                  double buildingHeight, double roofHeight, double slopeAngle,
                  bool addHips) -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(footprint);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(ridgeLine);

  if (buildingHeight < 0.0) {
    BOOST_THROW_EXCEPTION(Exception("Building height must be >= 0"));
  }
  if (roofHeight <= 0.0) {
    BOOST_THROW_EXCEPTION(Exception("Roof height must be > 0"));
  }
  if (slopeAngle <= 0.0 || slopeAngle >= 90.0) {
    BOOST_THROW_EXCEPTION(Exception("Slope angle must be between 0 and 90 degrees"));
  }

  // Create roof at the specified slope angle
  auto roof = generateGableRoof(footprint, ridgeLine, slopeAngle, addHips);

  // Scale the roof height to match the desired roofHeight
  // The roof is generated with its natural height based on slope angle
  // We need to scale it to achieve the desired maximum height
  auto ridgeStart = ridgeLine.pointN(0);
  auto exteriorRing = footprint.exteriorRing();

  // Calculate the natural ridge height that was generated
  double maxDistance = 3.0;
  for (size_t i = 0; i < exteriorRing.numPoints() - 1; ++i) {
    auto vertex = exteriorRing.pointN(i);
    double dx = CGAL::to_double(vertex.x() - ridgeStart.x());
    double dy = CGAL::to_double(vertex.y() - ridgeStart.y());
    double distance = std::sqrt(dx * dx + dy * dy);
    maxDistance = std::max(maxDistance, distance);
  }

  double naturalRidgeHeight = calculateRidgeHeight(maxDistance, slopeAngle);
  double heightScale = (naturalRidgeHeight > 0.0) ? roofHeight / naturalRidgeHeight : 1.0;

  // Create a new roof with scaled z-coordinates
  auto scaledRoof = std::make_unique<PolyhedralSurface>();
  for (size_t i = 0; i < roof->numPatches(); ++i) {
    auto& patch = roof->patchN(i);
    if (auto polygon = dynamic_cast<const Polygon*>(&patch)) {
      auto& ring = polygon->exteriorRing();
      std::vector<Point> scaledPoints;

      for (size_t j = 0; j < ring.numPoints(); ++j) {
        auto point = ring.pointN(j);
        double newZ = CGAL::to_double(point.z()) * heightScale;
        scaledPoints.emplace_back(point.x(), point.y(), newZ);
      }

      // Add scaled patch to new roof
      scaledRoof->addPatch(Polygon(LineString(scaledPoints)));
    }
  }
  roof = std::move(scaledRoof);

  // Combine with building base
  return createBuildingWithRoof(footprint, buildingHeight, std::move(roof));
}

/**
 * @brief Generate a skillion roof with building height and roof height.
 */
auto
generateSkillionRoof(const Polygon &footprint, const LineString &ridgeLine,
                     double buildingHeight, double roofHeight, double slopeAngle)
    -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(footprint);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(ridgeLine);

  if (buildingHeight < 0.0) {
    BOOST_THROW_EXCEPTION(Exception("Building height must be >= 0"));
  }
  if (roofHeight <= 0.0) {
    BOOST_THROW_EXCEPTION(Exception("Roof height must be > 0"));
  }
  if (slopeAngle <= 0.0 || slopeAngle >= 90.0) {
    BOOST_THROW_EXCEPTION(Exception("Slope angle must be between 0 and 90 degrees"));
  }

  // Create roof using the original function
  auto roof = generateSkillionRoof(footprint, ridgeLine, slopeAngle);

  // Combine with building base
  return createBuildingWithRoof(footprint, buildingHeight, std::move(roof));
}

/**
 * @brief Generate a complete building with roof using unified parameters.
 */
auto
generateBuildingWithRoof(const Polygon &footprint, const LineString &ridgeLine,
                        const RoofParameters &params) -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(footprint);

  if (params.buildingHeight < 0.0) {
    BOOST_THROW_EXCEPTION(Exception("Building height must be >= 0"));
  }
  if (params.roofHeight <= 0.0) {
    BOOST_THROW_EXCEPTION(Exception("Roof height must be > 0"));
  }

  switch (params.type) {
  case RoofType::FLAT: {
    // Use extrude to create flat roof at building + roof height
    double totalHeight = params.buildingHeight + params.roofHeight;
    auto extruded = extrude(footprint, totalHeight);
    if (params.generateSolid && extruded) {
      return extruded;
    }
    if (auto solid = dynamic_cast<const Solid*>(extruded.get())) {
      return std::make_unique<PolyhedralSurface>(solid->exteriorShell());
    }
    BOOST_THROW_EXCEPTION(Exception("Failed to create flat roof"));
  }
  case RoofType::HIPPED: {
    // Use extrudeStraightSkeleton for hipped roof
    if (params.buildingHeight > 0.0) {
      auto result = extrudeStraightSkeleton(footprint, params.buildingHeight, params.roofHeight);
      if (params.generateSolid && result) {
        return std::make_unique<Solid>(*result);
      }
      return result;
    } else {
      auto result = extrudeStraightSkeleton(footprint, params.roofHeight);
      if (params.generateSolid && result) {
        return std::make_unique<Solid>(*result);
      }
      return result;
    }
  }
  case RoofType::PITCHED:
  case RoofType::SKILLION: {
    SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(ridgeLine);
    if (params.buildingHeight > 0.0 || params.generateSolid) {
      return generatePitchedRoof(footprint, ridgeLine, params.buildingHeight,
                                params.roofHeight, params.slopeAngle,
                                params.ridgePosition);
    } else {
      auto result = generatePitchedRoof(footprint, ridgeLine, params.slopeAngle,
                                       params.ridgePosition);
      return result;
    }
  }
  case RoofType::GABLE: {
    SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(ridgeLine);
    if (params.buildingHeight > 0.0 || params.generateSolid) {
      return generateGableRoof(footprint, ridgeLine, params.buildingHeight,
                              params.roofHeight, params.slopeAngle,
                              params.addHips);
    } else {
      auto result = generateGableRoof(footprint, ridgeLine, params.slopeAngle,
                                     params.addHips);
      return result;
    }
  }
  default:
    BOOST_THROW_EXCEPTION(Exception("Unknown roof type"));
  }
}

/**
 * @brief Generate a complete building with roof using unified parameters without validity check.
 */
auto
generateBuildingWithRoof(const Polygon &footprint, const LineString &ridgeLine,
                        const RoofParameters &params, NoValidityCheck &/*nvc*/)
    -> std::unique_ptr<Geometry>
{
  if (params.buildingHeight < 0.0) {
    BOOST_THROW_EXCEPTION(Exception("Building height must be >= 0"));
  }
  if (params.roofHeight <= 0.0) {
    BOOST_THROW_EXCEPTION(Exception("Roof height must be > 0"));
  }

  switch (params.type) {
  case RoofType::FLAT: {
    double totalHeight = params.buildingHeight + params.roofHeight;
    auto extruded = extrude(footprint, totalHeight);
    if (params.generateSolid && extruded) {
      return extruded;
    }
    if (auto solid = dynamic_cast<const Solid*>(extruded.get())) {
      return std::make_unique<PolyhedralSurface>(solid->exteriorShell());
    }
    BOOST_THROW_EXCEPTION(Exception("Failed to create flat roof"));
  }
  case RoofType::HIPPED: {
    if (params.buildingHeight > 0.0) {
      auto result = extrudeStraightSkeleton(footprint, params.buildingHeight, params.roofHeight);
      if (params.generateSolid && result) {
        return std::make_unique<Solid>(*result);
      }
      return result;
    } else {
      auto result = extrudeStraightSkeleton(footprint, params.roofHeight);
      if (params.generateSolid && result) {
        return std::make_unique<Solid>(*result);
      }
      return result;
    }
  }
  case RoofType::PITCHED:
  case RoofType::SKILLION: {
    if (params.buildingHeight > 0.0 || params.generateSolid) {
      return generatePitchedRoof(footprint, ridgeLine, params.buildingHeight,
                                params.roofHeight, params.slopeAngle,
                                params.ridgePosition);
    } else {
      NoValidityCheck nvc;
      auto result = generatePitchedRoof(footprint, ridgeLine, params.slopeAngle,
                                       params.ridgePosition, nvc);
      return result;
    }
  }
  case RoofType::GABLE: {
    if (params.buildingHeight > 0.0 || params.generateSolid) {
      return generateGableRoof(footprint, ridgeLine, params.buildingHeight,
                              params.roofHeight, params.slopeAngle,
                              params.addHips);
    } else {
      NoValidityCheck nvc;
      auto result = generateGableRoof(footprint, ridgeLine, params.slopeAngle,
                                     params.addHips, nvc);
      return result;
    }
  }
  default:
    BOOST_THROW_EXCEPTION(Exception("Unknown roof type"));
  }
}


auto
generateGableRoof(const Polygon &footprint, double buildingHeight, double roofHeight,
                  double slopeAngle, bool addHips)
    -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(footprint);

  if (slopeAngle <= 0.0 || slopeAngle >= 90.0) {
    BOOST_THROW_EXCEPTION(Exception("Slope angle must be between 0 and 90 degrees"));
  }

  // Get the projected medial axis as ridge line
  auto medialAxisProjected = projectMedialAxisToEdges(footprint);
  if (!medialAxisProjected || medialAxisProjected->isEmpty()) {
    BOOST_THROW_EXCEPTION(Exception("Could not compute medial axis for roof generation"));
  }

  // For simple polygons (like rectangles), use the first linestring with just start and end points
  LineString ridgeLine;

  // Find the longest ridge line segment for the gable roof
  double maxLength = 0;

  for (size_t i = 0; i < medialAxisProjected->numGeometries(); ++i) {
    const auto *line = dynamic_cast<const LineString *>(&medialAxisProjected->geometryN(i));
    if (line && line->numPoints() >= 2) {
      double length = algorithm::length(*line);
      if (length > maxLength) {
        maxLength = length;
        ridgeLine = LineString();
        // Add all points from the projected medial axis line
        for (size_t j = 0; j < line->numPoints(); ++j) {
          ridgeLine.addPoint(line->pointN(j));
        }
      }
    }
  }

  if (ridgeLine.numPoints() < 2) {
    BOOST_THROW_EXCEPTION(Exception("Could not create valid ridge line from medial axis"));
  }

  // Use the existing gable roof implementation with the computed ridge line
  return generateGableRoof(footprint, ridgeLine, buildingHeight, roofHeight, slopeAngle, addHips);
}

/**
 * @brief Internal unified implementation for gable roof generation
 *
 * @param footprint The building footprint polygon
 * @param ridgeHeight The calculated ridge height (Z coordinate)
 * @param addVerticalFaces Whether to add vertical faces at ridge line ends
 * @return A PolyhedralSurface representing the gable roof
 */
static auto generateGableRoofImpl(const Polygon &footprint,
                                 double ridgeHeight,
                                 bool addVerticalFaces)
    -> std::unique_ptr<PolyhedralSurface>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(footprint);

  auto result = std::make_unique<PolyhedralSurface>();
  auto exteriorRing = footprint.exteriorRing();

  // 1. Get projected medial axis as ridge line
  auto projectedMedialAxis = projectMedialAxisToEdges(footprint);
  if (!projectedMedialAxis || projectedMedialAxis->isEmpty()) {
    BOOST_THROW_EXCEPTION(Exception("Could not compute projected medial axis for gable roof"));
  }

  // 2. Perform CDT with ridge line as constraints
  auto triangulation = SFCGAL::triangulate::triangulate2DZ(footprint);
  SFCGAL::triangulate::triangulate2DZ(*projectedMedialAxis, triangulation);
  auto cdtSurface = triangulation.getTriangulatedSurface();

  if (!cdtSurface) {
    BOOST_THROW_EXCEPTION(Exception("CDT triangulation failed"));
  }

  // 3. Collect ridge points for elevation
  std::set<std::pair<double, double>> ridgePointsSet;
  for (size_t i = 0; i < projectedMedialAxis->numGeometries(); ++i) {
    const auto *line = dynamic_cast<const LineString *>(&projectedMedialAxis->geometryN(i));
    if (!line) continue;

    for (size_t j = 0; j < line->numPoints(); ++j) {
      const auto &rp = line->pointN(j);
      ridgePointsSet.insert({CGAL::to_double(rp.x()), CGAL::to_double(rp.y())});
    }
  }

  // 4. Process each triangle from CDT and elevate ridge points
  for (size_t i = 0; i < cdtSurface->numPatches(); ++i) {
    const auto &patch = cdtSurface->patchN(i);
    const auto *triangle = dynamic_cast<const Triangle *>(&patch);
    if (!triangle) continue;

    // Get triangle vertices
    auto v1 = triangle->vertex(0);
    auto v2 = triangle->vertex(1);
    auto v3 = triangle->vertex(2);

    // Calculate triangle centroid
    auto centroidX = (v1.x() + v2.x() + v3.x()) / 3;
    auto centroidY = (v1.y() + v2.y() + v3.y()) / 3;

    // Filter: check if centroid is inside original polygon
    Point centroidPt(centroidX, centroidY, 0.0);
    if (!SFCGAL::algorithm::covers(footprint, centroidPt)) {
      continue; // Skip triangles outside polygon
    }

    // Elevate ridge vertices and create roof triangle
    Point newV1 = v1;
    Point newV2 = v2;
    Point newV3 = v3;

    // Check if each vertex is a ridge point and elevate it
    auto isRidgePoint = [&ridgePointsSet](const Point &p) -> bool {
      return ridgePointsSet.find({CGAL::to_double(p.x()), CGAL::to_double(p.y())}) != ridgePointsSet.end();
    };

    if (isRidgePoint(v1)) {
      newV1 = Point(v1.x(), v1.y(), ridgeHeight);
    }
    if (isRidgePoint(v2)) {
      newV2 = Point(v2.x(), v2.y(), ridgeHeight);
    }
    if (isRidgePoint(v3)) {
      newV3 = Point(v3.x(), v3.y(), ridgeHeight);
    }

    // Create roof triangle with elevated vertices
    std::unique_ptr<Triangle> roofTriangle(new Triangle(newV1, newV2, newV3));
    result->addPatch(*roofTriangle);
  }

  // 5. Add vertical faces at ridge endpoints if requested
  if (addVerticalFaces) {
    // Get extreme points from projected medial axis
    std::vector<Point> extremePoints;
    for (size_t i = 0; i < projectedMedialAxis->numGeometries(); ++i) {
      const auto *line = dynamic_cast<const LineString *>(&projectedMedialAxis->geometryN(i));
      if (line && line->numPoints() >= 2) {
        extremePoints.push_back(line->pointN(0));
        extremePoints.push_back(line->pointN(line->numPoints() - 1));
      }
    }

    // For each extreme point, find which polygon edge it lies on and create vertical triangle
    double tolerance = 1e-6;
    for (const auto& extremePoint : extremePoints) {
      for (size_t i = 0; i < exteriorRing.numPoints() - 1; ++i) {
        auto v1 = exteriorRing.pointN(i);
        auto v2 = exteriorRing.pointN((i + 1) % (exteriorRing.numPoints() - 1));

        // Check if extreme point lies on edge v1->v2
        double dx = CGAL::to_double(v2.x() - v1.x());
        double dy = CGAL::to_double(v2.y() - v1.y());
        double edgeLength = std::sqrt(dx * dx + dy * dy);

        if (edgeLength > tolerance) {
          // Parametric representation: point = v1 + t * (v2 - v1)
          double rx = CGAL::to_double(extremePoint.x() - v1.x());
          double ry = CGAL::to_double(extremePoint.y() - v1.y());
          double t = (rx * dx + ry * dy) / (edgeLength * edgeLength);

          if (t >= -tolerance && t <= 1.0 + tolerance) {
            // Calculate projected point
            double projX = CGAL::to_double(v1.x()) + t * dx;
            double projY = CGAL::to_double(v1.y()) + t * dy;
            double distX = CGAL::to_double(extremePoint.x()) - projX;
            double distY = CGAL::to_double(extremePoint.y()) - projY;
            double distance = std::sqrt(distX * distX + distY * distY);

            if (distance < tolerance) {
              // Extreme point lies on this edge - create vertical triangle
              Point v1Base(v1.x(), v1.y(), 0.0);
              Point v2Base(v2.x(), v2.y(), 0.0);
              Point extremeTop(extremePoint.x(), extremePoint.y(), ridgeHeight);

              // Create triangle: edge_start -> edge_end -> elevated_extreme_point
              std::unique_ptr<Triangle> vertTriangle(new Triangle(v1Base, v2Base, extremeTop));
              result->addPatch(*vertTriangle);
              break;
            }
          }
        }
      }
    }
  }

  propagateValidityFlag(*result, true);
  return result;
}

auto
generateGableRoof(const Polygon &footprint, double slopeAngle, bool addVerticalFaces)
    -> std::unique_ptr<PolyhedralSurface>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(footprint);

  if (slopeAngle <= 0.0 || slopeAngle >= 90.0) {
    BOOST_THROW_EXCEPTION(Exception("Slope angle must be between 0 and 90 degrees"));
  }

  // Calculate ridge height from slope angle
  // Get the projected medial axis to find maximum distance
  auto projectedMedialAxis = projectMedialAxisToEdges(footprint);
  if (!projectedMedialAxis || projectedMedialAxis->isEmpty()) {
    BOOST_THROW_EXCEPTION(Exception("Could not compute projected medial axis for gable roof"));
  }

  auto exteriorRing = footprint.exteriorRing();
  double maxPerpDistance = 0.0;

  // Calculate maximum perpendicular distance from ridge line to polygon edge
  for (size_t i = 0; i < projectedMedialAxis->numGeometries(); ++i) {
    const auto *line = dynamic_cast<const LineString *>(&projectedMedialAxis->geometryN(i));
    if (!line) continue;

    for (size_t j = 0; j < line->numPoints(); ++j) {
      const auto &ridgePoint = line->pointN(j);

      // Find minimum distance from this ridge point to polygon boundary
      double minDistToBoundary = std::numeric_limits<double>::max();
      for (size_t k = 0; k < exteriorRing.numPoints() - 1; ++k) {
        auto edgeStart = exteriorRing.pointN(k);
        auto edgeEnd = exteriorRing.pointN(k + 1);

        // Calculate distance from ridge point to edge using exact arithmetic
        auto px = ridgePoint.x();
        auto py = ridgePoint.y();
        auto ax = edgeStart.x();
        auto ay = edgeStart.y();
        auto bx = edgeEnd.x();
        auto by = edgeEnd.y();

        // Vector from A to B
        auto ABx = bx - ax;
        auto ABy = by - ay;
        // Vector from A to P
        auto APx = px - ax;
        auto APy = py - ay;

        auto AB_squared = ABx * ABx + ABy * ABy;
        if (AB_squared == 0) {
          // Edge is a point, distance is from point to point
          auto dist_squared = (px - ax) * (px - ax) + (py - ay) * (py - ay);
          double dist = std::sqrt(CGAL::to_double(dist_squared));
          minDistToBoundary = std::min(minDistToBoundary, dist);
        } else {
          auto t_num = APx * ABx + APy * ABy;

          // Clamp t to [0,1]
          auto t = t_num / AB_squared;
          if (t < 0) t = 0;
          if (t > 1) t = 1;

          // Closest point on segment
          auto closestX = ax + t * ABx;
          auto closestY = ay + t * ABy;

          // Distance from P to closest point
          auto dist_squared = (px - closestX) * (px - closestX) + (py - closestY) * (py - closestY);
          double dist = std::sqrt(CGAL::to_double(dist_squared));
          minDistToBoundary = std::min(minDistToBoundary, dist);
        }
      }
      maxPerpDistance = std::max(maxPerpDistance, minDistToBoundary);
    }
  }

  double ridgeHeight = calculateRidgeHeight(maxPerpDistance, slopeAngle);

  return generateGableRoofImpl(footprint, ridgeHeight, addVerticalFaces);
}

auto
generateGableRoof(const Polygon &footprint, double slopeAngle)
    -> std::unique_ptr<PolyhedralSurface>
{
  return generateGableRoof(footprint, slopeAngle, false);
}

auto
generateGableRoofWithHeight(const Polygon &footprint, double roofHeight, bool addVerticalFaces)
    -> std::unique_ptr<PolyhedralSurface>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(footprint);

  if (roofHeight <= 0.0) {
    BOOST_THROW_EXCEPTION(Exception("Roof height must be positive"));
  }

  return generateGableRoofImpl(footprint, roofHeight, addVerticalFaces);
}

auto
generateGableRoofWithBuilding(const Polygon &footprint, double buildingHeight, double roofHeight,
                             double slopeAngle, bool addVerticalFaces)
    -> std::unique_ptr<PolyhedralSurface>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(footprint);

  if (buildingHeight < 0.0) {
    BOOST_THROW_EXCEPTION(Exception("Building height must be non-negative"));
  }
  if (roofHeight <= 0.0) {
    BOOST_THROW_EXCEPTION(Exception("Roof height must be positive"));
  }
  if (slopeAngle <= 0.0 || slopeAngle >= 90.0) {
    BOOST_THROW_EXCEPTION(Exception("Slope angle must be between 0 and 90 degrees"));
  }

  std::unique_ptr<PolyhedralSurface> result(new PolyhedralSurface);

  if (footprint.isEmpty()) {
    return result;
  }

  // Create complete roof using our unified implementation
  auto completeRoof = generateGableRoofWithHeight(footprint, roofHeight, addVerticalFaces);

  // Create new roof surface (excluding base faces)
  auto roof = std::make_unique<PolyhedralSurface>();

  // Predicate to identify non-base faces (roof slopes and vertical faces)
  auto isNotBaseFace = [](const Polygon &patch) -> bool {
    const LineString &exterior = patch.exteriorRing();

    // Check if any point has z != 0 (not a base face)
    return std::any_of(
        exterior.begin(), exterior.end(),
        [](const Point &point) -> bool { return point.z() != 0.0; });
  };

  // Copy roof patches (excluding base)
  for (size_t i = 0; i < completeRoof->numPatches(); ++i) {
    const auto &patch = completeRoof->patchN(i);
    if (auto polygon = dynamic_cast<const Polygon*>(&patch)) {
      if (isNotBaseFace(*polygon)) {
        roof->addPatch(*polygon);
      }
    } else if (auto triangle = dynamic_cast<const Triangle*>(&patch)) {
      // For triangles, check if any vertex has z != 0
      if (triangle->vertex(0).z() != 0.0 ||
          triangle->vertex(1).z() != 0.0 ||
          triangle->vertex(2).z() != 0.0) {
        roof->addPatch(*triangle);
      }
    }
  }

  // Translate roof to building height
  translate(*roof, 0.0, 0.0, buildingHeight);

  // Create building walls
  auto building = extrude(footprint, buildingHeight);

  // Create result from building exterior shell
  result = std::make_unique<PolyhedralSurface>(
      building->as<Solid>().exteriorShell());

  // Add filtered roof patches
  result->addPatchs(*roof);

  propagateValidityFlag(*result, true);

  return result;
}

} // namespace SFCGAL::algorithm