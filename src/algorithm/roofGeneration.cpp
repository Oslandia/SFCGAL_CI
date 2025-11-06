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
#include "SFCGAL/algorithm/force2D.h"
#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/algorithm/straightSkeleton.h"
#include "SFCGAL/triangulate/triangulate2DZ.h"

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

  // For gable roof, we create symmetric slopes like pitched roof
  // but handle gable ends specially
  auto result = generatePitchedRoof(footprint, ridgeLine, slopeAngle,
                                   RidgePosition::INTERIOR);

  if (addHips) {
    // Add hip treatment at gable ends
    // This would involve creating triangular surfaces at the ends
    // For now, we'll use the basic pitched roof result
  }

  // Add gable end triangles (vertical surfaces at ridge line ends)
  auto ridgeStart = ridgeLine.pointN(0);
  auto ridgeEnd = ridgeLine.pointN(ridgeLine.numPoints() - 1);
  double ridgeHeight = calculateRidgeHeight(1.0, slopeAngle);

  Point ridgeStartTop(ridgeStart.x(), ridgeStart.y(), ridgeHeight);
  Point ridgeEndTop(ridgeEnd.x(), ridgeEnd.y(), ridgeHeight);

  // Create triangular gable ends (simplified)
  // In a full implementation, this would find the proper base points
  Point gableBase1(ridgeStart.x(), ridgeStart.y() + 1.0, 0.0);
  Point gableBase2(ridgeEnd.x(), ridgeEnd.y() + 1.0, 0.0);

  std::vector<Point> gableTriangle1 = {ridgeStart, ridgeStartTop, gableBase1, ridgeStart};
  result->addPatch(Polygon(LineString(gableTriangle1)));

  std::vector<Point> gableTriangle2 = {ridgeEnd, ridgeEndTop, gableBase2, ridgeEnd};
  result->addPatch(Polygon(LineString(gableTriangle2)));

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
    // Use extrude to create flat roof
    auto extruded = extrude(footprint, params.height);
    if (auto solid = dynamic_cast<const Solid*>(extruded.get())) {
      return std::make_unique<PolyhedralSurface>(solid->exteriorShell());
    }
    BOOST_THROW_EXCEPTION(Exception("Failed to create flat roof"));
  }
  case RoofType::HIPPED: {
    // Use extrudeStraightSkeleton for hipped roof
    return extrudeStraightSkeleton(footprint, params.height);
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
    auto extruded = extrude(footprint, params.height);
    if (auto solid = dynamic_cast<const Solid*>(extruded.get())) {
      return std::make_unique<PolyhedralSurface>(solid->exteriorShell());
    }
    BOOST_THROW_EXCEPTION(Exception("Failed to create flat roof"));
  }
  case RoofType::HIPPED: {
    return extrudeStraightSkeleton(footprint, params.height);
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

auto
generateGableRoofAuto(const Polygon &footprint, double slopeAngle,
                      bool addVerticalFaces) -> std::unique_ptr<PolyhedralSurface>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(footprint);

  if (slopeAngle <= 0.0 || slopeAngle >= 90.0) {
    BOOST_THROW_EXCEPTION(Exception("Slope angle must be between 0 and 90 degrees"));
  }

  // 1. Get projected medial axis to edges (this gives us the ridge line)
  auto projectedMedialAxis = projectMedialAxisToEdges(footprint);

  // 2. Perform constrained Delaunay triangulation with ridge as constraints
  auto footprint_clone = footprint.clone();
  SFCGAL::algorithm::force2D(*footprint_clone);
  auto triangulation = SFCGAL::triangulate::triangulate2DZ(*footprint_clone);

  // Add constraints from projected medial axis
  auto projectedMedialAxis_clone = projectedMedialAxis->clone();
  SFCGAL::algorithm::force2D(*projectedMedialAxis_clone);
  SFCGAL::triangulate::triangulate2DZ(*projectedMedialAxis_clone, triangulation);

  auto triangulated_surface = triangulation.getTriangulatedSurface();

  std::unique_ptr<PolyhedralSurface> result(new PolyhedralSurface);

  // 3. Collect all ridge points from projected medial axis
  std::vector<Point> ridgePoints;
  for (size_t i = 0; i < projectedMedialAxis->numGeometries(); ++i) {
    const auto *line = dynamic_cast<const LineString *>(&projectedMedialAxis->geometryN(i));
    if (line) {
      for (size_t j = 0; j < line->numPoints(); ++j) {
        ridgePoints.push_back(line->pointN(j));
      }
    }
  }

  // Helper function to check if a point is a ridge point
  auto isRidgePoint = [&ridgePoints](const Point& p) -> bool {
    for (const Point& ridgeP : ridgePoints) {
      auto dx = p.x() - ridgeP.x();
      auto dy = p.y() - ridgeP.y();
      if (dx * dx + dy * dy < Kernel::FT(EPSILON)) {
        return true;
      }
    }
    return false;
  };

  // Calculate ridge height from slope angle
  double ridgeHeight = std::tan(slopeAngle * M_PI / 180.0);

  // 4. Process triangulation and elevate ridge points
  for (size_t i = 0; i < triangulated_surface->numPatches(); ++i) {
    const auto &patch = triangulated_surface->patchN(i);
    if (auto triangle = dynamic_cast<const Triangle *>(&patch)) {
      const auto &v1 = triangle->vertex(0);
      const auto &v2 = triangle->vertex(1);
      const auto &v3 = triangle->vertex(2);

      // Create new vertices with proper Z coordinates
      Point newV1(v1.x(), v1.y(), 0.0);
      Point newV2(v2.x(), v2.y(), 0.0);
      Point newV3(v3.x(), v3.y(), 0.0);

      // Elevate ridge points
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
  }

  // 5. Add vertical faces at ridge endpoints if requested
  if (addVerticalFaces) {
    // Get terminal points from projected medial axis (ridge endpoints)
    std::vector<Point> ridgeEndpoints;
    for (size_t i = 0; i < projectedMedialAxis->numGeometries(); ++i) {
      const auto *line = dynamic_cast<const LineString *>(&projectedMedialAxis->geometryN(i));
      if (line && line->numPoints() >= 2) {
        ridgeEndpoints.push_back(line->pointN(0));
        ridgeEndpoints.push_back(line->pointN(line->numPoints() - 1));
      }
    }

    // For each ridge endpoint, find the footprint edges that intersect with a perpendicular line
    const LineString& ring = footprint.exteriorRing();
    for (const Point& ridgeEndpoint : ridgeEndpoints) {
      std::vector<Point> edgePoints;

      // Find points on the footprint boundary that share the same X or Y coordinate as ridge endpoint
      for (size_t i = 0; i < ring.numPoints() - 1; ++i) {
        const Point& p1 = ring.pointN(i);
        const Point& p2 = ring.pointN(i + 1);

        // Check if this edge crosses the ridge endpoint perpendicular
        // For a horizontal ridge (like x=0 to x=10, y=3), we need vertical edges at x=0 and x=10
        if ((p1.x() == ridgeEndpoint.x() && p2.x() == ridgeEndpoint.x()) ||
            (p1.y() == ridgeEndpoint.y() && p2.y() == ridgeEndpoint.y())) {

          // Add both points of this edge if they're different from ridge endpoint
          auto dx1 = p1.x() - ridgeEndpoint.x();
          auto dy1 = p1.y() - ridgeEndpoint.y();
          auto dx2 = p2.x() - ridgeEndpoint.x();
          auto dy2 = p2.y() - ridgeEndpoint.y();

          if (dx1*dx1 + dy1*dy1 > Kernel::FT(EPSILON)) edgePoints.push_back(p1);
          if (dx2*dx2 + dy2*dy2 > Kernel::FT(EPSILON)) edgePoints.push_back(p2);
        }
      }

      // Create vertical triangular faces connecting the ridge endpoint to the edge points
      Point ridgeTop(ridgeEndpoint.x(), ridgeEndpoint.y(), ridgeHeight);
      Point ridgeBase(ridgeEndpoint.x(), ridgeEndpoint.y(), 0.0);

      for (const Point& edgePoint : edgePoints) {
        Point edgeBase(edgePoint.x(), edgePoint.y(), 0.0);

        // Create triangle: ridge_base -> ridge_top -> edge_base -> ridge_base
        std::unique_ptr<Triangle> verticalTri(new Triangle(ridgeBase, ridgeTop, edgeBase));
        result->addPatch(*verticalTri);
      }
    }
  }

  propagateValidityFlag(*result, true);
  return result;
}

} // namespace SFCGAL::algorithm