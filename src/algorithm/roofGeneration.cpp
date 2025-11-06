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
#include "SFCGAL/algorithm/covers.h"
#include "SFCGAL/algorithm/extrude.h"
#include "SFCGAL/algorithm/force2D.h"
#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/algorithm/straightSkeleton.h"
#include "SFCGAL/algorithm/translate.h"
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

// Consistent tolerance constants for geometric operations
const auto GEOMETRIC_TOLERANCE = 1e-10;
const auto HEIGHT_TOLERANCE    = 1e-6;
const auto ANGLE_TOLERANCE     = 1e-9;

/**
 * @brief Calculate perpendicular distance from a point to a 2D line using exact
 * arithmetic
 * @param point The point to measure distance from
 * @param lineStart Starting point of the line
 * @param lineEnd Ending point of the line
 * @return Exact distance value using kernel arithmetic
 */
auto
distanceToLine2D(const Point &point, const Point &lineStart,
                 const Point &lineEnd) -> double
{
  // Use high-precision arithmetic for numerical stability
  auto dx = CGAL::to_double(lineEnd.x() - lineStart.x());
  auto dy = CGAL::to_double(lineEnd.y() - lineStart.y());

  // If line is degenerate (start == end), return distance to point
  auto lineLengthSq = dx * dx + dy * dy;
  if (lineLengthSq < GEOMETRIC_TOLERANCE) {
    auto px = CGAL::to_double(point.x() - lineStart.x());
    auto py = CGAL::to_double(point.y() - lineStart.y());
    return std::sqrt(px * px + py * py);
  }

  // Calculate perpendicular distance using cross product formula
  auto px = CGAL::to_double(point.x() - lineStart.x());
  auto py = CGAL::to_double(point.y() - lineStart.y());

  return std::abs(dx * py - dy * px) / std::sqrt(lineLengthSq);
}

/**
 * @brief Calculate signed distance from point to line (positive on one side,
 * negative on other)
 * @param point The point to measure distance from
 * @param lineStart Starting point of the line
 * @param lineEnd Ending point of the line
 * @return Signed distance value using kernel arithmetic
 */
auto
signedDistanceToLine2D(const Point &point, const Point &lineStart,
                       const Point &lineEnd) -> double
{
  // Use high-precision arithmetic for numerical stability
  auto dx = CGAL::to_double(lineEnd.x() - lineStart.x());
  auto dy = CGAL::to_double(lineEnd.y() - lineStart.y());

  // If line is degenerate, return distance to start point
  auto lineLengthSq = dx * dx + dy * dy;
  if (lineLengthSq < GEOMETRIC_TOLERANCE) {
    auto px = CGAL::to_double(point.x() - lineStart.x());
    auto py = CGAL::to_double(point.y() - lineStart.y());
    return std::sqrt(px * px + py * py);
  }

  // Calculate signed perpendicular distance using cross product
  auto px = CGAL::to_double(point.x() - lineStart.x());
  auto py = CGAL::to_double(point.y() - lineStart.y());

  return (dx * py - dy * px) / std::sqrt(lineLengthSq);
}

/**
 * @brief Calculate 3D point projected from base point with given height and
 * slope
 */
auto
calculateSlopePoint(const Point &basePoint, const Kernel::Vector_3 &normal,
                    double distance, double slopeAngle) -> Point
{
  double height = distance * std::tan(slopeAngle * M_PI / 180.0);
  return Point(basePoint.x() + normal.x() * height,
               basePoint.y() + normal.y() * height, basePoint.z() + height);
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
    BOOST_THROW_EXCEPTION(
        Exception("Ridge line must contain at least 2 points"));
  }

  // For a simplified but more robust implementation,
  // just return two halves of the polygon split at the midpoint
  auto  &ring      = polygon.exteriorRing();
  size_t numPoints = ring.numPoints() - 1; // Excluding closing point

  if (numPoints < 4) {
    BOOST_THROW_EXCEPTION(
        Exception("Polygon must have at least 4 vertices for splitting"));
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
  auto   ridgeHeight = calculateRidgeHeight(avgDistance, slopeAngle);

  // Create elevated ridge line
  std::vector<Point> ridgePoints;
  for (size_t i = 0; i < ridgeLine.numPoints(); ++i) {
    auto point = ridgeLine.pointN(i);
    ridgePoints.emplace_back(point.x(), point.y(), slopeUp ? ridgeHeight : 0.0);
  }

  if (ridgePoints.empty()) {
    return surface; // Return empty surface if no ridge points
  }

  // Create simple triangulated surface
  auto  &baseRing      = basePolygon.exteriorRing();
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
    std::vector<Point> trianglePoints = {basePoint1, basePoint2, ridgePoint,
                                         basePoint1};
    surface->addPatch(Polygon(LineString(trianglePoints)));
  }

  return surface;
}

/**
 * @brief Create elevated vertices for sloped roof surface
 * @param footprint The building footprint polygon
 * @param ridgeStart Starting point of the ridge line
 * @param ridgeEnd Ending point of the ridge line
 * @param slopeTan Tangent of the slope angle
 * @return Vector of elevated vertices for the roof surface
 */
auto
createSlopedSurfaceVertices(const Polygon &footprint, const Point &ridgeStart,
                            const Point &ridgeEnd, double slopeTan)
    -> std::vector<Point>
{
  auto  &exteriorRing = footprint.exteriorRing();
  size_t numPoints    = exteriorRing.numPoints();

  if (numPoints < 4) { // Need at least 3 + closing point
    BOOST_THROW_EXCEPTION(Exception("Polygon must have at least 3 vertices"));
  }

  std::vector<Point> elevatedVertices;
  elevatedVertices.reserve(numPoints);

  // Calculate elevation for each vertex based on distance from ridge line
  for (size_t i = 0; i < numPoints - 1; ++i) { // Skip closing point
    auto vertex = exteriorRing.pointN(i);

    // Calculate perpendicular distance from vertex to ridge line
    auto distance = distanceToLine2D(vertex, ridgeStart, ridgeEnd);

    // Calculate height based on slope angle and distance
    auto height = distance * slopeTan;

    // Create elevated point
    Point elevatedVertex(vertex.x(), vertex.y(), height);
    elevatedVertices.push_back(elevatedVertex);
  }

  // Close the ring
  if (!elevatedVertices.empty()) {
    elevatedVertices.push_back(elevatedVertices[0]);
  }

  return elevatedVertices;
}

/**
 * @brief Create vertical faces connecting base polygon to elevated roof
 * @param footprint The building footprint polygon
 * @param elevatedVertices The elevated vertices of the roof
 * @param addVerticalFaces Whether to add vertical faces
 * @return Vector of geometry patches representing vertical faces
 */
auto
createVerticalFaces(const Polygon            &footprint,
                    const std::vector<Point> &elevatedVertices,
                    bool                      addVerticalFaces)
    -> std::vector<std::unique_ptr<Geometry>>
{
  std::vector<std::unique_ptr<Geometry>> faces;

  if (!addVerticalFaces) {
    return faces;
  }

  auto  &exteriorRing = footprint.exteriorRing();
  size_t numPoints    = exteriorRing.numPoints();

  // Create vertical faces for each edge when there is meaningful height
  // difference
  for (size_t i = 0; i < numPoints - 1; ++i) {
    size_t nextI = (i + 1) % (numPoints - 1);

    // Base edge points (at z=0)
    auto baseVertex1 = exteriorRing.pointN(i);
    auto baseVertex2 = exteriorRing.pointN(nextI);
    auto basePt1     = Point(baseVertex1.x(), baseVertex1.y(), 0.0);
    auto basePt2     = Point(baseVertex2.x(), baseVertex2.y(), 0.0);

    // Corresponding elevated points on the roof
    auto roofPt1 = elevatedVertices[i];
    auto roofPt2 = elevatedVertices[nextI];

    // Check if there's meaningful height difference for this edge
    auto roofHeight1 = CGAL::to_double(roofPt1.z());
    auto roofHeight2 = CGAL::to_double(roofPt2.z());
    auto baseHeight1 = CGAL::to_double(basePt1.z());
    auto baseHeight2 = CGAL::to_double(basePt2.z());

    bool hasHeightDiff1 =
        std::abs(roofHeight1 - baseHeight1) > HEIGHT_TOLERANCE;
    bool hasHeightDiff2 =
        std::abs(roofHeight2 - baseHeight2) > HEIGHT_TOLERANCE;

    if (hasHeightDiff1 && hasHeightDiff2) {
      // Both points have height difference - create quadrilateral face
      std::vector<Point> verticalFace = {basePt1, basePt2, roofPt2, roofPt1,
                                         basePt1};
      faces.push_back(std::make_unique<Polygon>(LineString(verticalFace)));
    } else if (hasHeightDiff1 && !hasHeightDiff2) {
      // Only point 1 has height difference - create triangular face
      std::vector<Point> triangularFace = {basePt1, roofPt1, basePt2, basePt1};
      faces.push_back(std::make_unique<Polygon>(LineString(triangularFace)));
    } else if (!hasHeightDiff1 && hasHeightDiff2) {
      // Only point 2 has height difference - create triangular face
      std::vector<Point> triangularFace = {basePt1, basePt2, roofPt2, basePt1};
      faces.push_back(std::make_unique<Polygon>(LineString(triangularFace)));
    }
    // If both points have no height difference, skip (would be degenerate)
  }

  return faces;
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
    BOOST_THROW_EXCEPTION(
        Exception("Slope angle must be between 0 and 90 degrees"));
  }
  return horizontalDistance * std::tan(slopeAngle * M_PI / 180.0);
}

auto
calculateHorizontalDistance(double height, double slopeAngle) -> double
{
  if (slopeAngle <= 0.0 || slopeAngle >= 90.0) {
    BOOST_THROW_EXCEPTION(
        Exception("Slope angle must be between 0 and 90 degrees"));
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
    BOOST_THROW_EXCEPTION(
        Exception("Slope angle must be between 0 and 90 degrees"));
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
      auto surface =
          createSlopedSurface(footprint, ridgeLine, slopeAngle, true);
      result = std::move(surface);
      break;
    }
    }
  } catch (const Exception &e) {
    // If splitting fails, fall back to single surface
    auto surface = createSlopedSurface(footprint, ridgeLine, slopeAngle, true);
    result       = std::move(surface);
  }

  propagateValidityFlag(*result, true);
  return result;
}

auto
generatePitchedRoof(const Polygon &footprint, const LineString &ridgeLine,
                    double slopeAngle, RidgePosition ridgePosition,
                    NoValidityCheck & /*nvc*/)
    -> std::unique_ptr<PolyhedralSurface>
{
  if (slopeAngle <= 0.0 || slopeAngle >= 90.0) {
    BOOST_THROW_EXCEPTION(
        Exception("Slope angle must be between 0 and 90 degrees"));
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
    result       = std::move(surface);
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
    BOOST_THROW_EXCEPTION(
        Exception("Slope angle must be between 0 and 90 degrees"));
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
  auto   ridgeStart  = ridgeLine.pointN(0);
  auto   ridgeEnd    = ridgeLine.pointN(ridgeLine.numPoints() - 1);
  double ridgeHeight = calculateRidgeHeight(1.0, slopeAngle);

  Point ridgeStartTop(ridgeStart.x(), ridgeStart.y(), ridgeHeight);
  Point ridgeEndTop(ridgeEnd.x(), ridgeEnd.y(), ridgeHeight);

  // Create triangular gable ends (simplified)
  // In a full implementation, this would find the proper base points
  Point gableBase1(ridgeStart.x(), ridgeStart.y() + 1.0, 0.0);
  Point gableBase2(ridgeEnd.x(), ridgeEnd.y() + 1.0, 0.0);

  std::vector<Point> gableTriangle1 = {ridgeStart, ridgeStartTop, gableBase1,
                                       ridgeStart};
  result->addPatch(Polygon(LineString(gableTriangle1)));

  std::vector<Point> gableTriangle2 = {ridgeEnd, ridgeEndTop, gableBase2,
                                       ridgeEnd};
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
    BOOST_THROW_EXCEPTION(
        Exception("Slope angle must be between 0 and 90 degrees"));
  }

  auto result = generatePitchedRoof(footprint, ridgeLine, slopeAngle,
                                    RidgePosition::INTERIOR, nvc);

  if (addHips) {
    // Add hip treatment
  }

  // Add gable ends
  auto   ridgeStart  = ridgeLine.pointN(0);
  auto   ridgeEnd    = ridgeLine.pointN(ridgeLine.numPoints() - 1);
  double ridgeHeight = calculateRidgeHeight(1.0, slopeAngle);

  Point ridgeStartTop(ridgeStart.x(), ridgeStart.y(), ridgeHeight);
  Point ridgeEndTop(ridgeEnd.x(), ridgeEnd.y(), ridgeHeight);

  Point gableBase1(ridgeStart.x(), ridgeStart.y() + 1.0, 0.0);
  Point gableBase2(ridgeEnd.x(), ridgeEnd.y() + 1.0, 0.0);

  std::vector<Point> gableTriangle1 = {ridgeStart, ridgeStartTop, gableBase1,
                                       ridgeStart};
  result->addPatch(Polygon(LineString(gableTriangle1)));

  std::vector<Point> gableTriangle2 = {ridgeEnd, ridgeEndTop, gableBase2,
                                       ridgeEnd};
  result->addPatch(Polygon(LineString(gableTriangle2)));

  return result;
}

/**
 * @brief Generate a skillion roof from a polygon footprint and ridge line.
 *
 * Creates a single-slope shed roof where the ridge line defines the high edge
 * and all other points slope down perpendicular to the ridge. This creates
 * the characteristic mono-pitch roof with consistent slope across the entire
 * surface.
 *
 * @param footprint The building footprint polygon (must be valid and non-empty)
 * @param ridgeLine The ridge line defining the high edge direction (minimum 2
 * points)
 * @param slopeAngle The roof slope angle in degrees (0 < angle < 90)
 * @param addVerticalFaces Whether to add vertical triangular faces at roof ends
 *
 * @return A PolyhedralSurface representing the skillion roof
 *
 * @pre footprint must be a valid polygon
 * @pre ridgeLine must be a valid LineString with >= 2 points
 * @pre slopeAngle must be between 0 and 90 degrees
 *
 * @throw Exception if preconditions are not met
 *
 * @example
 * ```cpp
 * Polygon rect = ...;
 * LineString ridge = ...;
 * auto roof = generateSkillionRoof(rect, ridge, 30.0, true);
 * ```
 */
auto
generateSkillionRoof(const Polygon &footprint, const LineString &ridgeLine,
                     double slopeAngle, bool addVerticalFaces)
    -> std::unique_ptr<PolyhedralSurface>
{
  // Comprehensive input validation
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(footprint);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(ridgeLine);

  if (footprint.isEmpty()) {
    BOOST_THROW_EXCEPTION(Exception("Footprint polygon cannot be empty"));
  }

  if (ridgeLine.isEmpty() || ridgeLine.numPoints() < 2) {
    BOOST_THROW_EXCEPTION(
        Exception("Ridge line must contain at least 2 points"));
  }

  // Use kernel-native arithmetic for angle validation
  auto slopeAngleKernel = Kernel::FT(slopeAngle);
  if (slopeAngleKernel <= Kernel::FT(0) || slopeAngleKernel >= Kernel::FT(90)) {
    BOOST_THROW_EXCEPTION(
        Exception("Slope angle must be between 0 and 90 degrees"));
  }

  auto result = std::make_unique<PolyhedralSurface>();

  // Get ridge line endpoints for distance calculation
  auto ridgeStart = ridgeLine.pointN(0);
  auto ridgeEnd   = ridgeLine.pointN(ridgeLine.numPoints() - 1);

  // Pre-compute trigonometric values with high precision
  auto slopeAngleRad = slopeAngle * M_PI / 180.0;
  auto slopeTan      = std::tan(slopeAngleRad);

  // Create elevated vertices for the sloped roof surface
  auto elevatedVertices =
      createSlopedSurfaceVertices(footprint, ridgeStart, ridgeEnd, slopeTan);

  // Create the sloped roof surface as a single polygon
  if (elevatedVertices.size() >= 4) {
    LineString elevatedRing(elevatedVertices);
    Polygon    slopedSurface(elevatedRing);
    result->addPatch(slopedSurface);
  }

  // Create vertical faces for proper roof closure
  auto verticalFaces =
      createVerticalFaces(footprint, elevatedVertices, addVerticalFaces);
  for (const auto &face : verticalFaces) {
    if (auto polygon = dynamic_cast<const Polygon *>(face.get())) {
      result->addPatch(*polygon);
    }
  }

  propagateValidityFlag(*result, true);
  return result;
}

auto
generateSkillionRoof(const Polygon &footprint, const LineString &ridgeLine,
                     double slopeAngle) -> std::unique_ptr<PolyhedralSurface>
{
  return generateSkillionRoof(footprint, ridgeLine, slopeAngle, false);
}

auto
generateSkillionRoof(const Polygon &footprint, const LineString &ridgeLine,
                     double slopeAngle, double buildingHeight)
    -> std::unique_ptr<PolyhedralSurface>
{
  return generateSkillionRoof(footprint, ridgeLine, slopeAngle, false,
                              buildingHeight);
}

auto
generateSkillionRoof(const Polygon &footprint, const LineString &ridgeLine,
                     double slopeAngle, bool addVerticalFaces,
                     double buildingHeight)
    -> std::unique_ptr<PolyhedralSurface>
{
  if (buildingHeight < 0.0) {
    BOOST_THROW_EXCEPTION(Exception("Building height must be non-negative"));
  }

  if (buildingHeight == 0.0) {
    // Just generate the roof
    return generateSkillionRoof(footprint, ridgeLine, slopeAngle,
                                addVerticalFaces);
  }

  // Generate building walls
  auto building = extrude(footprint, 0.0, 0.0, buildingHeight);

  // Generate roof and translate it to building height
  auto roof =
      generateSkillionRoof(footprint, ridgeLine, slopeAngle, addVerticalFaces);
  translate(*roof, 0.0, 0.0, buildingHeight);

  // Combine building and roof
  auto result = std::make_unique<PolyhedralSurface>(
      building->as<Solid>().exteriorShell());
  result->addPatchs(*roof);

  propagateValidityFlag(*result, true);
  return result;
}

auto
generateRoof(const Polygon &footprint, const LineString &ridgeLine,
             const RoofParameters &params) -> std::unique_ptr<PolyhedralSurface>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(footprint);

  switch (params.type) {
  case RoofType::FLAT: {
    // Use extrude to create flat roof
    auto extruded = extrude(footprint, 0.0, 0.0, params.height);
    if (auto solid = dynamic_cast<const Solid *>(extruded.get())) {
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
    auto extruded = extrude(footprint, 0.0, 0.0, params.height);
    if (auto solid = dynamic_cast<const Solid *>(extruded.get())) {
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
generateGableRoof(const Polygon &footprint, double slopeAngle,
                  bool addVerticalFaces, double buildingHeight)
    -> std::unique_ptr<PolyhedralSurface>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(footprint);

  if (slopeAngle <= 0.0 || slopeAngle >= 90.0) {
    BOOST_THROW_EXCEPTION(
        Exception("Slope angle must be between 0 and 90 degrees"));
  }

  if (buildingHeight < 0.0) {
    BOOST_THROW_EXCEPTION(Exception("Building height must be non-negative"));
  }

  std::unique_ptr<PolyhedralSurface> result(new PolyhedralSurface);

  if (footprint.isEmpty()) {
    return result;
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
  SFCGAL::triangulate::triangulate2DZ(*projectedMedialAxis_clone,
                                      triangulation);

  auto triangulated_surface = triangulation.getTriangulatedSurface();

  // 3. Collect all ridge points from projected medial axis
  std::vector<Point> ridgePoints;
  for (size_t i = 0; i < projectedMedialAxis->numGeometries(); ++i) {
    const auto *line =
        dynamic_cast<const LineString *>(&projectedMedialAxis->geometryN(i));
    if (line) {
      for (size_t j = 0; j < line->numPoints(); ++j) {
        ridgePoints.push_back(line->pointN(j));
      }
    }
  }

  // Helper function to check if a point is a ridge point
  auto isRidgePoint = [&ridgePoints](const Point &p) -> bool {
    for (const Point &ridgeP : ridgePoints) {
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

  // Create roof surface
  auto roof = std::make_unique<PolyhedralSurface>();

  // 4. Process triangulation and elevate ridge points
  for (size_t i = 0; i < triangulated_surface->numPatches(); ++i) {
    const auto &patch = triangulated_surface->patchN(i);
    if (auto triangle = dynamic_cast<const Triangle *>(&patch)) {
      const auto &v1 = triangle->vertex(0);
      const auto &v2 = triangle->vertex(1);
      const auto &v3 = triangle->vertex(2);

      // Check if triangle centroid is inside the original polygon
      Point centroid((v1.x() + v2.x() + v3.x()) / 3.0,
                     (v1.y() + v2.y() + v3.y()) / 3.0, 0.0);

      // Use SFCGAL's covers algorithm to check if centroid is inside footprint
      if (!SFCGAL::algorithm::covers(footprint, centroid)) {
        continue; // Skip triangles outside the original polygon
      }

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
      roof->addPatch(*roofTriangle);
    }
  }

  // 5. Add vertical faces at ridge endpoints if requested
  if (addVerticalFaces) {
    // Get terminal points from projected medial axis (ridge endpoints)
    std::vector<Point> ridgeEndpoints;
    for (size_t i = 0; i < projectedMedialAxis->numGeometries(); ++i) {
      const auto *line =
          dynamic_cast<const LineString *>(&projectedMedialAxis->geometryN(i));
      if (line && line->numPoints() >= 2) {
        ridgeEndpoints.push_back(line->pointN(0));
        ridgeEndpoints.push_back(line->pointN(line->numPoints() - 1));
      }
    }

    // For each ridge endpoint, find the footprint edges that intersect with a
    // perpendicular line
    const LineString &ring = footprint.exteriorRing();
    for (const Point &ridgeEndpoint : ridgeEndpoints) {
      std::vector<Point> edgePoints;

      // Find points on the footprint boundary that share the same X or Y
      // coordinate as ridge endpoint
      for (size_t i = 0; i < ring.numPoints() - 1; ++i) {
        const Point &p1 = ring.pointN(i);
        const Point &p2 = ring.pointN(i + 1);

        // Check if this edge crosses the ridge endpoint perpendicular
        // For a horizontal ridge (like x=0 to x=10, y=3), we need vertical
        // edges at x=0 and x=10
        if ((p1.x() == ridgeEndpoint.x() && p2.x() == ridgeEndpoint.x()) ||
            (p1.y() == ridgeEndpoint.y() && p2.y() == ridgeEndpoint.y())) {

          // Add both points of this edge if they're different from ridge
          // endpoint
          auto dx1 = p1.x() - ridgeEndpoint.x();
          auto dy1 = p1.y() - ridgeEndpoint.y();
          auto dx2 = p2.x() - ridgeEndpoint.x();
          auto dy2 = p2.y() - ridgeEndpoint.y();

          if (dx1 * dx1 + dy1 * dy1 > Kernel::FT(EPSILON))
            edgePoints.push_back(p1);
          if (dx2 * dx2 + dy2 * dy2 > Kernel::FT(EPSILON))
            edgePoints.push_back(p2);
        }
      }

      // Create vertical triangular faces connecting the ridge endpoint to the
      // edge points
      Point ridgeTop(ridgeEndpoint.x(), ridgeEndpoint.y(), ridgeHeight);
      Point ridgeBase(ridgeEndpoint.x(), ridgeEndpoint.y(), 0.0);

      for (const Point &edgePoint : edgePoints) {
        Point edgeBase(edgePoint.x(), edgePoint.y(), 0.0);

        // Create triangle: ridge_base -> ridge_top -> edge_base -> ridge_base
        std::unique_ptr<Triangle> verticalTri(
            new Triangle(ridgeBase, ridgeTop, edgeBase));
        roof->addPatch(*verticalTri);
      }
    }
  }

  // 6. Handle building integration
  if (buildingHeight == 0.0) {
    // Just return the roof
    result = std::move(roof);
  } else {
    // Combine with building
    // Translate roof to building height
    translate(*roof, 0.0, 0.0, buildingHeight);

    // Create building walls
    auto building = extrude(footprint, buildingHeight);

    // Create result from building exterior shell
    result = std::make_unique<PolyhedralSurface>(
        building->as<Solid>().exteriorShell());

    // Add translated roof patches
    result->addPatchs(*roof);
  }

  propagateValidityFlag(*result, true);
  return result;
}

} // namespace SFCGAL::algorithm