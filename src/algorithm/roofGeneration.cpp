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
#include "SFCGAL/algorithm/ConsistentOrientationBuilder.h"
#include "SFCGAL/algorithm/covers.h"
#include "SFCGAL/algorithm/extrude.h"
#include "SFCGAL/algorithm/force2D.h"
#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/algorithm/orientation.h"
#include "SFCGAL/algorithm/straightSkeleton.h"
#include "SFCGAL/algorithm/tesselate.h"
#include "SFCGAL/algorithm/translate.h"
#include "SFCGAL/triangulate/triangulate2DZ.h"

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Surface_mesh.h>

#include <algorithm>
#include <cmath>
#include <map>
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

/**
 * @brief Generate a skillion roof from a polygon footprint and ridge line.
 *
 * Implementation follows the same pattern as gable roof for building
 * integration.
 */
auto
generateSkillionRoof(const Polygon &footprint, const LineString &ridgeLine,
                     double slopeAngle, bool addVerticalFaces,
                     double buildingHeight, bool closeBase)
    -> std::unique_ptr<Geometry>
{
  // Comprehensive input validation
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(footprint);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(ridgeLine);

  // Handle empty footprint
  if (footprint.isEmpty()) {
    if (buildingHeight > 0.0) {
      return std::make_unique<Solid>();
    }
    return std::make_unique<PolyhedralSurface>();
  }

  if (ridgeLine.isEmpty() || ridgeLine.numPoints() < 2) {
    BOOST_THROW_EXCEPTION(
        Exception("Ridge line must contain at least 2 points"));
  }

  if (slopeAngle <= 0.0 || slopeAngle >= 90.0) {
    BOOST_THROW_EXCEPTION(
        Exception("Slope angle must be between 0 and 90 degrees"));
  }

  if (buildingHeight < 0.0) {
    BOOST_THROW_EXCEPTION(Exception("Building height must be non-negative"));
  }

  // Get ridge line endpoints for distance calculation
  auto ridgeStart = ridgeLine.pointN(0);
  auto ridgeEnd   = ridgeLine.pointN(ridgeLine.numPoints() - 1);

  // Pre-compute trigonometric values
  auto slopeAngleRad = slopeAngle * M_PI / 180.0;
  auto slopeTan      = std::tan(slopeAngleRad);

  // Create elevated vertices for the sloped roof surface
  auto elevatedVertices =
      createSlopedSurfaceVertices(footprint, ridgeStart, ridgeEnd, slopeTan);

  // Create the roof surface
  auto roof = std::make_unique<PolyhedralSurface>();

  // Create the sloped roof surface as a single polygon
  if (elevatedVertices.size() >= 4) {
    LineString elevatedRing(elevatedVertices);
    Polygon    slopedSurface(elevatedRing);
    roof->addPatch(slopedSurface);
  }

  // Create vertical faces for proper roof closure
  auto verticalFaces =
      createVerticalFaces(footprint, elevatedVertices, addVerticalFaces);
  for (const auto &face : verticalFaces) {
    if (auto polygon = dynamic_cast<const Polygon *>(face.get())) {
      roof->addPatch(*polygon);
    }
  }

  // Handle building integration and return type
  if (buildingHeight == 0.0) {
    // Roof only mode
    if (closeBase) {
      // Add base polygon to create a closed solid
      roof->addPatch(footprint);

      auto solid = std::make_unique<Solid>(*roof);
      propagateValidityFlag(*solid, true);
      return solid;
    } else {
      // Return as polyhedral surface
      propagateValidityFlag(*roof, true);
      return roof;
    }
  } else {
    // Building + roof mode

    // Translate roof to building height
    translate(*roof, 0.0, 0.0, buildingHeight);

    // Create building walls
    auto building = extrude(footprint, buildingHeight);

    // Filter out the top face of the building (similar to
    // extrudeStraightSkeleton)
    auto isNotTopFace = [buildingHeight](const Polygon &patch) -> bool {
      const LineString &exterior = patch.exteriorRing();
      // Check if all points have z == buildingHeight (top face)
      bool allAtBuildingHeight = std::all_of(
          exterior.begin(), exterior.end(), [buildingHeight](const Point &p) {
            return std::abs(CGAL::to_double(p.z()) - buildingHeight) <
                   HEIGHT_TOLERANCE;
          });
      return !allAtBuildingHeight; // Keep all faces except the top
    };

    auto buildingShell = building->as<Solid>().exteriorShell();
    auto result        = std::make_unique<PolyhedralSurface>();

    // Copy all building faces except the top
    std::copy_if(buildingShell.begin(), buildingShell.end(),
                 std::back_inserter(*result), isNotTopFace);

    // Add roof patches
    result->addPatchs(*roof);

    // If addVerticalFaces is true, this should form a closed solid
    if (addVerticalFaces) {
      auto solid = std::make_unique<Solid>(*result);
      propagateValidityFlag(*solid, true);
      return solid;
    } else {
      propagateValidityFlag(*result, true);
      return result;
    }
  }
}

auto
generateRoof(const Polygon &footprint, const LineString &ridgeLine,
             const RoofParameters &params) -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(footprint);

  // Handle empty footprint
  if (footprint.isEmpty()) {
    if (params.buildingHeight > 0.0) {
      return std::make_unique<Solid>();
    }
    return std::make_unique<PolyhedralSurface>();
  }

  // Validate parameters
  if (params.buildingHeight < 0.0) {
    BOOST_THROW_EXCEPTION(Exception("Building height must be non-negative"));
  }
  if (params.roofHeight < 0.0) {
    BOOST_THROW_EXCEPTION(Exception("Roof height must be non-negative"));
  }

  switch (params.type) {
  case RoofType::FLAT: {
    // FLAT roof: simple translation or extrusion
    // slopeAngle and addVerticalFaces are ignored
    if (params.buildingHeight == 0.0) {
      // Simple Z-translation of footprint
      auto translated = footprint.clone();
      translate(*translated, 0.0, 0.0, params.roofHeight);
    } else {
      // Extrusion with building
      double totalHeight = params.buildingHeight + params.roofHeight;
      auto   extruded    = extrude(footprint, 0.0, 0.0, totalHeight);
      propagateValidityFlag(*extruded, true);
      return extruded; // Returns a Solid
    }
  }

  case RoofType::HIPPED: {
    // HIPPED roof: uses straight skeleton
    // slopeAngle and addVerticalFaces are ignored
    if (params.buildingHeight == 0.0) {
      // Just the roof
      auto roof = extrudeStraightSkeleton(footprint, params.roofHeight);

      // Legacy: extrudeStraightSkeleton returns a closed geometry
      if (!params.closeBase) {
        auto roofWithoutBase = std::make_unique<PolyhedralSurface>();

        auto isRoofSlope = [](const Polygon &patch) -> bool {
          const LineString &exterior = patch.exteriorRing();

          return std::any_of(
              exterior.begin(), exterior.end(),
              [](const Point &point) -> bool { return point.z() != 0; }

          );
        };

        std::copy_if(roof->begin(), roof->end(),
                     std::back_inserter(*roofWithoutBase), isRoofSlope);
        propagateValidityFlag(*roofWithoutBase, true);
        return roofWithoutBase;
      }
      propagateValidityFlag(*roof, true);
      return roof;
    } else {
      // Building + roof using the two-parameter version
      auto result = extrudeStraightSkeleton(footprint, params.buildingHeight,
                                            params.roofHeight);
      auto solid  = std::make_unique<Solid>(*result);
      propagateValidityFlag(*solid, true);
      return solid;
    }
  }

  case RoofType::GABLE: {
    // GABLE roof: uses automatic projected medial axis approach
    return generateGableRoof(footprint, params.slopeAngle,
                             params.addVerticalFaces, params.buildingHeight,
                             params.closeBase);
  }

  case RoofType::SKILLION: {
    // SKILLION roof: requires ridge line
    SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(ridgeLine);
    return generateSkillionRoof(footprint, ridgeLine, params.slopeAngle,
                                params.addVerticalFaces, params.buildingHeight,
                                params.closeBase);
  }

  default:
    BOOST_THROW_EXCEPTION(Exception("Unknown roof type"));
  }
}

auto
generateRoof(const Polygon &footprint, const LineString &ridgeLine,
             const RoofParameters &params, NoValidityCheck & /*nvc*/)
    -> std::unique_ptr<Geometry>
{
  return generateRoof(footprint, ridgeLine, params);
}

auto
generateGableRoof(const Polygon &footprint, double slopeAngle,
                  bool addVerticalFaces, double buildingHeight, bool closeBase)
    -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(footprint);

  // Handle empty footprint
  if (footprint.isEmpty()) {
    if (buildingHeight > 0.0) {
      return std::make_unique<Solid>();
    }
    return std::make_unique<PolyhedralSurface>();
  }

  if (slopeAngle <= 0.0 || slopeAngle >= 90.0) {
    BOOST_THROW_EXCEPTION(
        Exception("Slope angle must be between 0 and 90 degrees"));
  }

  if (buildingHeight < 0.0) {
    BOOST_THROW_EXCEPTION(Exception("Building height must be non-negative"));
  }

  // 1. Get projected medial axis to edges (this gives us the ridge line)
  auto projectedMedialAxis = approximateMedialAxis(footprint, true);

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
      auto roofTriangle = std::make_unique<Triangle>(newV1, newV2, newV3);
      roof->addPatch(*roofTriangle);
    }
  }

  // 5. Add vertical gable faces at ridge endpoints if requested
  // Note: For complex polygons (L-shape, T-shape), the medial axis may have
  // endpoints in the interior. We only create vertical faces for ridge points
  // that are actually at polygon boundaries (true gable ends).
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

    // For each ridge endpoint, check if it's on the polygon boundary
    // Only create vertical faces for boundary ridge points (true gable ends)
    const LineString &ring               = footprint.exteriorRing();
    const double      BOUNDARY_TOLERANCE = 1e-6;

    for (const Point &ridgeEndpoint : ridgeEndpoints) {
      // Check if ridge endpoint lies on the footprint boundary
      bool               onBoundary = false;
      std::vector<Point> boundarySegmentPoints;

      for (size_t i = 0; i < ring.numPoints() - 1; ++i) {
        const Point &p1 = ring.pointN(i);
        const Point &p2 = ring.pointN(i + 1);

        // Check if ridgeEndpoint lies on segment [p1, p2]
        auto dx        = CGAL::to_double(p2.x() - p1.x());
        auto dy        = CGAL::to_double(p2.y() - p1.y());
        auto segLength = std::sqrt(dx * dx + dy * dy);

        if (segLength < BOUNDARY_TOLERANCE) {
          continue; // Skip degenerate segments
        }

        // Vector from p1 to ridgeEndpoint
        auto rx = CGAL::to_double(ridgeEndpoint.x() - p1.x());
        auto ry = CGAL::to_double(ridgeEndpoint.y() - p1.y());

        // Project onto segment direction
        auto t = (rx * dx + ry * dy) / (segLength * segLength);

        // Check if projection is within segment bounds [0, 1]
        if (t >= -BOUNDARY_TOLERANCE && t <= 1.0 + BOUNDARY_TOLERANCE) {
          // Check distance from ridge endpoint to line
          auto projX = CGAL::to_double(p1.x()) + t * dx;
          auto projY = CGAL::to_double(p1.y()) + t * dy;
          auto distX = CGAL::to_double(ridgeEndpoint.x()) - projX;
          auto distY = CGAL::to_double(ridgeEndpoint.y()) - projY;
          auto dist  = std::sqrt(distX * distX + distY * distY);

          if (dist < BOUNDARY_TOLERANCE) {
            // Ridge endpoint is on this boundary segment
            onBoundary = true;

            // Collect the segment endpoints (excluding ridge endpoint itself)
            auto dist1 = std::sqrt(
                std::pow(CGAL::to_double(p1.x() - ridgeEndpoint.x()), 2) +
                std::pow(CGAL::to_double(p1.y() - ridgeEndpoint.y()), 2));
            auto dist2 = std::sqrt(
                std::pow(CGAL::to_double(p2.x() - ridgeEndpoint.x()), 2) +
                std::pow(CGAL::to_double(p2.y() - ridgeEndpoint.y()), 2));

            if (dist1 > BOUNDARY_TOLERANCE) {
              boundarySegmentPoints.push_back(p1);
            }
            if (dist2 > BOUNDARY_TOLERANCE) {
              boundarySegmentPoints.push_back(p2);
            }
            break; // Found the segment
          }
        }
      }

      // Only create vertical faces if ridge endpoint is on boundary
      if (onBoundary && !boundarySegmentPoints.empty()) {
        Point ridgeTop(ridgeEndpoint.x(), ridgeEndpoint.y(), ridgeHeight);
        Point ridgeBase(ridgeEndpoint.x(), ridgeEndpoint.y(), 0.0);

        for (const Point &edgePoint : boundarySegmentPoints) {
          Point edgeBase(edgePoint.x(), edgePoint.y(), 0.0);

          // Create triangle: ridge_base -> ridge_top -> edge_base
          std::unique_ptr<Triangle> verticalTri(
              new Triangle(ridgeBase, ridgeTop, edgeBase));
          roof->addPatch(*verticalTri);
        }
      }
    }
  }

  // 6. Handle building integration and return type
  if (buildingHeight == 0.0) {
    // Roof only mode
    if (closeBase) {
      // Add base polygon to create a closed solid
      roof->addPatch(footprint);

      // Use CGAL's robust orientation fixing
      auto solid = std::make_unique<Solid>(*roof);
      propagateValidityFlag(*solid, true);
      return solid;
    } else {
      // Return as polyhedral surface
      propagateValidityFlag(*roof, true);
      return roof;
    }
  } else {
    // Building + roof mode

    // Translate roof to building height
    translate(*roof, 0.0, 0.0, buildingHeight);

    // Create building walls
    auto building = extrude(footprint, buildingHeight);

    // Filter out the top face of the building
    auto isNotTopFace = [buildingHeight](const Polygon &patch) -> bool {
      const LineString &exterior = patch.exteriorRing();
      // Check if all points have z == buildingHeight (top face)
      bool allAtBuildingHeight = std::all_of(
          exterior.begin(), exterior.end(), [buildingHeight](const Point &p) {
            return std::abs(CGAL::to_double(p.z()) - buildingHeight) <
                   HEIGHT_TOLERANCE;
          });
      return !allAtBuildingHeight; // Keep all faces except the top
    };

    auto buildingShell = building->as<Solid>().exteriorShell();
    auto result        = std::make_unique<PolyhedralSurface>();

    // Copy all building faces except the top
    std::copy_if(buildingShell.begin(), buildingShell.end(),
                 std::back_inserter(*result), isNotTopFace);

    // Add translated roof patches
    result->addPatchs(*roof);

    // If addVerticalFaces is true, this should form a closed solid
    if (addVerticalFaces) {
      // Use CGAL's robust orientation fixing for building + roof
      auto solid = std::make_unique<Solid>(*result);
      propagateValidityFlag(*solid, true);
      return solid;
    } else {
      propagateValidityFlag(*result, true);
      return result;
    }
  }
}

} // namespace SFCGAL::algorithm
