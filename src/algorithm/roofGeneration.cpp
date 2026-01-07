// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/roofGeneration.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
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
#include "SFCGAL/detail/GetPointsVisitor.h"
#include "SFCGAL/detail/algorithm/FaceFilters.h"
#include "SFCGAL/triangulate/triangulate2DZ.h"

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/centroid.h>
#include <CGAL/squared_distance_2.h>

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

/**
 * @brief Orient a PolyhedralSurface to have consistent face orientations
 * Uses CGAL's Surface_mesh and polygon mesh processing to fix orientations
 */
// NOLINTBEGIN(readability-function-cognitive-complexity)
auto
orientPolyhedralSurface(PolyhedralSurface &surface) -> void
{
  if (surface.numPatches() == 0) {
    return;
  }

  // Convert to polygon soup
  using PointK      = CGAL::Point_3<Kernel>;
  using SurfaceMesh = CGAL::Surface_mesh<PointK>;

  std::vector<PointK>              points;
  std::vector<std::vector<size_t>> polygons;
  std::map<Coordinate, size_t>     coordToIndex;

  auto getPointIndex = [&](const Point &pt) -> size_t {
    const auto &coord = pt.coordinate();
    auto        it    = coordToIndex.find(coord);
    if (it != coordToIndex.end()) {
      return it->second;
    }
    size_t idx = points.size();
    points.emplace_back(pt.x(), pt.y(), pt.z());
    coordToIndex[coord] = idx;
    return idx;
  };

  for (size_t i = 0; i < surface.numPatches(); ++i) {
    const auto &patch = surface.patchN(i);
    if (const auto *triangle = dynamic_cast<const Triangle *>(&patch)) {
      std::vector<size_t> indices;
      indices.push_back(getPointIndex(triangle->vertex(0)));
      indices.push_back(getPointIndex(triangle->vertex(1)));
      indices.push_back(getPointIndex(triangle->vertex(2)));
      polygons.push_back(indices);
    } else {
      const auto         &ring = patch.exteriorRing();
      std::vector<size_t> indices;
      indices.reserve(ring.numPoints() - 1);
      for (size_t j = 0; j < ring.numPoints() - 1; ++j) {
        indices.push_back(getPointIndex(ring.pointN(j)));
      }
      polygons.push_back(indices);
    }
  }

  // Orient the polygon soup to be orientable
  bool orientable =
      CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);

  if (!orientable) {
    return; // Surface is not orientable, keep original
  }

  // Try to convert to Surface_mesh for better orientation
  SurfaceMesh mesh;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons,
                                                              mesh);

  if (CGAL::is_valid_polygon_mesh(mesh)) {
    // Orient to have consistent outward normals
    CGAL::Polygon_mesh_processing::orient(mesh);

    // Convert back to PolyhedralSurface
    surface = PolyhedralSurface();
    for (auto face : mesh.faces()) {
      std::vector<Point> pts;
      for (auto vertex : mesh.vertices_around_face(mesh.halfedge(face))) {
        const auto &meshPt = mesh.point(vertex);
        pts.emplace_back(meshPt.x(), meshPt.y(), meshPt.z());
      }
      if (pts.size() == 3) {
        surface.addPatch(Triangle(pts[0], pts[1], pts[2]));
      } else {
        pts.push_back(pts[0]); // Close the ring
        surface.addPatch(Polygon(LineString(pts)));
      }
    }
  } else {
    // Fallback: rebuild from oriented soup
    surface = PolyhedralSurface();
    for (const auto &polyIndices : polygons) {
      if (polyIndices.size() == 3) {
        surface.addPatch(Triangle(
            Point(points[polyIndices[0]].x(), points[polyIndices[0]].y(),
                  points[polyIndices[0]].z()),
            Point(points[polyIndices[1]].x(), points[polyIndices[1]].y(),
                  points[polyIndices[1]].z()),
            Point(points[polyIndices[2]].x(), points[polyIndices[2]].y(),
                  points[polyIndices[2]].z())));
      } else {
        std::vector<Point> pts;
        pts.reserve(polyIndices.size() + 1);
        for (size_t idx : polyIndices) {
          pts.emplace_back(points[idx].x(), points[idx].y(), points[idx].z());
        }
        pts.push_back(pts[0]); // Close the ring
        surface.addPatch(Polygon(LineString(pts)));
      }
    }
  }
}
// NOLINTEND(readability-function-cognitive-complexity)

/**
 * @brief Squared perpendicular distance from a point to a 2D line
 */
auto
squaredDistanceToLine2D(const Point &point, const Point &lineStart,
                        const Point &lineEnd) -> Kernel::FT
{
  Point_2 point2D(point.x(), point.y());
  Point_2 lineStart2D(lineStart.x(), lineStart.y());
  Point_2 lineEnd2D(lineEnd.x(), lineEnd.y());

  if (lineStart2D == lineEnd2D) {
    return CGAL::squared_distance(point2D, lineStart2D);
  }

  CGAL::Line_2<Kernel> line(lineStart2D, lineEnd2D);
  return CGAL::squared_distance(point2D, line);
}

/**
 * @brief Perpendicular distance from a point to a 2D line
 */
auto
distanceToLine2D(const Point &point, const Point &lineStart,
                 const Point &lineEnd) -> double
{
  return std::sqrt(
      CGAL::to_double(squaredDistanceToLine2D(point, lineStart, lineEnd)));
}

/**
 * @brief Minimum distance from a point to the polygon boundary (2D)
 *
 * Calculates the minimum perpendicular distance from a point to any edge
 * of the polygon's exterior ring.
 *
 * @param point The point to measure from
 * @param polygon The polygon whose boundary to measure to
 * @return Minimum distance to the polygon boundary
 */
auto
distanceToPolygonBoundary2D(const Point &point, const Polygon &polygon)
    -> Kernel::FT
{
  const auto &ring     = polygon.exteriorRing();
  size_t      numEdges = ring.numPoints() - 1;

  Point_2 point2D(point.x(), point.y());

  // Initialize with first edge distance
  const auto &firstStart = ring.pointN(0);
  const auto &firstEnd   = ring.pointN(1);
  Point_2     firstStart2D(firstStart.x(), firstStart.y());
  Point_2     firstEnd2D(firstEnd.x(), firstEnd.y());
  Segment_2   firstSegment(firstStart2D, firstEnd2D);
  auto        minSqDist = CGAL::squared_distance(point2D, firstSegment);

  for (size_t i = 1; i < numEdges; ++i) {
    const auto &edgeStart = ring.pointN(i);
    const auto &edgeEnd   = ring.pointN(i + 1);

    // Use segment distance, not line distance
    Point_2   start2D(edgeStart.x(), edgeStart.y());
    Point_2   end2D(edgeEnd.x(), edgeEnd.y());
    Segment_2 segment(start2D, end2D);
    auto      sqDist = CGAL::squared_distance(point2D, segment);

    if (sqDist < minSqDist) {
      minSqDist = sqDist;
    }
  }

  // sqrt requires conversion to double for exact kernels
  return {std::sqrt(CGAL::to_double(minSqDist))};
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
    const auto &point = ridgeLine.pointN(i);
    ridgePoints.emplace_back(point.x(), point.y(), slopeUp ? ridgeHeight : 0.0);
  }

  if (ridgePoints.empty()) {
    return surface; // Return empty surface if no ridge points
  }

  // Create simple triangulated surface
  const auto &baseRing      = basePolygon.exteriorRing();
  size_t      numBasePoints = baseRing.numPoints() - 1; // Exclude closing point

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
  const auto &exteriorRing = footprint.exteriorRing();
  size_t      numPoints    = exteriorRing.numPoints();

  if (numPoints < 4) { // Need at least 3 + closing point
    BOOST_THROW_EXCEPTION(Exception("Polygon must have at least 3 vertices"));
  }

  std::vector<Point> elevatedVertices;
  elevatedVertices.reserve(numPoints);

  // Calculate elevation for each vertex based on distance from ridge line
  for (size_t i = 0; i < numPoints - 1; ++i) { // Skip closing point
    const auto &vertex = exteriorRing.pointN(i);

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
 * @brief Compute polygon centroid (2D)
 */
auto
computePolygonCentroid2D(const Polygon &polygon) -> Point_2
{
  const auto          &ring = polygon.exteriorRing();
  std::vector<Point_2> points;
  points.reserve(ring.numPoints() - 1);
  for (size_t i = 0; i < ring.numPoints() - 1; ++i) {
    points.emplace_back(ring.pointN(i).x(), ring.pointN(i).y());
  }
  return CGAL::centroid(points.begin(), points.end(), CGAL::Dimension_tag<0>());
}

/**
 * @brief Create vertical faces connecting base polygon to elevated roof
 * with proper face orientations (normals pointing outward)
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

  const auto &exteriorRing = footprint.exteriorRing();
  size_t      numPoints    = exteriorRing.numPoints();

  // Get polygon centroid for orientation check
  Point_2 polyCentroid = computePolygonCentroid2D(footprint);

  for (size_t i = 0; i < numPoints - 1; ++i) {
    size_t nextI = (i + 1) % (numPoints - 1);

    // Base edge points (at z=0)
    const auto &baseVertex1 = exteriorRing.pointN(i);
    const auto &baseVertex2 = exteriorRing.pointN(nextI);
    Kernel::FT  zBase(0);
    Point       basePt1(baseVertex1.x(), baseVertex1.y(), zBase);
    Point       basePt2(baseVertex2.x(), baseVertex2.y(), zBase);

    // Corresponding elevated points on the roof
    Point roofPt1 = elevatedVertices[i];
    Point roofPt2 = elevatedVertices[nextI];

    // Check if there's height difference for this edge
    bool hasHeightDiff1 = roofPt1.z() != zBase;
    bool hasHeightDiff2 = roofPt2.z() != zBase;

    if (!hasHeightDiff1 && !hasHeightDiff2) {
      continue; // No height difference, skip
    }

    if (hasHeightDiff1 && hasHeightDiff2) {
      // Quadrilateral face: split into two triangles using diagonal
      // basePt2-roofPt1 This ensures roof-adjacent edge (roofPt2->roofPt1) is
      // reversed from roof edge T1: basePt1 -> basePt2 -> roofPt1 T2: basePt2
      // -> roofPt2 -> roofPt1
      faces.push_back(std::make_unique<Triangle>(basePt1, basePt2, roofPt1));
      faces.push_back(std::make_unique<Triangle>(basePt2, roofPt2, roofPt1));
    } else if (hasHeightDiff1) {
      // Triangle connecting base edge to roofPt1
      // Roof-adjacent edge: roofPt1->basePt1 is reversed from roof edge
      // basePt1->roofPt1 (Note: when hasHeightDiff2 is false, roofPt2.z ==
      // basePt2.z, so roof edge
      //  effectively goes from roofPt1 back towards the ridge)
      faces.push_back(std::make_unique<Triangle>(basePt1, basePt2, roofPt1));
    } else {
      // Triangle connecting base edge to roofPt2
      // Roof-adjacent edge: roofPt2->basePt1 is reversed from roof edge
      // basePt1->roofPt2 (Note: when hasHeightDiff1 is false, basePt1 ==
      // roofPt1, so this works)
      faces.push_back(std::make_unique<Triangle>(basePt1, basePt2, roofPt2));
    }
  }

  return faces;
}

} // anonymous namespace

/// @} end of private section

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
    if (const auto *triangle = dynamic_cast<const Triangle *>(face.get())) {
      roof->addPatch(*triangle);
    } else if (const auto *polygon =
                   dynamic_cast<const Polygon *>(face.get())) {
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
    }
    // Return as polyhedral surface
    propagateValidityFlag(*roof, true);
    return roof;
  }
  {
    // Building + roof mode

    // Translate roof to building height
    translate(*roof, 0.0, 0.0, buildingHeight);

    // Create building walls
    auto building = extrude(footprint, buildingHeight);

    // Filter out the top face of the building (similar to
    // extrudeStraightSkeleton)
    auto buildingShell = building->as<Solid>().exteriorShell();
    auto result        = std::make_unique<PolyhedralSurface>();

    // Copy all building faces except the top
    std::copy_if(buildingShell.begin(), buildingShell.end(),
                 std::back_inserter(*result),
                 [buildingHeight](const Polygon &patch) -> bool {
                   return detail::isNotFaceAtHeight(patch, buildingHeight);
                 });

    // Add roof patches
    result->addPatchs(*roof);

    // Orient the surface to have consistent face orientations
    orientPolyhedralSurface(*result);

    // If addVerticalFaces is true, this should form a closed solid
    if (addVerticalFaces) {
      auto solid = std::make_unique<Solid>(*result);
      propagateValidityFlag(*solid, true);
      return solid;
    }
    propagateValidityFlag(*result, true);
    return result;
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
      propagateValidityFlag(*translated, true);
      return translated;
    }
    // Extrusion with building
    double totalHeight = params.buildingHeight + params.roofHeight;
    auto   extruded    = extrude(footprint, 0.0, 0.0, totalHeight);
    propagateValidityFlag(*extruded, true);
    return extruded; // Returns a Solid
  }

  case RoofType::HIPPED: {
    // HIPPED roof: uses straight skeleton
    // slopeAngle and addVerticalFaces are ignored
    if (params.buildingHeight == 0.0) {
      // Just the roof
      auto roof =
          extrudeStraightSkeleton(footprint, params.roofHeight, params.angles);

      // Legacy: extrudeStraightSkeleton returns a closed geometry
      if (!params.closeBase) {
        auto roofWithoutBase = std::make_unique<PolyhedralSurface>();

        std::copy_if(roof->begin(), roof->end(),
                     std::back_inserter(*roofWithoutBase), detail::isRoofSlope);
        propagateValidityFlag(*roofWithoutBase, true);
        return roofWithoutBase;
      }
      propagateValidityFlag(*roof, true);
      return roof;
    }
    // Building + roof using the two-parameter version
    auto result = extrudeStraightSkeleton(footprint, params.buildingHeight,
                                          params.roofHeight, params.angles);
    auto solid  = std::make_unique<Solid>(*result);
    propagateValidityFlag(*solid, true);
    return solid;
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

// NOLINTBEGIN(readability-function-cognitive-complexity)
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

  // 2. Create a MultiPoint containing all points from footprint and
  // projectedMedialAxis
  auto multiPoint = std::make_unique<MultiPoint>();

  // Use GetPointsVisitor to collect all points from footprint and
  // projectedMedialAxis
  SFCGAL::detail::GetPointsVisitor visitor;

  // Collect points from footprint
  visitor.visit(footprint);

  // Collect points from projected medial axis
  visitor.visit(*projectedMedialAxis);

  // Add all collected points to the MultiPoint
  for (const Point *point : visitor.points) {
    multiPoint->addGeometry(*point);
  }

  // Perform constrained Delaunay triangulation on the MultiPoint
  auto triangulation = SFCGAL::triangulate::triangulate2DZ(*multiPoint);

  // Add constraints from projected medial axis
  auto projectedMedialAxis_clone = projectedMedialAxis->clone();
  SFCGAL::algorithm::force2D(*projectedMedialAxis_clone);
  SFCGAL::triangulate::triangulate2DZ(*projectedMedialAxis_clone,
                                      triangulation);

  auto triangulated_surface = triangulation.getTriangulatedSurface();

  // Calculate slope tangent for height computation
  // Height = distance_to_boundary * tan(angle)
  double slopeTan = std::tan(slopeAngle * M_PI / 180.0);

  // Helper to check if a point is on the polygon boundary
  auto isOnBoundary = [&footprint](const Point &pt) -> bool {
    const auto &ring = footprint.exteriorRing();
    Point_2     pt2D(pt.x(), pt.y());
    for (size_t i = 0; i < ring.numPoints() - 1; ++i) {
      Point_2   segPt1(ring.pointN(i).x(), ring.pointN(i).y());
      Point_2   segPt2(ring.pointN(i + 1).x(), ring.pointN(i + 1).y());
      Segment_2 segment(segPt1, segPt2);
      auto      sqDist = CGAL::squared_distance(pt2D, segment);
      if (sqDist < Kernel::FT(EPSILON_SQ)) {
        return true;
      }
    }
    return false;
  };

  // 3. Collect ridge points with proper heights:
  // - Interior points: use their actual distance to boundary
  // - Boundary points (gable ends): use height of nearest interior point on
  // same segment
  std::map<std::pair<Kernel::FT, Kernel::FT>, Kernel::FT> ridgePointHeights;

  for (size_t i = 0; i < projectedMedialAxis->numGeometries(); ++i) {
    const auto *line =
        dynamic_cast<const LineString *>(&projectedMedialAxis->geometryN(i));
    if (line == nullptr || line->numPoints() == 0) {
      continue;
    }

    // First pass: calculate heights for interior points
    std::vector<Kernel::FT> pointDistances(line->numPoints());
    std::vector<bool>       pointOnBoundary(line->numPoints());

    for (size_t j = 0; j < line->numPoints(); ++j) {
      const auto &pt     = line->pointN(j);
      pointDistances[j]  = distanceToPolygonBoundary2D(pt, footprint);
      pointOnBoundary[j] = isOnBoundary(pt);
    }

    // Second pass: assign heights
    // For boundary points, find the nearest interior point's distance
    for (size_t j = 0; j < line->numPoints(); ++j) {
      const auto &pt  = line->pointN(j);
      auto        key = std::make_pair(pt.x(), pt.y());

      Kernel::FT height(0);
      if (!pointOnBoundary[j]) {
        // Interior point: use its actual distance
        height = pointDistances[j] * slopeTan;
      } else {
        // Boundary point: find nearest interior point's height
        Kernel::FT nearestInteriorDist(0);
        // Search forward
        for (size_t k = j + 1; k < line->numPoints(); ++k) {
          if (!pointOnBoundary[k]) {
            nearestInteriorDist = pointDistances[k];
            break;
          }
        }
        // If not found forward, search backward
        if (nearestInteriorDist == Kernel::FT(0) && j > 0) {
          for (size_t k = j; k > 0; --k) {
            if (!pointOnBoundary[k - 1]) {
              nearestInteriorDist = pointDistances[k - 1];
              break;
            }
          }
        }
        height = nearestInteriorDist * slopeTan;
      }

      // Only update if we have a valid height or no existing entry
      auto existingIt = ridgePointHeights.find(key);
      if (existingIt == ridgePointHeights.end() || height > Kernel::FT(0)) {
        ridgePointHeights[key] = height;
      }
    }
  }

  // Helper function to check if a point is a ridge point and get its height
  auto getRidgeHeight = [&ridgePointHeights](const Point &point) -> Kernel::FT {
    auto key = std::make_pair(point.x(), point.y());
    auto it  = ridgePointHeights.find(key);
    if (it != ridgePointHeights.end()) {
      return it->second;
    }
    // Fallback: search with tolerance
    for (const auto &entry : ridgePointHeights) {
      auto dx = entry.first.first - point.x();
      auto dy = entry.first.second - point.y();
      if (dx * dx + dy * dy < Kernel::FT(EPSILON)) {
        return entry.second;
      }
    }
    return {0};
  };

  auto isRidgePoint = [&ridgePointHeights](const Point &point) -> bool {
    auto key = std::make_pair(point.x(), point.y());
    if (ridgePointHeights.find(key) != ridgePointHeights.end()) {
      return true;
    }
    // Fallback: search with tolerance
    return std::any_of(ridgePointHeights.begin(), ridgePointHeights.end(),
                       [&point](const auto &entry) -> bool {
                         auto dx = entry.first.first - point.x();
                         auto dy = entry.first.second - point.y();
                         return dx * dx + dy * dy < Kernel::FT(EPSILON);
                       });
  };

  // Create roof surface
  auto roof = std::make_unique<PolyhedralSurface>();

  // 4. Process triangulation and elevate ridge points
  for (size_t i = 0; i < triangulated_surface->numPatches(); ++i) {
    const auto &patch = triangulated_surface->patchN(i);
    if (const auto *triangle = dynamic_cast<const Triangle *>(&patch)) {
      const auto &vertex1 = triangle->vertex(0);
      const auto &vertex2 = triangle->vertex(1);
      const auto &vertex3 = triangle->vertex(2);

      // Check if triangle centroid is inside the original polygon
      Point_2 centroid2D = CGAL::centroid(Point_2(vertex1.x(), vertex1.y()),
                                          Point_2(vertex2.x(), vertex2.y()),
                                          Point_2(vertex3.x(), vertex3.y()));
      Point   centroid(centroid2D.x(), centroid2D.y(), Kernel::FT(0));

      // Use SFCGAL's covers algorithm to check if centroid is inside footprint
      if (!SFCGAL::algorithm::covers(footprint, centroid)) {
        continue; // Skip triangles outside the original polygon
      }

      // Create new vertices with proper Z coordinates
      // Ridge points are elevated based on their distance to polygon boundary
      Point newVertex1(vertex1.x(), vertex1.y(), Kernel::FT(0));
      Point newVertex2(vertex2.x(), vertex2.y(), Kernel::FT(0));
      Point newVertex3(vertex3.x(), vertex3.y(), Kernel::FT(0));

      if (isRidgePoint(vertex1)) {
        newVertex1 = Point(vertex1.x(), vertex1.y(), getRidgeHeight(vertex1));
      }
      if (isRidgePoint(vertex2)) {
        newVertex2 = Point(vertex2.x(), vertex2.y(), getRidgeHeight(vertex2));
      }
      if (isRidgePoint(vertex3)) {
        newVertex3 = Point(vertex3.x(), vertex3.y(), getRidgeHeight(vertex3));
      }

      // Create roof triangle - ensure normal points upward (positive Z)
      Vector_3 edge1(newVertex2.x() - newVertex1.x(),
                     newVertex2.y() - newVertex1.y(),
                     newVertex2.z() - newVertex1.z());
      Vector_3 edge2(newVertex3.x() - newVertex1.x(),
                     newVertex3.y() - newVertex1.y(),
                     newVertex3.z() - newVertex1.z());
      Vector_3 normal = CGAL::cross_product(edge1, edge2);

      std::unique_ptr<Triangle> roofTriangle;
      if (normal.z() >= Kernel::FT(0)) {
        roofTriangle =
            std::make_unique<Triangle>(newVertex1, newVertex2, newVertex3);
      } else {
        roofTriangle =
            std::make_unique<Triangle>(newVertex1, newVertex3, newVertex2);
      }
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
      if ((line != nullptr) && line->numPoints() >= 2) {
        ridgeEndpoints.push_back(line->pointN(0));
        ridgeEndpoints.push_back(line->pointN(line->numPoints() - 1));
      }
    }

    // For each ridge endpoint, check if it's on the polygon boundary.
    // Only create vertical faces for boundary ridge points (true gable ends).
    const LineString &ring = footprint.exteriorRing();

    const Kernel::FT TOLERANCE_SQ(1e-12);

    // Build polygon centroid
    std::vector<Point_2> ringPoints2D;
    ringPoints2D.reserve(ring.numPoints() - 1);
    for (size_t i = 0; i < ring.numPoints() - 1; ++i) {
      ringPoints2D.emplace_back(ring.pointN(i).x(), ring.pointN(i).y());
    }
    Point_2 polyCentroid2D = CGAL::centroid(
        ringPoints2D.begin(), ringPoints2D.end(), CGAL::Dimension_tag<0>());

    for (const Point &ridgeEndpoint : ridgeEndpoints) {
      bool               onBoundary = false;
      std::vector<Point> boundarySegmentPoints;

      Point_2 ridgePoint2D(ridgeEndpoint.x(), ridgeEndpoint.y());

      for (size_t i = 0; i < ring.numPoints() - 1; ++i) {
        const Point &segmentPoint1 = ring.pointN(i);
        const Point &segmentPoint2 = ring.pointN(i + 1);

        Point_2 segStart2D(segmentPoint1.x(), segmentPoint1.y());
        Point_2 segEnd2D(segmentPoint2.x(), segmentPoint2.y());

        if (segStart2D == segEnd2D) {
          continue;
        }

        Segment_2  segment(segStart2D, segEnd2D);
        Kernel::FT sqDist = CGAL::squared_distance(ridgePoint2D, segment);

        if (sqDist < TOLERANCE_SQ) {
          onBoundary = true;

          Kernel::FT sqDistToStart =
              CGAL::squared_distance(ridgePoint2D, segStart2D);
          Kernel::FT sqDistToEnd =
              CGAL::squared_distance(ridgePoint2D, segEnd2D);

          if (sqDistToStart > TOLERANCE_SQ) {
            boundarySegmentPoints.push_back(segmentPoint1);
          }
          if (sqDistToEnd > TOLERANCE_SQ) {
            boundarySegmentPoints.push_back(segmentPoint2);
          }
          break;
        }
      }

      // Create ONE gable triangle connecting both corners to ridge top
      if (onBoundary && boundarySegmentPoints.size() == 2) {
        // Calculate ridge height using the constant segment height
        auto ridgeZ = getRidgeHeight(ridgeEndpoint);

        Kernel::FT zBase(0);
        Point      corner1(boundarySegmentPoints[0].x(),
                           boundarySegmentPoints[0].y(), zBase);
        Point      corner2(boundarySegmentPoints[1].x(),
                           boundarySegmentPoints[1].y(), zBase);
        Point      ridgeTop(ridgeEndpoint.x(), ridgeEndpoint.y(), ridgeZ);

        // Gable centroid
        Kernel::FT gableCentroidX =
            (ridgeTop.x() + corner1.x() + corner2.x()) / Kernel::FT(3);
        Kernel::FT gableCentroidY =
            (ridgeTop.y() + corner1.y() + corner2.y()) / Kernel::FT(3);

        // Outward direction from polygon centroid to gable centroid
        Kernel::FT outwardX = gableCentroidX - polyCentroid2D.x();
        Kernel::FT outwardY = gableCentroidY - polyCentroid2D.y();

        // Normal of triangle corner1 -> ridgeTop -> corner2
        Vector_3 edge1(ridgeTop.x() - corner1.x(), ridgeTop.y() - corner1.y(),
                       ridgeTop.z() - corner1.z());
        Vector_3 edge2(corner2.x() - corner1.x(), corner2.y() - corner1.y(),
                       corner2.z() - corner1.z());
        Vector_3 normal = CGAL::cross_product(edge1, edge2);

        // Check if normal points outward (XY dot product)
        Kernel::FT dotProduct = normal.x() * outwardX + normal.y() * outwardY;

        std::unique_ptr<Triangle> gableTri;
        if (dotProduct >= Kernel::FT(0)) {
          gableTri = std::make_unique<Triangle>(corner1, ridgeTop, corner2);
        } else {
          gableTri = std::make_unique<Triangle>(corner1, corner2, ridgeTop);
        }
        roof->addPatch(*gableTri);
      }
    }
  }

  // 6. Handle building integration and return type
  if (buildingHeight == 0.0) {
    // Roof only mode
    if (closeBase) {
      // Add base polygon to create a closed solid
      // The base face normal should point downward (negative Z) for a valid
      // solid. A counter-clockwise polygon (viewed from above) has upward
      // normal, so we need to reverse it.
      Polygon baseFace(footprint);
      baseFace.reverse();
      roof->addPatch(baseFace);

      auto solid = std::make_unique<Solid>(*roof);
      propagateValidityFlag(*solid, true);
      return solid;
    }
    // Return as polyhedral surface
    propagateValidityFlag(*roof, true);
    return roof;
  }
  {
    // Building + roof mode

    // Translate roof to building height
    translate(*roof, 0.0, 0.0, buildingHeight);

    // Create building walls
    auto building = extrude(footprint, buildingHeight);

    // Filter out the top face of the building
    auto buildingShell = building->as<Solid>().exteriorShell();
    auto result        = std::make_unique<PolyhedralSurface>();

    // Copy all building faces except the top
    std::copy_if(buildingShell.begin(), buildingShell.end(),
                 std::back_inserter(*result),
                 [buildingHeight](const Polygon &patch) -> bool {
                   return detail::isNotFaceAtHeight(patch, buildingHeight);
                 });

    // Add translated roof patches
    result->addPatchs(*roof);

    // Orient the surface to have consistent face orientations
    orientPolyhedralSurface(*result);

    // If addVerticalFaces is true, this should form a closed solid
    if (addVerticalFaces) {
      auto solid = std::make_unique<Solid>(*result);
      propagateValidityFlag(*solid, true);
      return solid;
    }
    propagateValidityFlag(*result, true);
    return result;
  }
}
// NOLINTEND(readability-function-cognitive-complexity)

} // namespace SFCGAL::algorithm
