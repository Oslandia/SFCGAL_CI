// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/straightSkeleton.h"
#include "SFCGAL/algorithm/distance.h"

#include "SFCGAL/Envelope.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"

#include "SFCGAL/Exception.h"

#include "SFCGAL/detail/algorithm/FaceFilters.h"
#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/algorithm/orientation.h"
#include "SFCGAL/algorithm/tesselate.h"
#include "SFCGAL/algorithm/translate.h"

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Straight_skeleton_converter_2.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/extrude_skeleton.h>

#include <cmath>
#include <limits>
#include <memory>

namespace SFCGAL::algorithm {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

using Straight_skeleton_2 = CGAL::Straight_skeleton_2<Kernel>;
using Arrangement_2 = CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>>;

#if CGAL_VERSION_MAJOR < 6
template <class T>
using SHARED_PTR = boost::shared_ptr<T>;
#else
template <class T>
using SHARED_PTR = std::shared_ptr<T>;
#endif

namespace { // anonymous

template <class K, bool outputDistanceInM>
void
straightSkeletonToMultiLineString(const CGAL::Straight_skeleton_2<K> &skeleton,
                                  MultiLineString &result, bool innerOnly,
                                  Kernel::Vector_2 &translate,
                                  const double     &toleranceAbs)
{
  using Ss = CGAL::Straight_skeleton_2<K>;

  using Vertex_const_handle     = typename Ss::Vertex_const_handle;
  using Halfedge_const_handle   = typename Ss::Halfedge_const_handle;
  using Halfedge_const_iterator = typename Ss::Halfedge_const_iterator;

  Halfedge_const_handle const null_halfedge;
  Vertex_const_handle const   null_vertex;

  for (Halfedge_const_iterator it = skeleton.halfedges_begin();
       it != skeleton.halfedges_end(); ++it) {
    // skip contour edge
    if (!it->is_bisector()) {
      continue;
    }

    // Skip non-inner edges if requested
    if (innerOnly && !it->is_inner_bisector()) {
      continue;
    }

    // avoid duplicates
    if (it->opposite() < it) {
      continue;
    }

    LineString *lineString = nullptr;
    Point       pa(it->opposite()->vertex()->point());
    Point       pb(it->vertex()->point());
    // avoid degenerate cases.https://gitlab.com/Oslandia/SFCGAL/-/issues/143
    if (pa != pb && distancePointPoint(pa, pb) > toleranceAbs) {
      if (outputDistanceInM) {
        pa.setM(CGAL::to_double(it->opposite()->vertex()->time()));
        pb.setM(CGAL::to_double(it->vertex()->time()));
        lineString = new LineString(pa, pb);
      } else {
        lineString = new LineString(pa, pb);
      }

      algorithm::translate(*lineString, translate);
      result.addGeometry(lineString);
    }
  }
}

/**
 * Return abc angle in radians
 */
auto
angle(const Point &pointA, const Point &pointB, const Point &pointC) -> double
{
  Point const ab(CGAL::to_double(pointB.x() - pointA.x()),
                 CGAL::to_double(pointB.y() - pointA.y()));
  Point const cb(CGAL::to_double(pointB.x() - pointC.x()),
                 CGAL::to_double(pointB.y() - pointC.y()));

  double const dot =
      (CGAL::to_double(ab.x() * cb.x() + ab.y() * cb.y())); /* dot product */
  double const cross =
      (CGAL::to_double(ab.x() * cb.y() - ab.y() * cb.x())); /* cross product */

  double const alpha = std::atan2(cross, dot);

  return alpha;
}

/**
 * @brief Project free endpoint to boundary using edge midpoint method
 *
 * For a free endpoint P on medial axis:
 * 1. Find skeleton bisectors from contour incident to P (exclude inner
 * bisectors)
 * 2. Collect BOTH defining contour edges from each bisector
 * 3. Find the common/repeated edge (appears in multiple bisectors)
 * 4. Return midpoint of that edge's two vertices
 *
 * Example: Rectangle corner with medial axis endpoint at (1.5, 1.5)
 * - Two corner bisectors converge: from (0,0) and (0,3)
 * - From (0,0) bisector: left edge + bottom edge
 * - From (0,3) bisector: left edge + top edge
 * - Left edge appears twice → midpoint: (0, 1.5)
 */
template <class K>
auto
findEdgeMidpointProjection(
    typename CGAL::Straight_skeleton_2<K>::Vertex_const_handle vertex) -> Point
{
  using Ss                    = CGAL::Straight_skeleton_2<K>;
  using Halfedge_const_handle = typename Ss::Halfedge_const_handle;

  // Collect all bisectors incident to this vertex
  std::vector<Halfedge_const_handle> incidentBisectors;

  // Traverse circulator around vertex
  Halfedge_const_handle start   = vertex->halfedge();
  Halfedge_const_handle current = start;

  do {
    if (current->is_bisector()) {
      incidentBisectors.push_back(current);
    }
    current = current->opposite()->prev();
  } while (current != start);

  // Collect defining contour edges and count occurrences
  // Use sorted pair to normalize edge representation
  std::map<std::pair<Point, Point>, int> edgeCount;

  for (const auto &bisector : incidentBisectors) {
    // Skip inner bisectors - only use skeleton bisectors from the contour
    if (bisector->is_inner_bisector()) {
      continue;
    }

    // Each bisector is defined by TWO contour edges
    // Get first defining edge
    auto  de1 = bisector->defining_contour_edge();
    Point v1_1(de1->vertex()->point());
    Point v1_2(de1->opposite()->vertex()->point());
    auto  edge1 =
        (v1_1.x() < v1_2.x() || (v1_1.x() == v1_2.x() && v1_1.y() < v1_2.y()))
             ? std::make_pair(v1_1, v1_2)
             : std::make_pair(v1_2, v1_1);
    edgeCount[edge1]++;

    // Get second defining edge (from opposite halfedge)
    auto  de2 = bisector->opposite()->defining_contour_edge();
    Point v2_1(de2->vertex()->point());
    Point v2_2(de2->opposite()->vertex()->point());
    auto  edge2 =
        (v2_1.x() < v2_2.x() || (v2_1.x() == v2_2.x() && v2_1.y() < v2_2.y()))
             ? std::make_pair(v2_1, v2_2)
             : std::make_pair(v2_2, v2_1);
    edgeCount[edge2]++;
  }

  // Find the most common edge (appears in multiple bisectors)
  std::pair<Point, Point> dominantEdge;
  int                     maxCount = 0;

  for (const auto &[edge, count] : edgeCount) {
    if (count > maxCount) {
      maxCount     = count;
      dominantEdge = edge;
    }
  }

  // Return midpoint of the dominant edge
  Kernel::FT midX = (dominantEdge.first.x() + dominantEdge.second.x()) / 2;
  Kernel::FT midY = (dominantEdge.first.y() + dominantEdge.second.y()) / 2;

  return {midX, midY};
}

template <class K>
void
straightSkeletonToMedialAxis(const CGAL::Straight_skeleton_2<K> &skeleton,
                             MultiLineString                    &result,
                             Kernel::Vector_2 &translate, bool projectToEdges)
{
  using Ss = CGAL::Straight_skeleton_2<K>;

  using Vertex_const_handle     = typename Ss::Vertex_const_handle;
  using Halfedge_const_handle   = typename Ss::Halfedge_const_handle;
  using Halfedge_const_iterator = typename Ss::Halfedge_const_iterator;

  // Maximum angle for touching bisectors to be
  // retained in the ouptput. The value is in
  // radians and represents the angle between the
  // bisector and the defining edge.
  //
  // The angle between the incident edges will be double
  // this value, so CGAL_PI / 8.0 means 45 degrees between
  // defining edges
  //
  // TODO: take as parameter ?
  //
  // NOTE: I'm adding some tolerance here to include those angles
  //       of exactly 45 degrees that are otherwise cut out due
  //       to rounding precision
  //
  const double maxTouchingAngle = CGAL_PI / 8.0 + 1e-13;

  // Build connectivity graph to identify free endpoints (degree 1 vertices)
  std::map<Vertex_const_handle, int> vertexDegree;

  // First pass: identify medial axis bisectors and count vertex degrees
  std::vector<Halfedge_const_handle> medialAxisBisectors;

  for (Halfedge_const_iterator it = skeleton.halfedges_begin();
       it != skeleton.halfedges_end(); ++it) {
    // skip contour edge
    if (!it->is_bisector()) {
      continue;
    }

    // avoid duplicates
    if (it->opposite() < it) {
      continue;
    }

    // Skip non-inner edges unless they bisect
    // a low angle pair of segments
    if (!it->is_inner_bisector()) {

      const Halfedge_const_handle &de1 = it->defining_contour_edge();

      // We need to check the angle formed by:
      const Point &pointA = it->vertex()->point();
      const Point &pointB = de1->vertex()->point();
      const Point &pointC = de1->opposite()->vertex()->point();

      double const ang = angle(pointA, pointB, pointC);

      if (ang > maxTouchingAngle) {
        continue;
      }
    }

    // This is a medial axis bisector
    medialAxisBisectors.push_back(it);
    vertexDegree[it->vertex()]++;
    vertexDegree[it->opposite()->vertex()]++;
  }

  // Second pass: create linestrings with optional projection
  for (const auto &it : medialAxisBisectors) {
    Point startPoint(it->opposite()->vertex()->point());
    Point endPoint(it->vertex()->point());

    std::unique_ptr<LineString> lineString(new LineString);

    // Project start point if it's a free endpoint
    if (projectToEdges && vertexDegree[it->opposite()->vertex()] == 1) {
      Point proj = findEdgeMidpointProjection<K>(it->opposite()->vertex());
      algorithm::translate(proj, translate);
      lineString->addPoint(proj);
    }

    // Add original medial axis segment
    algorithm::translate(startPoint, translate);
    algorithm::translate(endPoint, translate);
    lineString->addPoint(startPoint);
    lineString->addPoint(endPoint);

    // Project end point if it's a free endpoint
    if (projectToEdges && vertexDegree[it->vertex()] == 1) {
      Point proj = findEdgeMidpointProjection<K>(it->vertex());
      algorithm::translate(proj, translate);
      lineString->addPoint(proj);
    }

    result.addGeometry(lineString.release());
  }
}

auto
straightSkeleton(const Polygon_with_holes_2 &poly)
    -> SHARED_PTR<Straight_skeleton_2>
{
  SHARED_PTR<CGAL::Straight_skeleton_2<CGAL::Epick>> const skeletonEpick =
      CGAL::create_interior_straight_skeleton_2(
          poly.outer_boundary().vertices_begin(),
          poly.outer_boundary().vertices_end(), poly.holes_begin(),
          poly.holes_end(), CGAL::Epick());
  SHARED_PTR<Straight_skeleton_2> ret;
  if (skeletonEpick) {
    ret =
        CGAL::convert_straight_skeleton_2<Straight_skeleton_2>(*skeletonEpick);
  }
  return ret;
}

// Throw an exception if any two polygon rings touch,
// since CGAL segfaults in that case.
// @todo find upstream reference to find out exact case
// See https://github.com/Oslandia/SFCGAL/issues/75
void
checkNoTouchingHoles(const Polygon &geom)
{
  const size_t numRings = geom.numRings();

  for (size_t ri = 0; ri < numRings - 1; ++ri) {
    for (size_t rj = ri + 1; rj < numRings; ++rj) {
      std::unique_ptr<Geometry> inter =
          geom.is3D() ? intersection3D(geom.ringN(ri), geom.ringN(rj))
                      : intersection(geom.ringN(ri), geom.ringN(rj));

      // @note this check would accept rings touching at
      //       more than a single point, which may be
      //       still dangerous. @todo find out !
      if (!inter->isEmpty() && inter->is<Point>()) {
        BOOST_THROW_EXCEPTION(NotImplementedException(
            std::string("straight skeleton of Polygon with point ") +
            "touching rings is not implemented. " + "Error at " +
            inter->asText()));
      }
    }
  }
}

auto
preparePolygon(const Polygon &poly, Kernel::Vector_2 &trans)
    -> Polygon_with_holes_2
{
  checkNoTouchingHoles(poly);
  Envelope const env = poly.envelope();
  trans              = Kernel::Vector_2(-env.xMin(), -env.yMin());

  // @todo: avoid cloning !
  std::unique_ptr<Polygon> cloned = poly.clone();
  algorithm::translate(*cloned, trans);
  Polygon_with_holes_2 ret = cloned->toPolygon_with_holes_2();
  trans                    = -trans;

  return ret;
}

void
extractPolygons(const Geometry &geom, std::vector<Polygon> &vect)
{
  switch (geom.geometryTypeId()) {
  case TYPE_TRIANGLE:
    vect.push_back(geom.as<Triangle>().toPolygon());
    break;
  case TYPE_POLYGON:
    vect.push_back(geom.as<Polygon>());
    break;
  case TYPE_MULTIPOLYGON: {
    const auto &multiPolygon = geom.as<MultiPolygon>();
    for (size_t i = 0; i < multiPolygon.numGeometries(); i++) {
      vect.push_back(multiPolygon.polygonN(i));
    }
    break;
  }
  default:
    break;
  }
}

} // namespace

/// @} end of private section

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

auto
straightSkeleton(const Geometry &geom, bool       autoOrientation,
                 NoValidityCheck /*unused*/, bool innerOnly,
                 bool outputDistanceInM, const double & /*toleranceAbs*/)
    -> std::unique_ptr<MultiLineString>
{
  switch (geom.geometryTypeId()) {
  case TYPE_TRIANGLE:
    return straightSkeleton(geom.as<Triangle>().toPolygon(), autoOrientation,
                            innerOnly, outputDistanceInM);

  case TYPE_POLYGON:
    return straightSkeleton(geom.as<Polygon>(), autoOrientation, innerOnly,
                            outputDistanceInM);

  case TYPE_MULTIPOLYGON:
    return straightSkeleton(geom.as<MultiPolygon>(), autoOrientation, innerOnly,
                            outputDistanceInM);

  default:
    return std::make_unique<MultiLineString>();
  }
}

auto
straightSkeleton(const Geometry &geom, bool autoOrientation, bool innerOnly,
                 bool outputDistanceInM, const double & /*toleranceAbs*/)
    -> std::unique_ptr<MultiLineString>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(geom);

  std::unique_ptr<MultiLineString> result(straightSkeleton(
      geom, autoOrientation, NoValidityCheck(), innerOnly, outputDistanceInM));
  propagateValidityFlag(*result, true);
  return result;
}
auto
straightSkeleton(const Polygon &geom, bool /*autoOrientation*/, bool innerOnly,
                 bool outputDistanceInM, const double &toleranceAbs)
    -> std::unique_ptr<MultiLineString>
{
  std::unique_ptr<MultiLineString> result(new MultiLineString);

  if (geom.isEmpty()) {
    return result;
  }

  Kernel::Vector_2                      trans;
  Polygon_with_holes_2 const            polygon  = preparePolygon(geom, trans);
  SHARED_PTR<Straight_skeleton_2> const skeleton = straightSkeleton(polygon);

  if (skeleton == nullptr) {
    BOOST_THROW_EXCEPTION(Exception("CGAL failed to create straightSkeleton"));
  }

  if (outputDistanceInM) {
    straightSkeletonToMultiLineString<Kernel, true>(
        *skeleton, *result, innerOnly, trans, toleranceAbs);
  } else {
    straightSkeletonToMultiLineString<Kernel, false>(
        *skeleton, *result, innerOnly, trans, toleranceAbs);
  }
  return result;
}

auto
straightSkeleton(const MultiPolygon &geom, bool /*autoOrientation*/,
                 bool innerOnly, bool outputDistanceInM,
                 const double &toleranceAbs) -> std::unique_ptr<MultiLineString>
{
  std::unique_ptr<MultiLineString> result(new MultiLineString);

  for (size_t i = 0; i < geom.numGeometries(); i++) {
    Kernel::Vector_2           trans;
    Polygon_with_holes_2 const polygon =
        preparePolygon(geom.polygonN(i), trans);
    SHARED_PTR<Straight_skeleton_2> const skeleton = straightSkeleton(polygon);

    if (skeleton == nullptr) {
      BOOST_THROW_EXCEPTION(
          Exception("CGAL failed to create straightSkeleton"));
    }

    if (outputDistanceInM) {
      straightSkeletonToMultiLineString<Kernel, true>(
          *skeleton, *result, innerOnly, trans, toleranceAbs);
    } else {
      straightSkeletonToMultiLineString<Kernel, false>(
          *skeleton, *result, innerOnly, trans, toleranceAbs);
    }
  }

  return result;
}

auto
approximateMedialAxis(const Geometry &geom, bool projectToEdges)
    -> std::unique_ptr<MultiLineString>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(geom);

  std::unique_ptr<MultiLineString> mx(new MultiLineString);

  std::vector<Polygon> polys;
  extractPolygons(geom, polys);

  for (auto &poly : polys) {
    Kernel::Vector_2                      trans;
    Polygon_with_holes_2 const            polygon = preparePolygon(poly, trans);
    SHARED_PTR<Straight_skeleton_2> const skeleton = straightSkeleton(polygon);

    if (skeleton == nullptr) {
      BOOST_THROW_EXCEPTION(
          Exception("CGAL failed to create straightSkeleton"));
    }

    straightSkeletonToMedialAxis(*skeleton, *mx, trans, projectToEdges);
  }

  propagateValidityFlag(*mx, true);
  return mx;
}

/// @private
auto
extrudeStraightSkeleton(const Polygon &geom, double height)
    -> std::unique_ptr<PolyhedralSurface>
{
  std::unique_ptr<PolyhedralSurface> polys(new PolyhedralSurface);
  if (geom.isEmpty()) {
    return polys;
  }
  Surface_mesh_3 sm;
  CGAL::extrude_skeleton(geom.toPolygon_with_holes_2(), sm,
                         CGAL::parameters::maximum_height(height));
  polys = std::make_unique<PolyhedralSurface>(sm);
  return polys;
}

/// @private
auto
extrudeStraightSkeleton(const Geometry &geom, double height)
    -> std::unique_ptr<PolyhedralSurface>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(geom);

  if (geom.geometryTypeId() != TYPE_POLYGON) {
    BOOST_THROW_EXCEPTION(Exception("Geometry must be a Polygon"));
  }
  std::unique_ptr<PolyhedralSurface> result(
      extrudeStraightSkeleton(geom.as<Polygon>(), height));
  propagateValidityFlag(*result, true);
  return result;
}

/// @private
auto
extrudeStraightSkeleton(const Geometry &geom, double building_height,
                        double roof_height)
    -> std::unique_ptr<PolyhedralSurface>
{
  std::unique_ptr<PolyhedralSurface> result(new PolyhedralSurface);

  if (geom.isEmpty()) {
    return result;
  }

  // Create complete roof with base
  auto completeRoof = extrudeStraightSkeleton(geom, roof_height);

  // Create new roof surface
  auto roof = std::make_unique<PolyhedralSurface>();

  // Copy all non-base faces (roof slopes)
  std::copy_if(completeRoof->begin(), completeRoof->end(),
               std::back_inserter(*roof), detail::isNotBaseFace);

  // Translate roof to building height
  translate(*roof, 0.0, 0.0, building_height);

  // Create building walls
  auto building = extrude(geom.as<Polygon>(), building_height);

  // Create result from building exterior shell
  result = std::make_unique<PolyhedralSurface>(
      building->as<Solid>().exteriorShell());

  // Add filtered roof patches
  result->addPatchs(*roof);

  propagateValidityFlag(*result, true);

  return result;
}

auto
straightSkeletonPartition(const Geometry &geom, bool autoOrientation)
    -> std::unique_ptr<PolyhedralSurface>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(geom);

  std::unique_ptr<PolyhedralSurface> result(new PolyhedralSurface);

  switch (geom.geometryTypeId()) {
  case TYPE_TRIANGLE:
    return straightSkeletonPartition(geom.as<Triangle>().toPolygon(),
                                     autoOrientation);
  case TYPE_POLYGON:
    return straightSkeletonPartition(geom.as<Polygon>(), autoOrientation);
  case TYPE_MULTIPOLYGON:
    return straightSkeletonPartition(geom.as<MultiPolygon>(), autoOrientation);
  default:
    BOOST_THROW_EXCEPTION(
        Exception("Geometry must be a Polygon or MultiPolygon"));
  }

  return result;
}

auto
straightSkeletonPartition(const MultiPolygon &geom, bool autoOrientation)
    -> std::unique_ptr<PolyhedralSurface>
{
  std::unique_ptr<PolyhedralSurface> result(new PolyhedralSurface);

  for (size_t i = 0; i < geom.numGeometries(); i++) {
    std::unique_ptr<PolyhedralSurface> partitioned =
        straightSkeletonPartition(geom.polygonN(i), autoOrientation);
    for (size_t j = 0; j < partitioned->numPatches(); j++) {
      result->addPatch(partitioned->patchN(j));
    }
  }

  return result;
}

auto
straightSkeletonPartition(const Polygon &geom, bool /*autoOrientation*/)
    -> std::unique_ptr<PolyhedralSurface>
{
  std::unique_ptr<PolyhedralSurface> result(new PolyhedralSurface);

  if (geom.isEmpty()) {
    return result;
  }

  Kernel::Vector_2                      trans;
  Polygon_with_holes_2 const            polygon  = preparePolygon(geom, trans);
  SHARED_PTR<Straight_skeleton_2> const skeleton = straightSkeleton(polygon);

  if (skeleton == nullptr) {
    BOOST_THROW_EXCEPTION(Exception("CGAL failed to create straightSkeleton"));
  }

  // Function to create a polygon from a face
  auto create_polygon_from_face =
      [&trans](const Straight_skeleton_2::Face_const_handle &face)
      -> SFCGAL::Polygon {
    std::vector<Point>                         points;
    Straight_skeleton_2::Halfedge_const_handle start   = face->halfedge();
    Straight_skeleton_2::Halfedge_const_handle current = start;
    do {
      const Point_2 &point = current->vertex()->point();
      points.emplace_back(CGAL::to_double(point.x()),
                          CGAL::to_double(point.y()));
      current = current->next();
    } while (current != start);

    SFCGAL::Polygon poly((SFCGAL::LineString(points)));
    algorithm::translate(poly, trans);
    return poly;
  };

  // Function to check if a face corresponds to a hole
  auto is_hole_face =
      [&polygon](const Straight_skeleton_2::Face_const_handle &face) -> bool {
    Point_2 test_point = face->halfedge()->vertex()->point();
    for (auto hit = polygon.holes_begin(); hit != polygon.holes_end(); ++hit) {
      if (hit->bounded_side(test_point) == CGAL::ON_BOUNDED_SIDE) {
        return true;
      }
    }
    return false;
  };

  // Iterate through the faces of the skeleton
  for (auto face = skeleton->faces_begin(); face != skeleton->faces_end();
       ++face) {
    // Skip the faces that correspond to holes
    if (is_hole_face(face)) {
      continue;
    }

    result->addPatch(create_polygon_from_face(face));
  }

  return result;
}
} // namespace SFCGAL::algorithm
