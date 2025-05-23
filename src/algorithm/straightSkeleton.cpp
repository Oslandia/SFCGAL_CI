// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
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
straightSkeletonToMultiLineString(const CGAL::Straight_skeleton_2<K> &ss,
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

  for (Halfedge_const_iterator it = ss.halfedges_begin();
       it != ss.halfedges_end(); ++it) {
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

    LineString *ls = nullptr;
    Point       pa(it->opposite()->vertex()->point());
    Point       pb(it->vertex()->point());
    // avoid degenerate cases.https://gitlab.com/Oslandia/SFCGAL/-/issues/143
    if (pa != pb && distancePointPoint(pa, pb) > toleranceAbs) {
      if (outputDistanceInM) {
        pa.setM(CGAL::to_double(it->opposite()->vertex()->time()));
        pb.setM(CGAL::to_double(it->vertex()->time()));
        ls = new LineString(pa, pb);
      } else {
        ls = new LineString(pa, pb);
      }

      algorithm::translate(*ls, translate);
      result.addGeometry(ls);
    }
  }
}

/**
 * Return abc angle in radians
 */
auto
angle(const Point &a, const Point &b, const Point &c) -> double
{
  Point const ab(CGAL::to_double(b.x() - a.x()),
                 CGAL::to_double(b.y() - a.y()));
  Point const cb(CGAL::to_double(b.x() - c.x()),
                 CGAL::to_double(b.y() - c.y()));

  double const dot =
      (CGAL::to_double(ab.x() * cb.x() + ab.y() * cb.y())); /* dot product */
  double const cross =
      (CGAL::to_double(ab.x() * cb.y() - ab.y() * cb.x())); /* cross product */

  double const alpha = std::atan2(cross, dot);

  return alpha;
}

template <class K>
void
straightSkeletonToMedialAxis(const CGAL::Straight_skeleton_2<K> &ss,
                             MultiLineString                    &result,
                             Kernel::Vector_2                   &translate)
{
  using Ss = CGAL::Straight_skeleton_2<K>;

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

  for (Halfedge_const_iterator it = ss.halfedges_begin();
       it != ss.halfedges_end(); ++it) {
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
      const Point &p1 = it->vertex()->point();
      const Point &p2 = de1->vertex()->point();
      const Point &p3 = de1->opposite()->vertex()->point();

      double const ang = angle(p1, p2, p3);

      if (ang > maxTouchingAngle) {
        continue;
      }
    }

    std::unique_ptr<LineString> ls(
        new LineString(Point(it->opposite()->vertex()->point()),
                       Point(it->vertex()->point())));
    algorithm::translate(*ls, translate);
    result.addGeometry(ls.release());
  }
}

auto
straightSkeleton(const Polygon_with_holes_2 &poly)
    -> SHARED_PTR<Straight_skeleton_2>
{
  SHARED_PTR<CGAL::Straight_skeleton_2<CGAL::Epick>> const sk =
      CGAL::create_interior_straight_skeleton_2(
          poly.outer_boundary().vertices_begin(),
          poly.outer_boundary().vertices_end(), poly.holes_begin(),
          poly.holes_end(), CGAL::Epick());
  SHARED_PTR<Straight_skeleton_2> ret;
  if (sk) {
    ret = CGAL::convert_straight_skeleton_2<Straight_skeleton_2>(*sk);
  }
  return ret;
}

// Throw an exception if any two polygon rings touch,
// since CGAL segfaults in that case.
// @todo find upstream reference to find out exact case
// See https://github.com/Oslandia/SFCGAL/issues/75
void
checkNoTouchingHoles(const Polygon &g)
{
  const size_t numRings = g.numRings();

  for (size_t ri = 0; ri < numRings - 1; ++ri) {
    for (size_t rj = ri + 1; rj < numRings; ++rj) {
      std::unique_ptr<Geometry> inter =
          g.is3D() ? intersection3D(g.ringN(ri), g.ringN(rj))
                   : intersection(g.ringN(ri), g.ringN(rj));

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
  std::unique_ptr<Polygon> cloned(poly.clone());
  algorithm::translate(*cloned, trans);
  Polygon_with_holes_2 ret = cloned->toPolygon_with_holes_2();
  trans                    = -trans;

  return ret;
}

void
extractPolygons(const Geometry &g, std::vector<Polygon> &vect)
{
  switch (g.geometryTypeId()) {
  case TYPE_TRIANGLE:
    vect.push_back(g.as<Triangle>().toPolygon());
    break;
  case TYPE_POLYGON:
    vect.push_back(g.as<Polygon>());
    break;
  case TYPE_MULTIPOLYGON: {
    const auto &mp = g.as<MultiPolygon>();
    for (size_t i = 0; i < mp.numGeometries(); i++) {
      vect.push_back(mp.polygonN(i));
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
straightSkeleton(const Geometry &g, bool          autoOrientation,
                 NoValidityCheck /*unused*/, bool innerOnly,
                 bool outputDistanceInM, const double & /*toleranceAbs*/)
    -> std::unique_ptr<MultiLineString>
{
  switch (g.geometryTypeId()) {
  case TYPE_TRIANGLE:
    return straightSkeleton(g.as<Triangle>().toPolygon(), autoOrientation,
                            innerOnly, outputDistanceInM);

  case TYPE_POLYGON:
    return straightSkeleton(g.as<Polygon>(), autoOrientation, innerOnly,
                            outputDistanceInM);

  case TYPE_MULTIPOLYGON:
    return straightSkeleton(g.as<MultiPolygon>(), autoOrientation, innerOnly,
                            outputDistanceInM);

  default:
    return std::make_unique<MultiLineString>();
  }
}

auto
straightSkeleton(const Geometry &g, bool autoOrientation, bool innerOnly,
                 bool outputDistanceInM, const double & /*toleranceAbs*/)
    -> std::unique_ptr<MultiLineString>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(g);

  std::unique_ptr<MultiLineString> result(straightSkeleton(
      g, autoOrientation, NoValidityCheck(), innerOnly, outputDistanceInM));
  propagateValidityFlag(*result, true);
  return result;
}
auto
straightSkeleton(const Polygon &g, bool /*autoOrientation*/, bool innerOnly,
                 bool outputDistanceInM, const double &toleranceAbs)
    -> std::unique_ptr<MultiLineString>
{
  std::unique_ptr<MultiLineString> result(new MultiLineString);

  if (g.isEmpty()) {
    return result;
  }

  Kernel::Vector_2                      trans;
  Polygon_with_holes_2 const            polygon  = preparePolygon(g, trans);
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
straightSkeleton(const MultiPolygon &g, bool /*autoOrientation*/,
                 bool innerOnly, bool outputDistanceInM,
                 const double &toleranceAbs) -> std::unique_ptr<MultiLineString>
{
  std::unique_ptr<MultiLineString> result(new MultiLineString);

  for (size_t i = 0; i < g.numGeometries(); i++) {
    Kernel::Vector_2           trans;
    Polygon_with_holes_2 const polygon = preparePolygon(g.polygonN(i), trans);
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
approximateMedialAxis(const Geometry &g) -> std::unique_ptr<MultiLineString>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(g);

  std::unique_ptr<MultiLineString> mx(new MultiLineString);

  std::vector<Polygon> polys;
  extractPolygons(g, polys);

  for (auto &poly : polys) {
    Kernel::Vector_2                      trans;
    Polygon_with_holes_2 const            polygon = preparePolygon(poly, trans);
    SHARED_PTR<Straight_skeleton_2> const skeleton = straightSkeleton(polygon);

    if (skeleton == nullptr) {
      BOOST_THROW_EXCEPTION(
          Exception("CGAL failed to create straightSkeleton"));
    }

    straightSkeletonToMedialAxis(*skeleton, *mx, trans);
  }

  propagateValidityFlag(*mx, true);
  return mx;
}

auto
extrudeStraightSkeleton(const Polygon &g, double height)
    -> std::unique_ptr<PolyhedralSurface>
{
  Surface_mesh_3 sm;
  CGAL::extrude_skeleton(g.toPolygon_with_holes_2(), sm,
                         CGAL::parameters::maximum_height(height));
  std::unique_ptr<PolyhedralSurface> polys(new PolyhedralSurface(sm));

  return polys;
}

auto
extrudeStraightSkeleton(const Geometry &g, double height)
    -> std::unique_ptr<PolyhedralSurface>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(g);

  if (g.geometryTypeId() != TYPE_POLYGON) {
    BOOST_THROW_EXCEPTION(Exception("Geometry must be a Polygon"));
  }
  std::unique_ptr<PolyhedralSurface> result(
      extrudeStraightSkeleton(g.as<Polygon>(), height));
  propagateValidityFlag(*result, true);
  return result;
}

auto
extrudeStraightSkeleton(const Geometry &geom, double building_height,
                        double roof_height)
    -> std::unique_ptr<PolyhedralSurface>
{
  // Create complete roof with base
  auto completeRoof = extrudeStraightSkeleton(geom, roof_height);

  // Create new roof surface
  auto roof = std::make_unique<PolyhedralSurface>();

  // Predicate to identify non-base faces (roof slopes)
  auto isNotBaseFace = [](const Polygon &patch) {
    const LineString &exterior = patch.exteriorRing();

    // Check if any point has z != 0 (not a base face)
    return std::any_of(exterior.begin(), exterior.end(),
                       [](const Point &point) { return point.z() != 0.0; });
  };

  std::copy_if(completeRoof->begin(), completeRoof->end(),
               std::back_inserter(*roof), isNotBaseFace);

  // Translate roof to building height
  translate(*roof, 0.0, 0.0, building_height);

  // Create building walls
  auto building = extrude(geom.as<Polygon>(), building_height);

  // Create result from building exterior shell
  auto result = std::make_unique<PolyhedralSurface>(
      building->as<Solid>().exteriorShell());

  // Add filtered roof patches
  result->addPatchs(*roof);

  propagateValidityFlag(*result, true);

  return result;
}

auto
straightSkeletonPartition(const Geometry &g, bool autoOrientation)
    -> std::unique_ptr<PolyhedralSurface>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(g);

  std::unique_ptr<PolyhedralSurface> result(new PolyhedralSurface);

  switch (g.geometryTypeId()) {
  case TYPE_TRIANGLE:
    return straightSkeletonPartition(g.as<Triangle>().toPolygon(),
                                     autoOrientation);
  case TYPE_POLYGON:
    return straightSkeletonPartition(g.as<Polygon>(), autoOrientation);
  case TYPE_MULTIPOLYGON:
    return straightSkeletonPartition(g.as<MultiPolygon>(), autoOrientation);
  default:
    BOOST_THROW_EXCEPTION(
        Exception("Geometry must be a Polygon or MultiPolygon"));
  }

  return result;
}

auto
straightSkeletonPartition(const MultiPolygon &g, bool autoOrientation)
    -> std::unique_ptr<PolyhedralSurface>
{
  std::unique_ptr<PolyhedralSurface> result(new PolyhedralSurface);

  for (size_t i = 0; i < g.numGeometries(); i++) {
    std::unique_ptr<PolyhedralSurface> partitioned =
        straightSkeletonPartition(g.polygonN(i), autoOrientation);
    for (size_t j = 0; j < partitioned->numPatches(); j++) {
      result->addPatch(partitioned->patchN(j));
    }
  }

  return result;
}

auto
straightSkeletonPartition(const Polygon &g, bool /*autoOrientation*/)
    -> std::unique_ptr<PolyhedralSurface>
{
  std::unique_ptr<PolyhedralSurface> result(new PolyhedralSurface);

  if (g.isEmpty()) {
    return result;
  }

  Kernel::Vector_2                      trans;
  Polygon_with_holes_2 const            polygon  = preparePolygon(g, trans);
  SHARED_PTR<Straight_skeleton_2> const skeleton = straightSkeleton(polygon);

  if (skeleton == nullptr) {
    BOOST_THROW_EXCEPTION(Exception("CGAL failed to create straightSkeleton"));
  }

  // Function to create a polygon from a face
  auto create_polygon_from_face =
      [&trans](const Straight_skeleton_2::Face_const_handle &face) {
        std::vector<Point>                         points;
        Straight_skeleton_2::Halfedge_const_handle start   = face->halfedge();
        Straight_skeleton_2::Halfedge_const_handle current = start;
        do {
          const Point_2 &p = current->vertex()->point();
          points.emplace_back(CGAL::to_double(p.x()), CGAL::to_double(p.y()));
          current = current->next();
        } while (current != start);

        SFCGAL::Polygon poly((SFCGAL::LineString(points)));
        algorithm::translate(poly, trans);
        return poly;
      };

  // Function to check if a face corresponds to a hole
  auto is_hole_face = [&polygon](
                          const Straight_skeleton_2::Face_const_handle &face) {
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
