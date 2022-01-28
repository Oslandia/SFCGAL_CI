// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/algorithm/distance.h>
#include <SFCGAL/algorithm/straightSkeleton.h>

#include <SFCGAL/Envelope.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/Triangle.h>

#include <SFCGAL/Exception.h>

#include <SFCGAL/algorithm/intersection.h>
#include <SFCGAL/algorithm/isValid.h>
#include <SFCGAL/algorithm/orientation.h>
#include <SFCGAL/algorithm/translate.h>

#include <CGAL/Straight_skeleton_converter_2.h>
#include <CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>

#include <memory>

namespace SFCGAL {
namespace algorithm {

using Point_2              = Kernel::Point_2;
using Polygon_2            = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;
using Straight_skeleton_2  = CGAL::Straight_skeleton_2<Kernel>;

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

  Halfedge_const_handle null_halfedge;
  Vertex_const_handle   null_vertex;

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
static auto
angle(const Point &a, const Point &b, const Point &c) -> double
{
  Point ab(to_double(b.x() - a.x()), to_double(b.y() - a.y()));
  Point cb(to_double(b.x() - c.x()), to_double(b.y() - c.y()));

  double dot = (to_double(ab.x() * cb.x() + ab.y() * cb.y())); /* dot product */
  double cross =
      (to_double(ab.x() * cb.y() - ab.y() * cb.x())); /* cross product */

  double alpha = std::atan2(cross, dot);

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

      double ang = angle(p1, p2, p3);

      if (ang > maxTouchingAngle)
        continue;
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
    -> boost::shared_ptr<Straight_skeleton_2>
{
  boost::shared_ptr<CGAL::Straight_skeleton_2<CGAL::Epick>> sk =
      CGAL::create_interior_straight_skeleton_2(
          poly.outer_boundary().vertices_begin(),
          poly.outer_boundary().vertices_end(), poly.holes_begin(),
          poly.holes_end(), CGAL::Epick());
  boost::shared_ptr<Straight_skeleton_2> ret;
  if (sk)
    ret = CGAL::convert_straight_skeleton_2<Straight_skeleton_2>(*sk);
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
        BOOST_THROW_EXCEPTION(
            NotImplementedException("straight skeleton of Polygon with point "
                                    "touching rings is not implemented."));
      }
    }
  }
}

auto
preparePolygon(const Polygon &poly, Kernel::Vector_2 &trans)
    -> Polygon_with_holes_2
{
  checkNoTouchingHoles(poly);
  Envelope env = poly.envelope();
  trans        = Kernel::Vector_2(-env.xMin(), -env.yMin());

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

///
///
///
auto
straightSkeleton(const Geometry &g, bool autoOrientation, NoValidityCheck,
                 bool innerOnly, bool outputDistanceInM,
                 const double &toleranceAbs) -> std::unique_ptr<MultiLineString>
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
                 bool outputDistanceInM, const double &toleranceAbs)
    -> std::unique_ptr<MultiLineString>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(g);

  std::unique_ptr<MultiLineString> result(straightSkeleton(
      g, autoOrientation, NoValidityCheck(), innerOnly, outputDistanceInM));
  propagateValidityFlag(*result, true);
  return result;
}
///
///
///
auto
straightSkeleton(const Polygon &g, bool /*autoOrientation*/, bool innerOnly,
                 bool outputDistanceInM, const double &toleranceAbs)
    -> std::unique_ptr<MultiLineString>
{
  std::unique_ptr<MultiLineString> result(new MultiLineString);

  if (g.isEmpty()) {
    return result;
  }

  Kernel::Vector_2                       trans;
  Polygon_with_holes_2                   polygon  = preparePolygon(g, trans);
  boost::shared_ptr<Straight_skeleton_2> skeleton = straightSkeleton(polygon);

  if (!skeleton.get()) {
    BOOST_THROW_EXCEPTION(Exception("CGAL failed to create straightSkeleton"));
  }

  if (outputDistanceInM)
    straightSkeletonToMultiLineString<Kernel, true>(
        *skeleton, *result, innerOnly, trans, toleranceAbs);
  else
    straightSkeletonToMultiLineString<Kernel, false>(
        *skeleton, *result, innerOnly, trans, toleranceAbs);
  return result;
}

///
///
///
auto
straightSkeleton(const MultiPolygon &g, bool /*autoOrientation*/,
                 bool innerOnly, bool outputDistanceInM,
                 const double &toleranceAbs) -> std::unique_ptr<MultiLineString>
{
  std::unique_ptr<MultiLineString> result(new MultiLineString);

  for (size_t i = 0; i < g.numGeometries(); i++) {
    Kernel::Vector_2     trans;
    Polygon_with_holes_2 polygon = preparePolygon(g.polygonN(i), trans);
    boost::shared_ptr<Straight_skeleton_2> skeleton = straightSkeleton(polygon);

    if (!skeleton.get()) {
      BOOST_THROW_EXCEPTION(
          Exception("CGAL failed to create straightSkeleton"));
    }

    if (outputDistanceInM)
      straightSkeletonToMultiLineString<Kernel, true>(
          *skeleton, *result, innerOnly, trans, toleranceAbs);
    else
      straightSkeletonToMultiLineString<Kernel, false>(
          *skeleton, *result, innerOnly, trans, toleranceAbs);
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
    Kernel::Vector_2     trans;
    Polygon_with_holes_2 polygon = preparePolygon(poly, trans);
    boost::shared_ptr<Straight_skeleton_2> skeleton = straightSkeleton(polygon);

    if (!skeleton.get()) {
      BOOST_THROW_EXCEPTION(
          Exception("CGAL failed to create straightSkeleton"));
    }

    straightSkeletonToMedialAxis(*skeleton, *mx, trans);
  }

  propagateValidityFlag(*mx, true);
  return mx;
}

} // namespace algorithm
} // namespace SFCGAL
