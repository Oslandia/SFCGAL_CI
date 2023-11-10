// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/algorithm/convexHull.h>

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>

#include <SFCGAL/detail/GetPointsVisitor.h>

#include <SFCGAL/Exception.h>
#include <SFCGAL/Kernel.h>
#include <boost/format.hpp>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/point_generators_3.h>
#include <vector>

namespace SFCGAL::algorithm {

using Point_3      = CGAL::Point_3<Kernel>;
using Segment_3    = CGAL::Segment_3<Kernel>;
using Triangle_3   = CGAL::Triangle_3<Kernel>;
using Polyhedron_3 = CGAL::Polyhedron_3<Kernel>;

using Point_2 = CGAL::Point_2<Kernel>;

///
///
///
auto
convexHull(const Geometry &g) -> std::unique_ptr<Geometry>
{
  using CGAL::object_cast;

  if (g.isEmpty()) {
    return std::unique_ptr<Geometry>(g.clone());
  }

  SFCGAL::detail::GetPointsVisitor getPointVisitor;
  const_cast<Geometry &>(g).accept(getPointVisitor);

  // collect points

  if (getPointVisitor.points.empty()) {
    return std::unique_ptr<Geometry>(new GeometryCollection());
  }

  std::vector<Point_2> points;

  points.reserve(getPointVisitor.points.size());
for (auto &point : getPointVisitor.points) {
    points.push_back(point->toPoint_2());
  }

  // resulting extreme points
  std::list<Point_2> epoints;
  CGAL::convex_hull_2(points.begin(), points.end(),
                      std::back_inserter(epoints));

  if (epoints.size() == 1) {
    return std::unique_ptr<Geometry>(new Point(*epoints.begin()));
  } if (epoints.size() == 2) {
    std::list<Point_2>::const_iterator it = epoints.begin();
    return std::unique_ptr<Geometry>(
        new LineString(Point(*it++), Point(*it++)));
  }
  // GEOS does not seem to return triangles
  else if (epoints.size() == 3) {
    std::list<Point_2>::const_iterator it = epoints.begin();
    Point_2                            p(*it++);
    Point_2                            q(*it++);
    Point_2                            r(*it++);
    return std::unique_ptr<Geometry>(new Triangle(p, q, r));
  } else if (epoints.size() > 3) {
    auto *poly = new Polygon;

    for (std::list<Point_2>::const_iterator it = epoints.begin();
         it != epoints.end(); ++it) {
      poly->exteriorRing().addPoint(*it);
    }

    // add back the first point to close the ring
    poly->exteriorRing().addPoint(*epoints.begin());
    return std::unique_ptr<Geometry>(poly);
  } else {
    BOOST_THROW_EXCEPTION(
        Exception("unexpected CGAL output type in CGAL::convex_hull_2"));
  }
}

///
///
///
auto
convexHull3D(const Geometry &g) -> std::unique_ptr<Geometry>
{
  using CGAL::object_cast;

  SFCGAL::detail::GetPointsVisitor getPointVisitor;
  const_cast<Geometry &>(g).accept(getPointVisitor);

  // collect points

  std::vector<Point_3> points;

  points.reserve(getPointVisitor.points.size());
for (auto &point : getPointVisitor.points) {
    points.push_back(point->toPoint_3());
  }

  /*
   * http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Convex_hull_3/Chapter_main.html
   *
   * handles all degenerate cases and returns a CGAL::Object,
   * which may be a point, a segment, a triangle, or a polyhedron.
   */
  CGAL::Object hull;
  CGAL::convex_hull_3(points.begin(), points.end(), hull);

  if (hull.empty()) {
    return std::unique_ptr<Geometry>(new GeometryCollection());
  } if (const auto *point = object_cast<Point_3>(&hull)) {
    return std::unique_ptr<Geometry>(new Point(*point));
  } else if (const auto *segment = object_cast<Segment_3>(&hull)) {
    return std::unique_ptr<Geometry>(
        new LineString(Point(segment->start()), Point(segment->end())));
  } else if (const auto *triangle = object_cast<Triangle_3>(&hull)) {
    return std::unique_ptr<Geometry>(new Triangle(Point(triangle->vertex(0)),
                                                  Point(triangle->vertex(1)),
                                                  Point(triangle->vertex(2))));
  } else if (const auto *polyhedron = object_cast<Polyhedron_3>(&hull)) {
    std::unique_ptr<PolyhedralSurface> result(new PolyhedralSurface());

    for (Polyhedron_3::Facet_const_iterator it_facet =
             polyhedron->facets_begin();
         it_facet != polyhedron->facets_end(); ++it_facet) {
      Polyhedron_3::Halfedge_around_facet_const_circulator it =
          it_facet->facet_begin();

      std::vector<Point> ring;

      do {
        ring.emplace_back(it->vertex()->point());
      } while (++it != it_facet->facet_begin());

      ring.push_back(ring.front());

      result->addPolygon(Polygon(ring));
    }

    return std::unique_ptr<Geometry>(result.release());
  } else {
    BOOST_THROW_EXCEPTION(
        Exception("unexpected CGAL output type in CGAL::convex_hull_3"));
  }
}

} // namespace SFCGAL
