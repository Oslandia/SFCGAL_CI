// Copyright (c) 2023-2023, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later
#include <SFCGAL/Polygon.h>
#include <SFCGAL/algorithm/visibility.h>

#include <SFCGAL/algorithm/isValid.h>

#include <SFCGAL/Exception.h>

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Triangular_expansion_visibility_2.h>

#include <CGAL/Arr_naive_point_location.h>
#include <memory>

namespace SFCGAL {
namespace algorithm {

using Point_2               = Kernel::Point_2;
using Polygon_2             = CGAL::Polygon_2<Kernel>;
using Segment_2             = Kernel::Segment_2;
using Traits_2              = CGAL::Arr_segment_traits_2<Kernel>;
using Arrangement_2         = CGAL::Arrangement_2<Traits_2>;
using Face_handle           = Arrangement_2::Face_handle;
using Halfedge_const_handle = Arrangement_2::Halfedge_const_handle;
using TEV =
    CGAL::Triangular_expansion_visibility_2<Arrangement_2, CGAL::Tag_true>;
using PolygonWithHoles = CGAL::Polygon_with_holes_2<Kernel>;
///
///
///
auto
visibility(const Geometry &polygon, const Geometry &point)
    -> std::unique_ptr<Polygon>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(polygon);

  std::unique_ptr<Polygon> result(
      visibility(polygon, point, NoValidityCheck()));
  propagateValidityFlag(*result, true);
  return result;
}

auto
visibility(const Geometry &polygon, const Geometry &point, NoValidityCheck)
    -> std::unique_ptr<Polygon>
{

  Point_2                     queryPoint{point.as<Point>().toPoint_2()};
  std::unique_ptr<LineString> extRing{new LineString()};

  // insert geometry into the arrangement
  CGAL::Polygon_with_holes_2 pwh{
      polygon.as<Polygon>().toPolygon_with_holes_2()};
  Arrangement_2 arr;

  CGAL::insert(arr, pwh.outer_boundary().edges_begin(),
               pwh.outer_boundary().edges_end());
  // the holes
  for (PolygonWithHoles::Hole_const_iterator hit = pwh.holes_begin();
       hit != pwh.holes_end(); ++hit)
    CGAL::insert(arr, hit->edges_begin(), hit->edges_end());

  // Find the face
  Arrangement_2::Face_const_handle                    *face;
  CGAL::Arr_naive_point_location<Arrangement_2>        pl(arr);
  CGAL::Arr_point_location_result<Arrangement_2>::Type obj =
      pl.locate(queryPoint);

  // The query point locates in the interior of a face
  face = boost::get<Arrangement_2::Face_const_handle>(&obj);
  Arrangement_2         output_arr;
  Face_handle           fh;
  Halfedge_const_handle he = Halfedge_const_handle();

  // Create Triangular Expansion Visibility object.
  TEV tev(arr);

  // If the point is within a face, we can compute the visibility that way
  if (face != nullptr) {
    fh = tev.compute_visibility(queryPoint, *face, output_arr);
  } else {
    // If the point in a boundary segment, find the corresponding half edge
    he        = arr.halfedges_begin();
    bool cont = !Segment_2(he->source()->point(), he->target()->point())
                     .has_on(queryPoint) ||
                he->source()->point() == queryPoint ||
                he->face()->is_unbounded();
    // While we are not in the right half edge, or while q is the source,
    // continue
    while (cont) {
      he++;
      if (he == arr.halfedges_end()) {
        throw(std::exception());
      }

      cont = !Segment_2(he->source()->point(), he->target()->point())
                  .has_on(queryPoint) ||
             he->source()->point() == queryPoint || he->face()->is_unbounded();
    }

    // Use the half edge to compute the visibility
    fh = tev.compute_visibility(queryPoint, he, output_arr);
  }
  // Make sure the visibility polygon we find has an outer boundary
  if (fh->has_outer_ccb()) {
    Arrangement_2::Ccb_halfedge_circulator curr = fh->outer_ccb();

    // find the right halfedge first
    if (he != Halfedge_const_handle())
      while (++curr != fh->outer_ccb())
        if (curr->source()->point() == he->source()->point())
          break;

    Arrangement_2::Ccb_halfedge_circulator first = curr;
    extRing->addPoint(Point(curr->source()->point()));

    // Save the points from the visibility polygon
    while (++curr != first) {
      extRing->addPoint(Point(curr->source()->point()));
    }
  }

  extRing->closes();
  std::unique_ptr<Polygon> result{new Polygon(extRing.release())};

  return result;
}
///
///
///
auto
visibility(const Geometry &polygon, const Geometry &pointA,
           const Geometry &pointB) -> std::unique_ptr<Polygon>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(polygon);

  std::unique_ptr<Polygon> result(
      visibility(polygon, pointA, pointB, NoValidityCheck()));
  propagateValidityFlag(*result, true);
  return result;
}

auto
visibility(const Geometry &polygon, const Geometry &pointA,
           const Geometry &pointB, NoValidityCheck) -> std::unique_ptr<Polygon>
{

  Point_2 startPoint{pointA.as<Point>().toPoint_2()};
  Point_2 endPoint{pointB.as<Point>().toPoint_2()};

  // insert geometry into the arrangement
  CGAL::Polygon_with_holes_2 pwh{
      polygon.as<Polygon>().toPolygon_with_holes_2()};
  Arrangement_2 arr;

  CGAL::insert(arr, pwh.outer_boundary().edges_begin(),
               pwh.outer_boundary().edges_end());
  for (PolygonWithHoles::Hole_const_iterator hit = pwh.holes_begin();
       hit != pwh.holes_end(); ++hit)
    CGAL::insert(arr, hit->edges_begin(), hit->edges_end());

  // If the point is in a boundary segment, find the corresponding half edge
  Halfedge_const_handle he = arr.halfedges_begin();
  while (he->source()->point() != startPoint ||
         he->target()->point() != endPoint) {
    he++;
    if (he == arr.halfedges_end()) {
      throw(std::exception());
    }
  }

  // visibility query
  Arrangement_2 output_arr;
  TEV           tev(arr);
  Face_handle   fh = tev.compute_visibility(endPoint, he, output_arr);

  // export the visibility region.
  Arrangement_2::Ccb_halfedge_circulator curr = fh->outer_ccb();
  std::unique_ptr<LineString>            extRing{new LineString()};

  // Make sure the visibility polygon we find has an outer boundary
  if (fh->has_outer_ccb()) {
    extRing->addPoint(Point(curr->target()->point()));
    while (++curr != fh->outer_ccb()) {
      extRing->addPoint(Point(curr->target()->point()));
    }
  }
  extRing->closes();
  std::unique_ptr<Polygon> result{new Polygon(extRing.release())};

  return result;
}

} // namespace algorithm
} // namespace SFCGAL
