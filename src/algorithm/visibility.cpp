// Copyright (c) 2023-2023, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/visibility.h"
#include "SFCGAL/Polygon.h"

#include "SFCGAL/algorithm/isValid.h"

#include "SFCGAL/Exception.h"

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Triangular_expansion_visibility_2.h>

#include <CGAL/Arr_naive_point_location.h>
#include <memory>

namespace SFCGAL::algorithm {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

using Traits_2              = CGAL::Arr_segment_traits_2<Kernel>;
using Arrangement_2         = CGAL::Arrangement_2<Traits_2>;
using Face_handle           = Arrangement_2::Face_handle;
using Face_const_handle     = Arrangement_2::Face_const_handle;
using Halfedge_const_handle = Arrangement_2::Halfedge_const_handle;
using TEV =
    CGAL::Triangular_expansion_visibility_2<Arrangement_2, CGAL::Tag_true>;

static auto
query_visibility(Face_handle fh, Halfedge_const_handle he)
    -> std::unique_ptr<Polygon>
{

  std::unique_ptr<LineString> extRing{new LineString()};
  // Make sure the visibility polygon we find has an outer boundary
  if (fh->has_outer_ccb()) {
    Arrangement_2::Ccb_halfedge_circulator curr = fh->outer_ccb();

    // find the right halfedge first
    if (he != Halfedge_const_handle()) {
      while (++curr != fh->outer_ccb()) {
        if (curr->source()->point() == he->source()->point()) {
          break;
        }
      }
    }

    Arrangement_2::Ccb_halfedge_circulator const first = curr;
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

/// @} end of private section

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

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
visibility(const Geometry &polygon, const Geometry &point,
           NoValidityCheck /*unused*/) -> std::unique_ptr<Polygon>
{

  Point_2 const                     queryPoint{point.as<Point>().toPoint_2()};
  std::unique_ptr<LineString> const extRing{new LineString()};

  // insert geometry into the arrangement
  CGAL::Polygon_with_holes_2 pwh{
      polygon.as<Polygon>().toPolygon_with_holes_2(true)};
  Arrangement_2 arr;

  CGAL::insert(arr, pwh.outer_boundary().edges_begin(),
               pwh.outer_boundary().edges_end());
  // the holes
  for (auto hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit) {
    CGAL::insert(arr, hit->edges_begin(), hit->edges_end());
  }

  // Find the face
  CGAL::Arr_naive_point_location<Arrangement_2> const  pl(arr);
  CGAL::Arr_point_location_result<Arrangement_2>::Type obj =
      pl.locate(queryPoint);

  Arrangement_2 output_arr;
  Face_handle   fh;

  // Create Triangular Expansion Visibility object.
  TEV const tev(arr);
#if CGAL_VERSION_MAJOR < 6
  switch (obj.which())
#else
  switch (obj.index())
#endif
  {
  case 0: {
    Halfedge_const_handle he = Halfedge_const_handle();

    // If the point is in a boundary segment, find the corresponding half edge
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
        BOOST_THROW_EXCEPTION(
            Exception("Can not find corresponding half edge (from vertex)."));
      }

      cont = !Segment_2(he->source()->point(), he->target()->point())
                  .has_on(queryPoint) ||
             he->source()->point() == queryPoint || he->face()->is_unbounded();
    }

    // Use the half edge to compute the visibility
    fh = tev.compute_visibility(queryPoint, he, output_arr);
    break;
  }
  case 1: {
    auto *he =
#if CGAL_VERSION_MAJOR < 6
        boost::get<Arrangement_2::Halfedge_const_handle>(&obj);
#else
        std::get_if<Arrangement_2::Halfedge_const_handle>(&obj);
#endif
    if (he != nullptr) {
      fh = tev.compute_visibility(queryPoint, *he, output_arr);
    } else {
      BOOST_THROW_EXCEPTION(Exception("Can not find corresponding hedge."));
    }
    break;
  }
  case 2: {
    auto *face =
#if CGAL_VERSION_MAJOR < 6
        boost::get<Arrangement_2::Face_const_handle>(&obj);
#else
        std::get_if<Arrangement_2::Face_const_handle>(&obj);
#endif
    if ((face != nullptr) && !((*face)->is_unbounded())) {
      fh = tev.compute_visibility(queryPoint, *face, output_arr);
    } else {
      BOOST_THROW_EXCEPTION(Exception("Can not find corresponding face."));
    }
    break;
  }
  default:
    break;
  }

  return query_visibility(fh, fh->outer_ccb());
}
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
           const Geometry &pointB, NoValidityCheck /*unused*/)
    -> std::unique_ptr<Polygon>
{

  Point_2 const startPoint{pointA.as<Point>().toPoint_2()};
  Point_2 const endPoint{pointB.as<Point>().toPoint_2()};
  Point_2 const queryPoint{pointB.as<Point>().toPoint_2()};

  // insert geometry into the arrangement
  CGAL::Polygon_with_holes_2 pwh{
      polygon.as<Polygon>().toPolygon_with_holes_2(true)};
  Arrangement_2 arr;

  CGAL::insert(arr, pwh.outer_boundary().edges_begin(),
               pwh.outer_boundary().edges_end());
  for (auto hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit) {
    CGAL::insert(arr, hit->edges_begin(), hit->edges_end());
  }

  // If the point is in a boundary segment, find the corresponding half edge
  Halfedge_const_handle he = arr.halfedges_begin();
  bool cont = !Segment_2(he->source()->point(), he->target()->point())
                   .has_on(queryPoint) ||
              he->source()->point() == startPoint ||
              he->target()->point() == endPoint || he->face()->is_unbounded();
  // While we are not in the right half edge, or while q is the source,
  // continue
  while (cont) {
    he++;
    if (he == arr.halfedges_end()) {
      BOOST_THROW_EXCEPTION(Exception("Can not find corresponding half edge."));
    }

    cont = !Segment_2(he->source()->point(), he->target()->point())
                .has_on(queryPoint) ||
           he->source()->point() == queryPoint || he->face()->is_unbounded();
  }

  // visibility query
  Arrangement_2     output_arr;
  TEV const         tev(arr);
  Face_handle const fh = tev.compute_visibility(endPoint, he, output_arr);

  return query_visibility(fh, he);
}

} // namespace SFCGAL::algorithm
