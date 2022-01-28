// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/Geometry.h>
#include <SFCGAL/Kernel.h>
#include <SFCGAL/algorithm/covers.h>
#include <SFCGAL/algorithm/intersection.h>
#include <SFCGAL/algorithm/intersects.h>
#include <SFCGAL/detail/GeometrySet.h>
#include <SFCGAL/detail/TypeForDimension.h>

#include <CGAL/box_intersection_d.h>

using namespace SFCGAL::detail;

namespace SFCGAL {
namespace algorithm {

auto
covers(const PrimitiveHandle<3> &, const PrimitiveHandle<3> &) -> bool
{
  return false;
}

template <int Dim>
auto
segmentsLength(const GeometrySet<Dim> &gs) -> double
{
  double result = 0.0;

  for (auto it = gs.segments().begin(); it != gs.segments().end(); ++it) {
    result = result + sqrt(CGAL::to_double(it->primitive().squared_length()));
  }

  return result;
}

auto
solidsVolume(const GeometrySet<3> &gs, bool planarSurface = false) -> double
{
  double result = 0.0;

  for (const auto &it : gs.volumes()) {
    // TODO : we use areas of surfaces here instead of volumes
    const MarkedPolyhedron &poly = it.primitive();

    if (poly.is_closed() && planarSurface) {
      continue;
    }

    if (!poly.is_closed() && !planarSurface) {
      continue;
    }

    BOOST_ASSERT(poly.is_pure_triangle());

    CGAL::Point_3<Kernel> p1, p2, p3;

    for (MarkedPolyhedron::Facet_const_iterator fit = poly.facets_begin();
         fit != poly.facets_end(); ++fit) {
      MarkedPolyhedron::Halfedge_around_facet_const_circulator cit =
          fit->facet_begin();
      p1 = cit->vertex()->point();
      cit++;
      p2 = cit->vertex()->point();
      cit++;
      p3 = cit->vertex()->point();
      CGAL::Triangle_3<Kernel> tri(p1, p2, p3);
      result = result + sqrt(CGAL::to_double(tri.squared_area()));
    }
  }

  return result;
}

auto
surfacesArea(const GeometrySet<2> &gs) -> double
{
  Kernel::FT result = 0.0;

  for (const auto &it : gs.surfaces()) {
    const CGAL::Polygon_with_holes_2<Kernel> &polygon = it.primitive();
    result = result + CGAL::abs(polygon.outer_boundary().area());

    for (auto hit = polygon.holes_begin(); hit != polygon.holes_end(); ++hit) {
      result = result - CGAL::abs(hit->area());
    }
  }

  return CGAL::to_double(result);
}

auto
surfacesArea(const GeometrySet<3> &gs) -> double
{
  double result = 0.0;

  if (gs.surfaces().empty() && !gs.volumes().empty()) {
    result = solidsVolume(gs, /* planarSurface = */ true);
  }

  for (const auto &it : gs.surfaces()) {
    result = result + sqrt(CGAL::to_double(it.primitive().squared_area()));
  }

  return result;
}

auto
solidsVolume(const GeometrySet<2> &) -> double
{
  return 0.0;
}

template <int Dim>
auto
equalLength(const GeometrySet<Dim> &a, const GeometrySet<Dim> &b, int dim)
    -> bool
{
  // compare 'length' of primitives in A with 'length' of primitives in B
  // 'length' is :
  // - number of elements for points
  // - length for segments
  // - area for surfaces
  // - should be volume for volumes. We use area here

  double tol = 1e-9;

  switch (dim) {
  case 0: {
    if (a.points().size() != b.points().size()) {
      return false;
    }
  } break;

  case 1: {

    //
    // Compare lengths
    double lengthA = segmentsLength(a);
    double lengthB = segmentsLength(b);
    double cmp     = (lengthA - lengthB) * (lengthA - lengthB);

    if (cmp > tol) {
      return false;
    }
  } break;

  case 2: {
    //
    // Compare areas
    double areaA = surfacesArea(a);
    double areaB = surfacesArea(b);
    double cmp   = (areaA - areaB) * (areaA - areaB);

    if (cmp > tol) {
      return false;
    }
  } break;

  case 3: {
    // Compare volumes
    double volA = solidsVolume(a);
    double volB = solidsVolume(b);
    double cmp  = (volA - volB) * (volA - volB);

    if (cmp > tol) {
      return false;
    }
  } break;
  }

  return true;
}

template <int Dim>
auto
covers(const GeometrySet<Dim> &a, const GeometrySet<Dim> &b) -> bool
{
  int dimA = a.dimension();
  int dimB = b.dimension();

  if (dimA == -1 || dimB == -1) {
    return false;
  }

  if (dimB > dimA) {
    return false;
  }

  //
  // This is a very naive (not efficient) implementation of covers() !
  //
  // covers(A,B) <=> A inter B == B
  // '==' is here implemented with comparison of length, area and volumes
  // TODO use only predicates if possible
  GeometrySet<Dim> inter;
  algorithm::intersection(a, b, inter);

  if (b.hasPoints() && !equalLength(b, inter, 0)) {
    return false;
  }

  if (b.hasSegments() && !equalLength(b, inter, 1)) {
    return false;
  }

  if (b.hasSurfaces() && !equalLength(b, inter, 2)) {
    return false;
  }

  if (b.hasVolumes() && !equalLength(b, inter, 3)) {
    return false;
  }

  return true;
}

template bool
covers<2>(const GeometrySet<2> &a, const GeometrySet<2> &b);
template bool
covers<3>(const GeometrySet<3> &a, const GeometrySet<3> &b);

auto
covers(const Geometry &ga, const Geometry &gb) -> bool
{
  if (ga.isEmpty() || gb.isEmpty()) {
    return false;
  }

  GeometrySet<2> gsa(ga);
  GeometrySet<2> gsb(gb);

  return covers(gsa, gsb);
}

auto
covers3D(const Geometry &ga, const Geometry &gb) -> bool
{
  if (ga.isEmpty() || gb.isEmpty()) {
    return false;
  }

  GeometrySet<3> gsa(ga);
  GeometrySet<3> gsb(gb);

  return covers(gsa, gsb);
}
} // namespace algorithm
} // namespace SFCGAL
