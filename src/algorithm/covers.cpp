// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/covers.h"
#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/algorithm/intersects.h"
#include "SFCGAL/detail/GeometrySet.h"
#include "SFCGAL/detail/TypeForDimension.h"

#include <CGAL/box_intersection_d.h>

using namespace SFCGAL::detail;

namespace SFCGAL::algorithm {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

auto
covers(const PrimitiveHandle<3> & /*unused*/,
       const PrimitiveHandle<3> & /*unused*/) -> bool
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

    CGAL::Point_3<Kernel> p1;
    CGAL::Point_3<Kernel> p2;
    CGAL::Point_3<Kernel> p3;

    for (MarkedPolyhedron::Facet_const_iterator fit = poly.facets_begin();
         fit != poly.facets_end(); ++fit) {
      MarkedPolyhedron::Halfedge_around_facet_const_circulator cit =
          fit->facet_begin();
      p1 = cit->vertex()->point();
      cit++;
      p2 = cit->vertex()->point();
      cit++;
      p3 = cit->vertex()->point();
      CGAL::Triangle_3<Kernel> const tri(p1, p2, p3);
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
solidsVolume(const GeometrySet<2> & /*unused*/) -> double
{
  return 0.0;
}

template <int Dim>
auto
equalLength(const GeometrySet<Dim> &geometrySet1,
            const GeometrySet<Dim> &geometrySet2, int dim) -> bool
{
  // compare 'length' of primitives in A with 'length' of primitives in B
  // 'length' is :
  // - number of elements for points
  // - length for segments
  // - area for surfaces
  // - should be volume for volumes. We use area here

  double const tol = 1e-9;

  switch (dim) {
  case 0: {
    if (geometrySet1.points().size() != geometrySet2.points().size()) {
      return false;
    }
  } break;

  case 1: {

    //
    // Compare lengths
    double const lengthA = segmentsLength(geometrySet1);
    double const lengthB = segmentsLength(geometrySet2);
    double const cmp     = (lengthA - lengthB) * (lengthA - lengthB);

    if (cmp > tol) {
      return false;
    }
  } break;

  case 2: {
    //
    // Compare areas
    double const areaA = surfacesArea(geometrySet1);
    double const areaB = surfacesArea(geometrySet2);
    double const cmp   = (areaA - areaB) * (areaA - areaB);

    if (cmp > tol) {
      return false;
    }
  } break;

  case 3: {
    // Compare volumes
    double const volA = solidsVolume(geometrySet1);
    double const volB = solidsVolume(geometrySet2);
    double const cmp  = (volA - volB) * (volA - volB);

    if (cmp > tol) {
      return false;
    }
  } break;
  }

  return true;
}

template <int Dim>
auto
covers(const GeometrySet<Dim> &geometrySet1,
       const GeometrySet<Dim> &geometrySet2) -> bool
{
  int const dimA = geometrySet1.dimension();
  int const dimB = geometrySet2.dimension();

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
  algorithm::intersection(geometrySet1, geometrySet2, inter);

  if (geometrySet2.hasPoints() && !equalLength(geometrySet2, inter, 0)) {
    return false;
  }

  if (geometrySet2.hasSegments() && !equalLength(geometrySet2, inter, 1)) {
    return false;
  }

  if (geometrySet2.hasSurfaces() && !equalLength(geometrySet2, inter, 2)) {
    return false;
  }

  if (geometrySet2.hasVolumes() && !equalLength(geometrySet2, inter, 3)) {
    return false;
  }

  return true;
}

template bool
covers<2>(const GeometrySet<2> &geometrySet1,
          const GeometrySet<2> &geometrySet2);
template bool
covers<3>(const GeometrySet<3> &geometrySet1,
          const GeometrySet<3> &geometrySet2);

/// @} end of private section

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

/// @private
auto
covers(const Geometry &geometry1, const Geometry &geometry2) -> bool
{
  if (geometry1.isEmpty() || geometry2.isEmpty()) {
    return false;
  }

  GeometrySet<2> const gsa(geometry1);
  GeometrySet<2> const gsb(geometry2);

  return covers(gsa, gsb);
}

/// @private
auto
covers3D(const Geometry &geometry1, const Geometry &geometry2) -> bool
{
  if (geometry1.isEmpty() || geometry2.isEmpty()) {
    return false;
  }

  GeometrySet<3> const gsa(geometry1);
  GeometrySet<3> const gsb(geometry2);

  return covers(gsa, gsb);
}
} // namespace SFCGAL::algorithm
