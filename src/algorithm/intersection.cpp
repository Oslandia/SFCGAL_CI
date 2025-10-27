// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/algorithm/collect.h"
#include "SFCGAL/algorithm/collectionHomogenize.h"
#include "SFCGAL/algorithm/intersects.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/detail/GeometrySet.h"
#include "SFCGAL/detail/tools/Registry.h"

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/box_intersection_d.h>

//
// Intersection kernel

using namespace SFCGAL::detail;

namespace SFCGAL::algorithm {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

// see Intersection3D.cpp
void
intersection(const PrimitiveHandle<3> &pa, const PrimitiveHandle<3> &pb,
             GeometrySet<3> &output, dim_t<3>);
// see Intersection2D.cpp
void
intersection(const PrimitiveHandle<2> &pa, const PrimitiveHandle<2> &pb,
             GeometrySet<2> &output, dim_t<2>);
//
// We deal here with symmetric call
template <int Dim>
void
dispatch_intersection_sym(const PrimitiveHandle<Dim> &pa,
                          const PrimitiveHandle<Dim> &pb,
                          GeometrySet<Dim>           &output)
{
  // assume types are ordered by dimension within the boost::variant
  if (pa.handle.which() >= pb.handle.which()) {
    intersection(pa, pb, output, dim_t<Dim>());
  } else {
    intersection(pb, pa, output, dim_t<Dim>());
  }
}

template <int Dim>
void
intersection(const PrimitiveHandle<Dim> &primitiveHandle1,
             const PrimitiveHandle<Dim> &primitiveHandle2,
             GeometrySet<Dim>           &output)
{
  dispatch_intersection_sym(primitiveHandle1, primitiveHandle2, output);
}

template void
intersection<2>(const PrimitiveHandle<2> &primitiveHandle1,
                const PrimitiveHandle<2> &primitiveHandle2, GeometrySet<2> &);
template void
intersection<3>(const PrimitiveHandle<3> &primitiveHandle1,
                const PrimitiveHandle<3> &primitiveHandle2, GeometrySet<3> &);

template <int Dim>
struct intersection_cb {
  intersection_cb(GeometrySet<Dim> &out) : output(out) {}

  void
  operator()(const typename PrimitiveBox<Dim>::Type &a,
             const typename PrimitiveBox<Dim>::Type &b)
  {
    dispatch_intersection_sym<Dim>(*a.handle(), *b.handle(), output);
  }

  GeometrySet<Dim> &output;
};

/**
 * intersection post processing
 */
void
post_intersection(const GeometrySet<2> &input, GeometrySet<2> &output)
{
  //
  // reverse orientation of polygons if needed
  for (const auto &it : input.surfaces()) {
    const CGAL::Polygon_with_holes_2<Kernel> &p     = it.primitive();
    CGAL::Polygon_2<Kernel>                   outer = p.outer_boundary();

    if (outer.orientation() == CGAL::CLOCKWISE) {
      outer.reverse_orientation();
    }

    std::list<CGAL::Polygon_2<Kernel>> rings;

    for (auto hit = p.holes_begin(); hit != p.holes_end(); ++hit) {
      rings.push_back(*hit);

      if (hit->orientation() == CGAL::COUNTERCLOCKWISE) {
        rings.back().reverse_orientation();
      }
    }

    output.surfaces().emplace_back(
        CGAL::Polygon_with_holes_2<Kernel>(outer, rings.begin(), rings.end()));
  }

  output.points()   = input.points();
  output.segments() = input.segments();
  output.volumes()  = input.volumes();
}

void
post_intersection(const GeometrySet<3> &input, GeometrySet<3> &output)
{
  // nothing special to do
  output = input;
}
/// @} end of private section

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

template <int Dim>
auto
intersection(const GeometrySet<Dim> &geometrySet1,
             const GeometrySet<Dim> &geometrySet2, GeometrySet<Dim> &output)
    -> void
{
  typename SFCGAL::detail::HandleCollection<Dim>::Type ahandles;
  typename SFCGAL::detail::HandleCollection<Dim>::Type bhandles;
  typename SFCGAL::detail::BoxCollection<Dim>::Type    aboxes;
  typename SFCGAL::detail::BoxCollection<Dim>::Type    bboxes;
  geometrySet1.computeBoundingBoxes(ahandles, aboxes);
  geometrySet2.computeBoundingBoxes(bhandles, bboxes);

  GeometrySet<Dim>           temp;
  GeometrySet<Dim>           temp2;
  intersection_cb<Dim> const cb(temp);
  CGAL::box_intersection_d(aboxes.begin(), aboxes.end(), bboxes.begin(),
                           bboxes.end(), cb);

  post_intersection(temp, temp2);
  output.merge(temp2);
}

/// @private
template void
intersection<2>(const GeometrySet<2> &geometrySet1,
                const GeometrySet<2> &geometrySet2, GeometrySet<2> &);
/// @private
template void
intersection<3>(const GeometrySet<3> &geometrySet1,
                const GeometrySet<3> &geometrySet2, GeometrySet<3> &);

/// @private
auto
intersection(const Geometry &geometry1, const Geometry &geometry2,
             NoValidityCheck /*unused*/) -> std::unique_ptr<Geometry>
{
  GeometrySet<2> const gsa(geometry1);
  GeometrySet<2> const gsb(geometry2);
  GeometrySet<2>       output;
  algorithm::intersection(gsa, gsb, output);

  GeometrySet<2> filtered;
  output.filterCovered(filtered);
  return filtered.recompose();
}

/// @private
auto
intersection(const Geometry &geometry1, const Geometry &geometry2)
    -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(geometry1);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(geometry2);

  return intersection(geometry1, geometry2, NoValidityCheck());
}

/// @private
auto
intersection3D(const Geometry &geometry1, const Geometry &geometry2,
               NoValidityCheck /*unused*/) -> std::unique_ptr<Geometry>
{
  GeometrySet<3> const gsa(geometry1);
  GeometrySet<3> const gsb(geometry2);
  GeometrySet<3>       output;
  algorithm::intersection(gsa, gsb, output);

  GeometrySet<3> filtered;
  output.filterCovered(filtered);

  return filtered.recompose();
}

/// @private
auto
intersection3D(const Geometry &geometry1, const Geometry &geometry2)
    -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(geometry1);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(geometry2);

  return intersection3D(geometry1, geometry2, NoValidityCheck());
}
} // namespace SFCGAL::algorithm
