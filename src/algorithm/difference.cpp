// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/Exception.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/TriangulatedSurface.h>
#include <SFCGAL/algorithm/difference.h>
#include <SFCGAL/algorithm/differencePrimitives.h>
#include <SFCGAL/algorithm/isValid.h>
#include <SFCGAL/detail/GeometrySet.h>
#include <SFCGAL/triangulate/triangulatePolygon.h>

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/box_intersection_d.h>

//
// Intersection kernel

using namespace SFCGAL::detail;

namespace SFCGAL::algorithm {

template <int Dim>
struct CollisionMapper {
  using PrimitiveHandleSet = std::vector<PrimitiveHandle<Dim> *>;
  using Map = std::map<PrimitiveHandle<Dim> *, PrimitiveHandleSet>;
  CollisionMapper(Map &map) : _map(map){};
  void
  operator()(const typename PrimitiveBox<Dim>::Type &a,
             const typename PrimitiveBox<Dim>::Type &b)
  {
    _map[a.handle()].push_back(b.handle());
  }

private:
  Map &_map;
};

template <typename OutputIteratorType>
auto
difference(const Point_2 &primitive, const PrimitiveHandle<2> &pb,
           OutputIteratorType out) -> OutputIteratorType
{
  switch (pb.handle.which()) {
  case PrimitivePoint:
    difference(primitive, *pb.as<Point_2>(), out);
    break;

  case PrimitiveSegment:
    difference(primitive, *pb.as<Segment_2>(), out);
    break;

  case PrimitiveSurface:
    difference(primitive, *pb.as<PolygonWH_2>(), out);
    break;
  }

  return out;
}

template <typename OutputIteratorType>
auto
difference(const Segment_2 &primitive, const PrimitiveHandle<2> &pb,
           OutputIteratorType out) -> OutputIteratorType
{
  switch (pb.handle.which()) {
  case PrimitivePoint:
    *out++ = primitive;
    break;

  case PrimitiveSegment:
    difference(primitive, *pb.as<Segment_2>(), out);
    break;

  case PrimitiveSurface:
    difference(primitive, *pb.as<PolygonWH_2>(), out);
    break;
  }

  return out;
}

template <typename OutputIteratorType>
auto
difference(const PolygonWH_2 &primitive, const PrimitiveHandle<2> &pb,
           OutputIteratorType out) -> OutputIteratorType
{
  switch (pb.handle.which()) {
  case PrimitivePoint:
    *out++ = primitive;
    break;

  case PrimitiveSegment:
    *out++ = primitive;
    break;

  case PrimitiveSurface:
    difference(primitive, *pb.as<PolygonWH_2>(), out);
    break;
  }

  return out;
}

template <typename OutputIteratorType>
auto
difference(const Point_3 &primitive, const PrimitiveHandle<3> &pb,
           OutputIteratorType out) -> OutputIteratorType
{
  switch (pb.handle.which()) {
  case PrimitivePoint:
    difference(primitive, *pb.as<Point_3>(), out);
    break;

  case PrimitiveSegment:
    difference(primitive, *pb.as<Segment_3>(), out);
    break;

  case PrimitiveSurface:
    difference(primitive, *pb.as<Triangle_3>(), out);
    break;

  case PrimitiveVolume:
    difference(primitive, *pb.as<MarkedPolyhedron>(), out);
    break;
  }

  return out;
}

template <typename OutputIteratorType>
auto
difference(const Segment_3 &primitive, const PrimitiveHandle<3> &pb,
           OutputIteratorType out) -> OutputIteratorType
{
  switch (pb.handle.which()) {
  case PrimitivePoint:
    *out++ = primitive;
    break;

  case PrimitiveSegment:
    difference(primitive, *pb.as<Segment_3>(), out);
    break;

  case PrimitiveSurface:
    difference(primitive, *pb.as<Triangle_3>(), out);
    break;

  case PrimitiveVolume:
    difference(primitive, *pb.as<MarkedPolyhedron>(), out);
    break;
  }

  return out;
}

template <typename OutputIteratorType>
auto
difference(const Triangle_3 &primitive, const PrimitiveHandle<3> &pb,
           OutputIteratorType out) -> OutputIteratorType
{
  switch (pb.handle.which()) {
  case PrimitivePoint:
    *out++ = primitive;
    break;

  case PrimitiveSegment:
    *out++ = primitive;
    break;

  case PrimitiveSurface:
    difference(primitive, *pb.as<Triangle_3>(), out);
    break;

  case PrimitiveVolume:
    difference(primitive, *pb.as<MarkedPolyhedron>(), out);
    break;
  }

  return out;
}

template <typename OutputIteratorType>
auto
difference(const MarkedPolyhedron &primitive, const PrimitiveHandle<3> &pb,
           OutputIteratorType out) -> OutputIteratorType
{
  switch (pb.handle.which()) {
  case PrimitivePoint:
    *out++ = primitive;
    break;

  case PrimitiveSegment:
    *out++ = primitive;
    break;

  case PrimitiveSurface:
    *out++ = primitive;
    break;

  case PrimitiveVolume:
    difference(primitive, *pb.as<MarkedPolyhedron>(), out);
    break;
  }

  return out;
}

template <typename Primitive, typename PrimitiveHandleConstIterator>
auto
difference(const Primitive &primitive, PrimitiveHandleConstIterator begin,
           PrimitiveHandleConstIterator end) -> std::vector<Primitive>
{
  std::vector<Primitive> primitives;
  primitives.push_back(primitive);

  for (PrimitiveHandleConstIterator b = begin; b != end; ++b) {
    std::vector<Primitive> new_primitives;

    for (auto a = primitives.begin();
         a != primitives.end(); ++a) {
      difference(*a, *(*b), std::back_inserter(new_primitives));
    }

    primitives.swap(new_primitives);
  }

  return primitives;
}

// just performs the type switch for the primitive to substract from
//
void
appendDifference(const PrimitiveHandle<2>                              &pa,
                 CollisionMapper<2>::PrimitiveHandleSet::const_iterator begin,
                 CollisionMapper<2>::PrimitiveHandleSet::const_iterator end,
                 GeometrySet<2>                                        &output)
{
  switch (pa.handle.which()) {
  case PrimitivePoint: {
    std::vector<Point_2> res = difference(*pa.as<Point_2>(), begin, end);
    output.addPoints(res.begin(), res.end());
    return;
  }

  case PrimitiveSegment: {
    std::vector<Segment_2> res = difference(*pa.as<Segment_2>(), begin, end);
    output.addSegments(res.begin(), res.end());
    return;
  }

  case PrimitiveSurface: {
    std::vector<PolygonWH_2> res =
        difference(*pa.as<PolygonWH_2>(), begin, end);
    output.addSurfaces(res.begin(), res.end());
    return;
  }
  }
}

void
appendDifference(const PrimitiveHandle<3>                              &pa,
                 CollisionMapper<3>::PrimitiveHandleSet::const_iterator begin,
                 CollisionMapper<3>::PrimitiveHandleSet::const_iterator end,
                 GeometrySet<3>                                        &output)
{
  switch (pa.handle.which()) {
  case PrimitivePoint: {
    std::vector<Point_3> res = difference(*pa.as<Point_3>(), begin, end);
    output.addPoints(res.begin(), res.end());
    return;
  }

  case PrimitiveSegment: {
    std::vector<Segment_3> res = difference(*pa.as<Segment_3>(), begin, end);
    output.addSegments(res.begin(), res.end());
    break;
  }

  case PrimitiveSurface: {
    std::vector<Triangle_3> res = difference(*pa.as<Triangle_3>(), begin, end);
    output.addSurfaces(res.begin(), res.end());
    break;
  }

  case PrimitiveVolume: {
    std::vector<MarkedPolyhedron> res =
        difference(*pa.as<MarkedPolyhedron>(), begin, end);
    output.addVolumes(res.begin(), res.end());
    break;
  }
  }
}

/**
 * difference post processing
 */
void
post_difference(const GeometrySet<2> &input, GeometrySet<2> &output)
{
  //
  // reverse orientation of polygons if needed
  for (const auto &it : input.surfaces()) {
    const PolygonWH_2 &p     = it.primitive();
    Polygon_2          outer = p.outer_boundary();

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

    output.surfaces().push_back(PolygonWH_2(outer, rings.begin(), rings.end()));
  }

  output.points()   = input.points();
  output.segments() = input.segments();
  output.volumes()  = input.volumes();
}

void
post_difference(const GeometrySet<3> &input, GeometrySet<3> &output)
{
  // nothing special to do
  output = input;
}

template <int Dim>
void
difference(const GeometrySet<Dim> &a, const GeometrySet<Dim> &b,
           GeometrySet<Dim> &output)
{
  using BoxCollection = typename SFCGAL::detail::BoxCollection<Dim>::Type;
  typename SFCGAL::detail::HandleCollection<Dim>::Type ahandles;
  typename SFCGAL::detail::HandleCollection<Dim>::Type bhandles;
  BoxCollection                                        aboxes;
  BoxCollection                                        bboxes;
  a.computeBoundingBoxes(ahandles, aboxes);
  b.computeBoundingBoxes(bhandles, bboxes);

  // here we use box_intersection_d to build the list of operations
  // that actually need to be performed
  GeometrySet<Dim>                   temp;
  GeometrySet<Dim>                   temp2;
  typename CollisionMapper<Dim>::Map map;
  CollisionMapper<Dim>               const cb(map);
  CGAL::box_intersection_d(aboxes.begin(), aboxes.end(), bboxes.begin(),
                           bboxes.end(), cb);

  // now we have in cb a map of operations to perform
  // we can put in the result right away all the primitives
  // that are not keys in this map
  {
    auto       ait = aboxes.begin();
    const auto end = aboxes.end();

    for (; ait != end; ++ait) {
      if (map.find(ait->handle()) == map.end()) {
        temp.addPrimitive(*ait->handle());
      }
    }
  }

  // then we delegate the operations according to type
  {
    auto       cbit = map.begin();
    const auto end  = map.end();

    for (; cbit != end; ++cbit) {
      appendDifference(*cbit->first, cbit->second.begin(), cbit->second.end(),
                       temp);
    }
  }

  post_difference(temp, temp2);
  output.merge(temp2);
}

template void
difference<2>(const GeometrySet<2> &a, const GeometrySet<2> &b,
              GeometrySet<2> &);
template void
difference<3>(const GeometrySet<3> &a, const GeometrySet<3> &b,
              GeometrySet<3> &);

auto
difference(const Geometry &ga, const Geometry &gb, NoValidityCheck /*unused*/)
    -> std::unique_ptr<Geometry>
{
  GeometrySet<2> gsa(ga);
  GeometrySet<2> gsb(gb);
  GeometrySet<2> output;
  algorithm::difference(gsa, gsb, output);

  GeometrySet<2> filtered;
  output.filterCovered(filtered);
  return filtered.recompose();
}

auto
difference(const Geometry &ga, const Geometry &gb) -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(ga);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(gb);

  return difference(ga, gb, NoValidityCheck());
}

auto
difference3D(const Geometry &ga, const Geometry &gb, NoValidityCheck /*unused*/)
    -> std::unique_ptr<Geometry>
{
  GeometrySet<3> gsa(ga);
  GeometrySet<3> gsb(gb);
  GeometrySet<3> output;
  algorithm::difference(gsa, gsb, output);

  GeometrySet<3> filtered;
  output.filterCovered(filtered);

  return filtered.recompose();
}

auto
difference3D(const Geometry &ga, const Geometry &gb)
    -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(ga);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(gb);

  return difference3D(ga, gb, NoValidityCheck());
}
} // namespace SFCGAL
