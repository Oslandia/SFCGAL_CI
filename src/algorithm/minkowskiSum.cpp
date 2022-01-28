// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/algorithm/isValid.h>
#include <SFCGAL/algorithm/minkowskiSum.h>

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>

#include <SFCGAL/detail/polygonSetToMultiPolygon.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/minkowski_sum_2.h>

#include <CGAL/Aff_transformation_2.h>

using Polygon_2            = CGAL::Polygon_2<SFCGAL::Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<SFCGAL::Kernel>;
using Polygon_set_2        = CGAL::Polygon_set_2<SFCGAL::Kernel>;

namespace SFCGAL {
namespace algorithm {

//-- private interface

/**
 * dispatch gA+gB sum
 */
void
minkowskiSum(const Geometry &gA, const Polygon_2 &gB,
             CGAL::Polygon_set_2<Kernel> &polygonSet);
/*
 * append gA+gB into the polygonSet
 */
void
minkowskiSum(const Point &gA, const Polygon_2 &gB, Polygon_set_2 &polygonSet);
/*
 * append gA+gB into the polygonSet
 */
void
minkowskiSum(const LineString &gA, const Polygon_2 &gB,
             Polygon_set_2 &polygonSet);
/*
 * append gA+gB into the polygonSet
 */
void
minkowskiSum(const Polygon &gA, const Polygon_2 &gB, Polygon_set_2 &polygonSet);
/*
 * append gA+gB into the polygonSet
 */
void
minkowskiSum(const Solid &gA, const Polygon_2 &gB, Polygon_set_2 &polygonSet);
/*
 * append gA+gB into the polygonSet
 */
void
minkowskiSumCollection(const Geometry &gA, const Polygon_2 &gB,
                       Polygon_set_2 &polygonSet);

//-- private interface implementation

///
///
///
void
minkowskiSum(const Geometry &gA, const Polygon_2 &gB,
             CGAL::Polygon_set_2<Kernel> &polygonSet)
{
  if (gA.isEmpty()) {
    return;
  }

  switch (gA.geometryTypeId()) {
  case TYPE_POINT:
    return minkowskiSum(gA.as<Point>(), gB, polygonSet);

  case TYPE_LINESTRING:
    return minkowskiSum(gA.as<LineString>(), gB, polygonSet);

  case TYPE_POLYGON:
    return minkowskiSum(gA.as<Polygon>(), gB, polygonSet);

  case TYPE_TRIANGLE:
    return minkowskiSum(gA.as<Triangle>().toPolygon(), gB, polygonSet);

  case TYPE_SOLID:
    return minkowskiSum(gA.as<Solid>(), gB, polygonSet);

  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
    return minkowskiSumCollection(gA, gB, polygonSet);
  }

  BOOST_THROW_EXCEPTION(
      Exception((boost::format("minkowskiSum( %s, 'Polygon' ) is not defined") %
                 gA.geometryType())
                    .str()));
}

/*
 * append gA+gB into the polygonSet
 */
void
minkowskiSum(const Point &gA, const Polygon_2 &gB, Polygon_set_2 &polygonSet)
{
  BOOST_ASSERT(gB.size());

  CGAL::Aff_transformation_2<Kernel> translate(CGAL::TRANSLATION,
                                               gA.toVector_2());

  Polygon_2 sum;

  for (auto it = gB.vertices_begin(); it != gB.vertices_end(); ++it) {
    sum.push_back(translate.transform(*it));
  }

  if (sum.is_clockwise_oriented()) {
    sum.reverse_orientation();
  }

  if (polygonSet.is_empty()) {
    polygonSet.insert(sum);
  } else {
    polygonSet.join(sum);
  }
}

///
///
///
void
minkowskiSum(const LineString &gA, const Polygon_2 &gB,
             Polygon_set_2 &polygonSet)
{
  if (gA.isEmpty()) {
    return;
  }

  int npt = gA.numPoints();

  for (int i = 0; i < npt - 1; i++) {
    Polygon_2 P;
    P.push_back(gA.pointN(i).toPoint_2());
    P.push_back(gA.pointN(i + 1).toPoint_2());

    //
    // We want to compute the "minkowski sum" on each segment of the line string
    // This is not very well defined. But it appears CGAL supports it.
    // However we must use the explicit "full convolution" method for that
    // particular case in CGAL >= 4.7
#if CGAL_VERSION_NR < 1040701000 // version 4.7
    Polygon_with_holes_2 part = minkowski_sum_2(P, gB);
#else
    Polygon_with_holes_2 part = minkowski_sum_by_full_convolution_2(P, gB);
#endif

    // merge into a polygon set
    if (polygonSet.is_empty()) {
      polygonSet.insert(part);
    } else {
      polygonSet.join(part);
    }
  }
}

///
///
///
void
minkowskiSum(const Polygon &gA, const Polygon_2 &gB, Polygon_set_2 &polygonSet)
{
  if (gA.isEmpty()) {
    return;
  }

  /*
   * Invoke minkowski_sum_2 for exterior ring
   */
  {
    Polygon_with_holes_2 sum =
        minkowski_sum_2(gA.exteriorRing().toPolygon_2(), gB);

    if (polygonSet.is_empty()) {
      polygonSet.insert(sum);
    } else {
      polygonSet.join(sum);
    }
  }

  /*
   * Compute the Minkowski sum for each segment of the interior rings
   * and perform the union of the result. The result is a polygon, and its holes
   * correspond to the inset.
   *
   */
  if (gA.hasInteriorRings()) {
    Polygon_set_2 sumInteriorRings;

    for (size_t i = 0; i < gA.numInteriorRings(); i++) {
      minkowskiSum(gA.interiorRingN(i), gB, sumInteriorRings);
    }

    /*
     * compute the difference for each hole of the resulting polygons
     */
    std::list<Polygon_with_holes_2> interiorPolygons;
    sumInteriorRings.polygons_with_holes(std::back_inserter(interiorPolygons));

    for (auto &interiorPolygon : interiorPolygons) {

      for (auto it_hole = interiorPolygon.holes_begin();
           it_hole != interiorPolygon.holes_end(); ++it_hole) {

        it_hole->reverse_orientation();
        polygonSet.difference(*it_hole);
      } // foreach hole
    }   // foreach polygon
  }
}

///
///
///
void
minkowskiSum(const Solid &gA, const Polygon_2 &gB, Polygon_set_2 &polygonSet)
{
  // use only the projection of exterior shell
  minkowskiSumCollection(gA.exteriorShell(), gB, polygonSet);
}

///
///
///
void
minkowskiSumCollection(const Geometry &gA, const Polygon_2 &gB,
                       Polygon_set_2 &polygonSet)
{
  for (size_t i = 0; i < gA.numGeometries(); i++) {
    minkowskiSum(gA.geometryN(i), gB, polygonSet);
  }
}

auto
minkowskiSum(const Geometry &gA, const Polygon &gB, NoValidityCheck)
    -> std::unique_ptr<Geometry>
{
  if (gB.isEmpty()) {
    return std::unique_ptr<Geometry>(gA.clone());
  }

  Polygon_set_2 polygonSet;
  minkowskiSum(gA, gB.toPolygon_2(), polygonSet);
  return std::unique_ptr<Geometry>(
      detail::polygonSetToMultiPolygon(polygonSet).release());
}

//-- public interface implementation

///
///
///
auto
minkowskiSum(const Geometry &gA, const Polygon &gB) -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(gA);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(gB);

  std::unique_ptr<Geometry> result(minkowskiSum(gA, gB, NoValidityCheck()));
  propagateValidityFlag(*result, true);
  return result;
}

} // namespace algorithm
} // namespace SFCGAL
