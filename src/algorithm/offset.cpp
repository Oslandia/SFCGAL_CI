// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/offset.h"

#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"

#include "SFCGAL/Exception.h"

#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/detail/polygonSetToMultiPolygon.h"

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Polygon_with_holes_2.h>

#include <CGAL/approximated_offset_2.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/offset_polygon_2.h>

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

using Gps_traits_2   = CGAL::Gps_circle_segment_traits_2<SFCGAL::Kernel>;
using Offset_curve_2 = Gps_traits_2::Curve_2;
using Offset_x_monotone_curve_2   = Gps_traits_2::X_monotone_curve_2;
using Offset_polygon_2            = Gps_traits_2::Polygon_2;
using Offset_polygon_with_holes_2 = Gps_traits_2::Polygon_with_holes_2;
using Offset_polygon_set_2        = CGAL::General_polygon_set_2<Gps_traits_2>;

#define SFCGAL_OFFSET_ACCURACY 0.0001

#define SFCGAL_OFFSET_ASSERT_FINITE_RADIUS(r)                                  \
  if (!std::isfinite(r))                                                       \
    BOOST_THROW_EXCEPTION(NonFiniteValueException("radius is non finite"));
namespace SFCGAL::algorithm {

/**
 * @brief dispatch a geometry
 */
void
offset(const Geometry &g, const double &radius,
       Offset_polygon_set_2 &polygonSet);

/**
 * @brief offset for a Point
 */
void
offset(const Point &g, const double &radius, Offset_polygon_set_2 &polygonSet);
/**
 * @brief offset for a LineString
 */
void
offset(const LineString &g, const double &radius,
       Offset_polygon_set_2 &polygonSet);
/**
 * @brief offset for a Polygon
 */
void
offset(const Polygon &g, const double &radius,
       Offset_polygon_set_2 &polygonSet);
/**
 * @brief offset for MultiPoint, MultiLineString, MultiPolygon,
 * TriangulatedSurface, PolyhedralSurface
 */
void
offsetCollection(const Geometry &g, const double &radius,
                 Offset_polygon_set_2 &polygonSet);

//-- helpers

/**
 * @brief approximate an Offset_polygon_2 (filter null segments)
 */
auto
approximate(const Offset_polygon_2 &polygon, const int &n = 0) -> Polygon_2
{
  std::list<std::pair<double, double>> pair_list;

  /*
   * iterate X_monotone_curve_2 components
   */
  for (auto it = polygon.curves_begin(); it != polygon.curves_end(); ++it) {
    it->approximate(std::back_inserter(pair_list), n);
  }

  // remove duplicated last point
  if (!pair_list.empty()) {
    pair_list.pop_back();
  }

  /*
   * convertr to polygon
   */
  Polygon_2 result;

  bool            isFirst = true;
  Kernel::Point_2 last;

  for (auto &it : pair_list) {
    Kernel::Point_2 const point(it.first, it.second);

    if (isFirst) {
      isFirst = false;
    } else if (point == last) {
      continue;
    }

    result.push_back(point);
    last = point;
  }

  return result;
}

/**
 * @brief approximate an Offset
 */
auto
approximate(const Offset_polygon_with_holes_2 &polygon, const int &n = 0)
    -> Polygon_with_holes_2
{
  Polygon_with_holes_2 result(approximate(polygon.outer_boundary(), n));

  for (auto it = polygon.holes_begin(); it != polygon.holes_end(); ++it) {
    result.add_hole(approximate(*it, n));
  }

  return result;
}

/**
 * @brief convert Offset_polygon_set_2 to MultiPolygon
 */
auto
polygonSetToMultiPolygon(const Offset_polygon_set_2 &polygonSet, const int &n)
    -> std::unique_ptr<MultiPolygon>
{
  std::list<Offset_polygon_with_holes_2> res;
  polygonSet.polygons_with_holes(std::back_inserter(res));

  std::unique_ptr<MultiPolygon> result(new MultiPolygon);

  for (auto &re : res) {
    result->addGeometry(new Polygon(approximate(re, n)));
  }

  return result;
}

/**
 * @brief helper to create a polygon from a circle
 */
auto
circleToPolygon(const Kernel::Circle_2 &circle) -> Offset_polygon_2
{
  /*
   * convert the circle into Offset_x_monotone_curve_2 (exactly 2)
   */
  Gps_traits_2 const   traits;
  Offset_curve_2 const curve(circle);

#if CGAL_VERSION_MAJOR < 6
  std::list<CGAL::Object> parts;
  traits.make_x_monotone_2_object()(curve, std::back_inserter(parts));
  BOOST_ASSERT(parts.size() == 2U);

  // Construct the polygon.
  Offset_polygon_2 result;

  for (auto &part : parts) {
    Offset_x_monotone_curve_2 arc;
    CGAL::assign(arc, part);
    result.push_back(arc);
  }
#else
  Offset_polygon_2 result;

  traits.make_x_monotone_2_object()(
      curve, CGAL::dispatch_or_drop_output<Offset_x_monotone_curve_2>(
                 std::back_inserter(result)));
#endif

  return result;
}

/**
 * @brief build Point offset
 */
void
offset(const Point &gA, const double &radius, Offset_polygon_set_2 &polygonSet)
{
  SFCGAL_OFFSET_ASSERT_FINITE_RADIUS(radius);
  Kernel::Circle_2 const circle(gA.toPoint_2(), radius * radius);

  if (polygonSet.is_empty()) {
    polygonSet.insert(circleToPolygon(circle));
  } else {
    polygonSet.join(circleToPolygon(circle));
  }
}

/**
 * @brief build LineString offset
 */
void
offset(const LineString &lineString, const double &radius,
       Offset_polygon_set_2 &polygonSet)
{
  SFCGAL_OFFSET_ASSERT_FINITE_RADIUS(radius);

  for (size_t i = 0; i < lineString.numSegments(); i++) {
    Polygon_2 P;
    P.push_back(lineString.pointN(i).toPoint_2());
    P.push_back(lineString.pointN(i + 1).toPoint_2());
    Offset_polygon_with_holes_2 const offset =
        CGAL::approximated_offset_2(P, radius, SFCGAL_OFFSET_ACCURACY);

    if (polygonSet.is_empty()) {
      polygonSet.insert(offset);
    } else {
      polygonSet.join(offset);
    }
  }
}

void
offset(const Polygon &g, const double &radius, Offset_polygon_set_2 &polygonSet)
{
  SFCGAL_OFFSET_ASSERT_FINITE_RADIUS(radius);

  if (g.isEmpty()) {
    return;
  }

  /*
   * Invoke minkowski_sum_2 for exterior ring
   */
  {
    Offset_polygon_with_holes_2 const offset = CGAL::approximated_offset_2(
        g.exteriorRing().toPolygon_2(), radius, SFCGAL_OFFSET_ACCURACY);

    if (polygonSet.is_empty()) {
      polygonSet.insert(offset);
    } else {
      polygonSet.join(offset);
    }
  }

  /*
   * Compute the Minkowski sum for each segment of the interior rings
   * and perform the union of the result. The result is a polygon, and its holes
   * correspond to the inset.
   *
   */
  if (g.hasInteriorRings()) {
    Offset_polygon_set_2 sumInteriorRings;

    for (size_t i = 0; i < g.numInteriorRings(); i++) {
      offset(g.interiorRingN(i), radius, sumInteriorRings);
    }

    /*
     * compute the difference for each hole of the resulting polygons
     */
    std::list<Offset_polygon_with_holes_2> interiorPolygons;
    sumInteriorRings.polygons_with_holes(std::back_inserter(interiorPolygons));

    for (auto &interiorPolygon : interiorPolygons) {

      for (auto it_hole = interiorPolygon.holes_begin();
           it_hole != interiorPolygon.holes_end(); ++it_hole) {

        it_hole->reverse_orientation();
        polygonSet.difference(*it_hole);
      } // foreach hole
    } // foreach polygon
  }
}

void
offsetCollection(const Geometry &g, const double &radius,
                 Offset_polygon_set_2 &polygonSet)
{
  SFCGAL_OFFSET_ASSERT_FINITE_RADIUS(radius);

  for (size_t i = 0; i < g.numGeometries(); i++) {
    offset(g.geometryN(i), radius, polygonSet);
  }
}

void
offset(const Geometry &g, const double &radius,
       Offset_polygon_set_2 &polygonSet)
{
  SFCGAL_OFFSET_ASSERT_FINITE_RADIUS(radius);

  if (g.isEmpty()) {
    return;
  }

  switch (g.geometryTypeId()) {
  case TYPE_POINT:
    return offset(g.as<Point>(), radius, polygonSet);

  case TYPE_LINESTRING:
    return offset(g.as<LineString>(), radius, polygonSet);

  case TYPE_POLYGON:
    return offset(g.as<Polygon>(), radius, polygonSet);

  case TYPE_TRIANGLE:
    return offset(g.as<Triangle>().toPolygon(), radius, polygonSet);

  case TYPE_SOLID:
    return offset(g.as<Solid>().exteriorShell(), radius, polygonSet);

  case TYPE_MULTISOLID:
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
    return offsetCollection(g, radius, polygonSet);
  }
}

/// @} end of private section

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

auto
offset(const Geometry &g, const double &r, NoValidityCheck /*unused*/)
    -> std::unique_ptr<MultiPolygon>
{
  SFCGAL_OFFSET_ASSERT_FINITE_RADIUS(r);
  Offset_polygon_set_2 polygonSet;
  offset(g, r, polygonSet);
  return polygonSetToMultiPolygon(polygonSet, 8);
}

auto
offset(const Geometry &g, const double &r) -> std::unique_ptr<MultiPolygon>
{
  SFCGAL_OFFSET_ASSERT_FINITE_RADIUS(r);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY(g);
  return offset(g, r, NoValidityCheck());
}

} // namespace SFCGAL::algorithm
