// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/area.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"

#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/algorithm/plane.h"

#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangle_2.h>

#include <CGAL/Plane_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Triangle_3.h>

#include "SFCGAL/Exception.h"
#include <boost/format.hpp>

namespace SFCGAL::algorithm {

using Point_2    = CGAL::Point_2<SFCGAL::Kernel>;
using Triangle_2 = CGAL::Triangle_2<SFCGAL::Kernel>;
using Polygon_2  = CGAL::Polygon_2<SFCGAL::Kernel>;

using Point_3    = CGAL::Point_3<SFCGAL::Kernel>;
using Triangle_3 = CGAL::Triangle_3<SFCGAL::Kernel>;
using Plane_3    = CGAL::Plane_3<SFCGAL::Kernel>;

///
///
///
auto
area(const Geometry &g, NoValidityCheck /*unused*/) -> double
{
  switch (g.geometryTypeId()) {
  case TYPE_POINT:
  case TYPE_LINESTRING:
    return 0;

  case TYPE_POLYGON:
    return area(g.as<Polygon>());

  case TYPE_TRIANGLE:
    return area(g.as<Triangle>());

  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_GEOMETRYCOLLECTION:
    return area(g.as<GeometryCollection>());

  case TYPE_TRIANGULATEDSURFACE:
    return area(g.as<TriangulatedSurface>());

  case TYPE_POLYHEDRALSURFACE:
    return area(g.as<PolyhedralSurface>());

  case TYPE_SOLID:
  case TYPE_MULTISOLID:
    return 0;
  }

  BOOST_THROW_EXCEPTION(Exception(
      (boost::format(
           "Unexpected geometry type (%s) in SFCGAL::algorithm::area") %
       g.geometryType())
          .str()));
}

auto
area(const Geometry &g) -> double
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(g);
  return area(g, NoValidityCheck());
}

///
///
///
auto
signedArea(const Triangle &g) -> Kernel::FT
{
  Triangle_2 const triangle = g.toTriangle_2();
  return triangle.area();
}

///
///
///
auto
signedArea(const LineString &g) -> Kernel::FT
{
  return g.toPolygon_2(false).area();
}

///
///
///
auto
area(const Triangle &g) -> double
{
  return CGAL::to_double(CGAL::abs(signedArea(g)));
}

///
///
///
auto
area(const Polygon &g) -> double
{
  Kernel::RT result = 0.0;

  for (size_t i = 0; i < g.numRings(); i++) {
    Kernel::FT const ringArea = CGAL::abs(signedArea(g.ringN(i)));

    if (i == 0) {
      // exterior ring
      result += CGAL::abs(ringArea);
    } else {
      // interior ring
      result -= CGAL::abs(ringArea);
    }
  }

  return CGAL::to_double(result);
}

///
///
///
auto
area(const GeometryCollection &g) -> double
{
  double result = 0.0;

  for (size_t i = 0; i < g.numGeometries(); i++) {
    result += area(g.geometryN(i));
  }

  return result;
}

///
///
///
auto
area(const TriangulatedSurface &g) -> double
{
  double result = 0.0;

  for (size_t i = 0; i < g.numTriangles(); i++) {
    result += area(g.triangleN(i));
  }

  return result;
}

///
///
///
auto
area(const PolyhedralSurface &g) -> double
{
  double result = 0.0;

  for (size_t i = 0; i < g.numPolygons(); i++) {
    result += area(g.polygonN(i));
  }

  return result;
}

// ----------------------------------------------------------------------------------
//  -- area3D
// ----------------------------------------------------------------------------------

///
///
///
auto
area3D(const Geometry &g, NoValidityCheck /*unused*/) -> double
{
  switch (g.geometryTypeId()) {
  case TYPE_POINT:
  case TYPE_LINESTRING:
    return 0;

  case TYPE_POLYGON:
    return area3D(g.as<Polygon>());

  case TYPE_TRIANGLE:
    return area3D(g.as<Triangle>());

  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_GEOMETRYCOLLECTION:
    return area3D(g.as<GeometryCollection>());

  case TYPE_TRIANGULATEDSURFACE:
    return area3D(g.as<TriangulatedSurface>());

  case TYPE_POLYHEDRALSURFACE:
    return area3D(g.as<PolyhedralSurface>());

  case TYPE_SOLID:
  case TYPE_MULTISOLID:
    return 0;
  }

  BOOST_THROW_EXCEPTION(Exception("missing case in SFCGAL::algorithm::area3D"));
}

auto
area3D(const Geometry &g) -> double
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(g);
  return area3D(g, NoValidityCheck());
}

///
///
///
auto
area3D(const Polygon &g) -> double
{
  double result = 0.0;

  if (g.isEmpty()) {
    return result;
  }

  CGAL::Point_3<Kernel> a;
  CGAL::Point_3<Kernel> b;
  CGAL::Point_3<Kernel> c;
  algorithm::plane3D<Kernel>(g, a, b, c);

  /*
   * compute polygon basis (CGAL doesn't build an orthonormal basis so that
   * computing the 2D area in this basis would lead to scale effects) ux = bc uz
   * = bc^ba uy = uz^ux
   *
   * Note that the basis is rounded to double (CGAL::sqrt)
   */
  CGAL::Vector_3<Kernel> ux = c - b;
  CGAL::Vector_3<Kernel> uz = CGAL::cross_product(ux, a - b);
  ux = ux / CGAL::sqrt(CGAL::to_double(ux.squared_length()));
  uz = uz / CGAL::sqrt(CGAL::to_double(uz.squared_length()));
  CGAL::Vector_3<Kernel> const uy = CGAL::cross_product(uz, ux);

  /*
   * compute the area for each ring in the local basis
   */
  for (size_t i = 0; i < g.numRings(); i++) {
    const LineString &ring = g.ringN(i);

    CGAL::Polygon_2<Kernel> projectedPolygon;

    for (size_t j = 0; j < ring.numPoints() - 1; j++) {
      CGAL::Point_3<Kernel> const point = ring.pointN(j).toPoint_3();
      CGAL::Point_2<Kernel> const projectedPoint((point - b) * ux,
                                                 (point - b) * uy);
      projectedPolygon.push_back(projectedPoint);
    }

    if (i == 0) {
      // exterior ring
      result += CGAL::to_double(CGAL::abs(projectedPolygon.area()));
    } else {
      // interior ring
      result -= CGAL::to_double(CGAL::abs(projectedPolygon.area()));
    }
  }

  return result;
}

///
///
///
auto
area3D(const Triangle &g) -> double
{
  CGAL::Triangle_3<Kernel> const triangle(g.vertex(0).toPoint_3(),
                                          g.vertex(1).toPoint_3(),
                                          g.vertex(2).toPoint_3());
  return sqrt(CGAL::to_double(triangle.squared_area()));
}

///
///
///
auto
area3D(const GeometryCollection &g) -> double
{
  double result = 0.0;

  for (size_t i = 0; i < g.numGeometries(); i++) {
    result += area3D(g.geometryN(i));
  }

  return result;
}

///
///
///
auto
area3D(const PolyhedralSurface &g) -> double
{
  double area = 0.0;

  for (size_t i = 0; i < g.numPolygons(); i++) {
    area += area3D(g.polygonN(i));
  }

  return area;
}

///
///
///
auto
area3D(const TriangulatedSurface &g) -> double
{
  double result = 0.0;

  for (size_t i = 0; i < g.numGeometries(); i++) {
    result += area3D(g.geometryN(i));
  }

  return result;
}

} // namespace SFCGAL::algorithm
