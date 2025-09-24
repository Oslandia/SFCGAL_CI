// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
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

auto
area(const Geometry &geom, NoValidityCheck /*unused*/) -> double
{
  switch (geom.geometryTypeId()) {
  case TYPE_POINT:
  case TYPE_LINESTRING:
    return 0;

  case TYPE_POLYGON:
    return area(geom.as<Polygon>());

  case TYPE_TRIANGLE:
    return area(geom.as<Triangle>());

  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_GEOMETRYCOLLECTION:
    return area(geom.as<GeometryCollection>());

  case TYPE_TRIANGULATEDSURFACE:
    return area(geom.as<TriangulatedSurface>());

  case TYPE_POLYHEDRALSURFACE:
    return area(geom.as<PolyhedralSurface>());

  case TYPE_SOLID:
  case TYPE_MULTISOLID:
    return 0;
  }

  BOOST_THROW_EXCEPTION(Exception(
      (boost::format(
           "Unexpected geometry type (%s) in SFCGAL::algorithm::area") %
       geom.geometryType())
          .str()));
}

auto
area(const Geometry &geom) -> double
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(geom);
  return area(geom, NoValidityCheck());
}

auto
signedArea(const Triangle &triangle) -> Kernel::FT
{
  Triangle_2 const triangle2D = triangle.toTriangle_2();
  return triangle2D.area();
}

auto
signedArea(const LineString &lineString) -> Kernel::FT
{
  return lineString.toPolygon_2(false).area();
}

auto
area(const Triangle &triangle) -> double
{
  return CGAL::to_double(CGAL::abs(signedArea(triangle)));
}

auto
area(const Polygon &polygon) -> double
{
  Kernel::RT result = 0.0;

  for (size_t i = 0; i < polygon.numRings(); i++) {
    Kernel::FT const ringArea = CGAL::abs(signedArea(polygon.ringN(i)));

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

auto
area(const GeometryCollection &collection) -> double
{
  double result = 0.0;

  for (size_t i = 0; i < collection.numGeometries(); i++) {
    result += area(collection.geometryN(i));
  }

  return result;
}

auto
area(const TriangulatedSurface &tin) -> double
{
  double result = 0.0;

  for (size_t i = 0; i < tin.numPatches(); i++) {
    result += area(tin.patchN(i));
  }

  return result;
}

auto
area(const PolyhedralSurface &surface) -> double
{
  double result = 0.0;

  for (size_t i = 0; i < surface.numPatches(); i++) {
    result += area(surface.patchN(i));
  }

  return result;
}

// ----------------------------------------------------------------------------------
//  -- area3D
// ----------------------------------------------------------------------------------

auto
area3D(const Geometry &geom, [[maybe_unused]] NoValidityCheck noCheck) -> double
{
  switch (geom.geometryTypeId()) {
  case TYPE_POINT:
  case TYPE_LINESTRING:
    return 0;

  case TYPE_POLYGON:
    return area3D(geom.as<Polygon>());

  case TYPE_TRIANGLE:
    return area3D(geom.as<Triangle>());

  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_GEOMETRYCOLLECTION:
    return area3D(geom.as<GeometryCollection>());

  case TYPE_TRIANGULATEDSURFACE:
    return area3D(geom.as<TriangulatedSurface>());

  case TYPE_POLYHEDRALSURFACE:
    return area3D(geom.as<PolyhedralSurface>());

  case TYPE_SOLID:
  case TYPE_MULTISOLID:
    return 0;
  }

  BOOST_THROW_EXCEPTION(Exception("missing case in SFCGAL::algorithm::area3D"));
}

auto
area3D(const Geometry &geom) -> double
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(geom);
  return area3D(geom, NoValidityCheck());
}

auto
area3D(const Polygon &polygon) -> double
{
  double result = 0.0;

  if (polygon.isEmpty()) {
    return result;
  }

  CGAL::Point_3<Kernel> a;
  CGAL::Point_3<Kernel> b;
  CGAL::Point_3<Kernel> c;
  algorithm::plane3D<Kernel>(polygon, a, b, c);

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
  for (size_t i = 0; i < polygon.numRings(); i++) {
    const LineString &ring = polygon.ringN(i);

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

auto
area3D(const Triangle &triangle) -> double
{
  CGAL::Triangle_3<Kernel> const triangle3D(triangle.vertex(0).toPoint_3(),
                                            triangle.vertex(1).toPoint_3(),
                                            triangle.vertex(2).toPoint_3());
  return sqrt(CGAL::to_double(triangle3D.squared_area()));
}

auto
area3D(const GeometryCollection &collection) -> double
{
  double result = 0.0;

  for (size_t i = 0; i < collection.numGeometries(); i++) {
    result += area3D(collection.geometryN(i));
  }

  return result;
}

auto
area3D(const PolyhedralSurface &surface) -> double
{
  double area = 0.0;

  for (size_t i = 0; i < surface.numPatches(); i++) {
    area += area3D(surface.patchN(i));
  }

  return area;
}

auto
area3D(const TriangulatedSurface &tin) -> double
{
  double result = 0.0;

  for (size_t i = 0; i < tin.numPatches(); i++) {
    result += area3D(tin.patchN(i));
  }

  return result;
}

} // namespace SFCGAL::algorithm
