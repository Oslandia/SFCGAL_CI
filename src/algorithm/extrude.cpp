// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/algorithm/extrude.h>

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/MultiSolid.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>

#include <SFCGAL/Exception.h>

#include <SFCGAL/algorithm/force3D.h>
#include <SFCGAL/algorithm/isValid.h>
#include <SFCGAL/algorithm/normal.h>
#include <SFCGAL/algorithm/translate.h>

#include <SFCGAL/detail/tools/Log.h>

#include <utility>

namespace SFCGAL::algorithm {

//-- private interface

auto
extrude(const Point &g, const Kernel::Vector_3 &v) -> LineString *;
auto
extrude(const LineString &g, const Kernel::Vector_3 &v) -> PolyhedralSurface *;
auto
extrude(const Polygon &g, const Kernel::Vector_3 &v, bool addTop = true)
    -> Solid *;
auto
extrude(const Triangle &g, const Kernel::Vector_3 &v) -> Solid *;

auto
extrude(const MultiPoint &g, const Kernel::Vector_3 &v) -> MultiLineString *;
auto
extrude(const MultiLineString &g, const Kernel::Vector_3 &v)
    -> PolyhedralSurface *;
auto
extrude(const MultiPolygon &g, const Kernel::Vector_3 &v) -> MultiSolid *;

/**
 * @warning suppose that the TriangulatedSurface is connected
 * @todo take orientation in account
 */
auto
extrude(const TriangulatedSurface &g, const Kernel::Vector_3 &v) -> Solid *;
/**
 * @warning doesn't take orientation in account
 * @todo take orientation in account
 */
auto
extrude(const PolyhedralSurface &g, const Kernel::Vector_3 &v) -> Solid *;

/**
 * extrude each geometry in a GeometryCollection
 */
auto
extrude(const GeometryCollection &g, const Kernel::Vector_3 &v)
    -> GeometryCollection *;

///
///
///
auto
extrude(const Point &g, const Kernel::Vector_3 &v) -> LineString *
{
  if (g.isEmpty()) {
    return new LineString();
  }

  Kernel::Point_3 const a = g.toPoint_3();
  Kernel::Point_3 const b = a + v;

  return new LineString(Point(a), Point(b));
}

///
///
///
auto
extrude(const LineString &g, const Kernel::Vector_3 &v) -> PolyhedralSurface *
{

  std::unique_ptr<PolyhedralSurface> polyhedralSurface(new PolyhedralSurface());

  if (g.isEmpty()) {
    return polyhedralSurface.release();
  }

  for (size_t i = 0; i < g.numPoints() - 1; i++) {
    std::unique_ptr<LineString> ring(new LineString);

    Kernel::Point_3 const a = g.pointN(i).toPoint_3();
    Kernel::Point_3 const b = g.pointN(i + 1).toPoint_3();
    ring->addPoint(new Point(a));
    ring->addPoint(new Point(b));
    ring->addPoint(new Point(b + v));
    ring->addPoint(new Point(a + v));
    ring->addPoint(new Point(a));

    polyhedralSurface->addPolygon(new Polygon(ring.release()));
  }

  return polyhedralSurface.release();
}

///
///
///
auto
extrude(const Polygon &g, const Kernel::Vector_3 &v, bool addTop) -> Solid *
{
  if (g.isEmpty()) {
    return new Solid();
  }

  bool const reverseOrientation = (v * normal3D<Kernel>(g)) > 0;

  // resulting shell
  PolyhedralSurface polyhedralSurface;

  // "bottom"
  Polygon bottom(g);
  force3D(bottom);

  if (reverseOrientation) {
    bottom.reverse();
  }

  polyhedralSurface.addPolygon(bottom);

  // "top"
  if (addTop) {
    Polygon top(bottom);
    top.reverse();
    translate(top, v);
    polyhedralSurface.addPolygon(top);
  }
  // exterior ring and interior rings extruded
  for (size_t i = 0; i < bottom.numRings(); i++) {
    std::unique_ptr<PolyhedralSurface> boundaryExtruded(
        extrude(bottom.ringN(i), v));

    for (size_t j = 0; j < boundaryExtruded->numPolygons(); j++) {
      boundaryExtruded->polygonN(j).reverse();
      polyhedralSurface.addPolygon(boundaryExtruded->polygonN(j));
    }
  }

  return new Solid(polyhedralSurface);
}

///
///
///
auto
extrude(const Triangle &g, const Kernel::Vector_3 &v) -> Solid *
{
  return extrude(g.toPolygon(), v);
}

///
///
///
auto
extrude(const MultiPoint &g, const Kernel::Vector_3 &v) -> MultiLineString *
{
  std::unique_ptr<MultiLineString> result(new MultiLineString());

  if (g.isEmpty()) {
    return result.release();
  }

  for (size_t i = 0; i < g.numGeometries(); i++) {
    result->addGeometry(extrude(g.pointN(i), v));
  }

  return result.release();
}

///
///
///
auto
extrude(const MultiLineString &g, const Kernel::Vector_3 &v)
    -> PolyhedralSurface *
{
  std::unique_ptr<PolyhedralSurface> result(new PolyhedralSurface());

  if (g.isEmpty()) {
    return result.release();
  }

  for (size_t i = 0; i < g.numGeometries(); i++) {
    std::unique_ptr<PolyhedralSurface> extruded(extrude(g.lineStringN(i), v));

    for (size_t j = 0; j < extruded->numPolygons(); j++) {
      result->addPolygon(extruded->polygonN(j));
    }
  }

  return result.release();
}

///
///
///
auto
extrude(const MultiPolygon &g, const Kernel::Vector_3 &v) -> MultiSolid *
{
  std::unique_ptr<MultiSolid> result(new MultiSolid());

  if (g.isEmpty()) {
    return result.release();
  }

  for (size_t i = 0; i < g.numGeometries(); i++) {
    result->addGeometry(extrude(g.polygonN(i), v));
  }

  return result.release();
}

///
///
///
auto
extrude(const TriangulatedSurface &g, const Kernel::Vector_3 &v) -> Solid *
{
  std::unique_ptr<Solid> result(new Solid());

  if (g.isEmpty()) {
    return result.release();
  }

  // bottom and top
  for (size_t i = 0; i < g.numGeometries(); i++) {
    Triangle bottomPart(g.geometryN(i));
    force3D(bottomPart);
    bottomPart.reverse();
    result->exteriorShell().addPolygon(bottomPart);

    Triangle topPart(g.geometryN(i));
    force3D(topPart);
    translate(topPart, v);
    result->exteriorShell().addPolygon(topPart);
  }

  // boundary
  std::unique_ptr<Geometry> boundary(g.boundary());
  BOOST_ASSERT(boundary.get() != NULL);

  // closed surface extruded
  if (!boundary->isEmpty()) {
    std::unique_ptr<Geometry> extrudedBoundary(extrude(*boundary, v));
    BOOST_ASSERT(extrudedBoundary->is<PolyhedralSurface>());
    result->exteriorShell().addPolygons(
        extrudedBoundary->as<PolyhedralSurface>());
  }

  return result.release();
}

///
///
///
auto
extrude(const PolyhedralSurface &g, const Kernel::Vector_3 &v) -> Solid *
{
  if (g.isEmpty()) {
    return new Solid();
  }

  TriangulatedSurface triangulatedSurface;
  triangulate::triangulatePolygon3D(g, triangulatedSurface);
  return extrude(triangulatedSurface, v);
}

///
///
///
auto
extrude(const GeometryCollection &g, const Kernel::Vector_3 &v)
    -> GeometryCollection *
{
  std::unique_ptr<GeometryCollection> result(new GeometryCollection());

  if (g.isEmpty()) {
    return result.release();
  }

  for (size_t i = 0; i < g.numGeometries(); i++) {
    result->addGeometry(extrude(g.geometryN(i), v).release());
  }

  return result.release();
}

//-- public interface

///
///
///
auto
extrude(const Geometry &g, const Kernel::Vector_3 &v)
    -> std::unique_ptr<Geometry>
{
  switch (g.geometryTypeId()) {
  case TYPE_POINT:
    return std::unique_ptr<Geometry>(extrude(g.as<Point>(), v));

  case TYPE_LINESTRING:
    return std::unique_ptr<Geometry>(extrude(g.as<LineString>(), v));

  case TYPE_POLYGON:
    return std::unique_ptr<Geometry>(extrude(g.as<Polygon>(), v));

  case TYPE_TRIANGLE:
    return std::unique_ptr<Geometry>(extrude(g.as<Triangle>(), v));

  case TYPE_GEOMETRYCOLLECTION:
    return std::unique_ptr<Geometry>(extrude(g.as<GeometryCollection>(), v));

  case TYPE_MULTIPOINT:
    return std::unique_ptr<Geometry>(extrude(g.as<MultiPoint>(), v));

  case TYPE_MULTILINESTRING:
    return std::unique_ptr<Geometry>(extrude(g.as<MultiLineString>(), v));

  case TYPE_MULTIPOLYGON:
    return std::unique_ptr<Geometry>(extrude(g.as<MultiPolygon>(), v));

  case TYPE_TRIANGULATEDSURFACE:
    return std::unique_ptr<Geometry>(extrude(g.as<TriangulatedSurface>(), v));

  case TYPE_POLYHEDRALSURFACE:
    return std::unique_ptr<Geometry>(extrude(g.as<PolyhedralSurface>(), v));

  case TYPE_SOLID:
  case TYPE_MULTISOLID:
    // extrusion not available
    break;
  }

  BOOST_THROW_EXCEPTION(InappropriateGeometryException(
      (boost::format("Extrusion of %s is not supported") % g.geometryType())
          .str()));
}

///
///
///
auto
extrude(const Geometry &g, const Kernel::FT &dx, const Kernel::FT &dy,
        const Kernel::FT &dz, NoValidityCheck /*unused*/) -> std::unique_ptr<Geometry>
{
  return extrude(g, Kernel::Vector_3(dx, dy, dz));
}

auto
extrude(const Geometry &g, const Kernel::FT& dx, const Kernel::FT& dy, const Kernel::FT& dz)
    -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY(g);
  std::unique_ptr<Geometry> result(extrude(g, dx, dy,
                                           dz, NoValidityCheck()));
  propagateValidityFlag(*result, true);
  return result;
}

SFCGAL_API auto
extrude(const Geometry &g, const double &dx, const double &dy, const double &dz)
    -> std::unique_ptr<Geometry>
{
  if (!std::isfinite(dx) || !std::isfinite(dy) || !std::isfinite(dz)) {
    BOOST_THROW_EXCEPTION(NonFiniteValueException(
        "trying to extrude with non finite value in direction"));
  }

  return extrude(g, Kernel::FT(dx), Kernel::FT(dy), Kernel::FT(dz));
}

SFCGAL_API auto
extrude(const Polygon &g, const double &height) -> std::unique_ptr<Geometry>
{

  if (!std::isfinite(height)) {
    BOOST_THROW_EXCEPTION(NonFiniteValueException(
        "trying to extrude with non finite value in direction"));
  }

  return std::unique_ptr<Geometry>(
      extrude(g, Kernel::Vector_3(0.0, 0.0, height), false));
}
} // namespace SFCGAL
