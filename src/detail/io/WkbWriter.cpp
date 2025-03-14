// Copyright (c) 2023-2023, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <array>
#include <boost/endian/conversion.hpp> // don't use bit, since it requires c++20
#include <cstddef>

#include "SFCGAL/detail/io/WkbWriter.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"

#include <boost/exception/all.hpp>
#include <exception>

namespace SFCGAL::detail::io {

void
WkbWriter::writeRec(const Geometry &g, boost::endian::order wkbOrder)
{
  switch (g.geometryTypeId()) {
  case TYPE_POINT:
    writeInner(g.as<Point>(), wkbOrder);
    return;

  case TYPE_LINESTRING:
    writeInner(g.as<LineString>(), wkbOrder);
    return;

  case TYPE_POLYGON:
    writeInner(g.as<Polygon>(), wkbOrder);
    return;

  case TYPE_GEOMETRYCOLLECTION:
    writeInner(g.as<GeometryCollection>(), wkbOrder);
    return;

  case TYPE_MULTIPOINT:
    writeInner<MultiPoint, Point>(g.as<MultiPoint>(), wkbOrder);
    return;

  case TYPE_MULTILINESTRING:
    writeInner<MultiLineString, LineString>(g.as<MultiLineString>(), wkbOrder);
    return;

  case TYPE_MULTIPOLYGON:
    writeInner<MultiPolygon, Polygon>(g.as<MultiPolygon>(), wkbOrder);
    return;

  case TYPE_TRIANGLE:
    writeInner(g.as<Triangle>(), wkbOrder);
    return;

  case TYPE_TRIANGULATEDSURFACE:
    writeInner<TriangulatedSurface, Triangle>(g.as<TriangulatedSurface>(),
                                              wkbOrder);
    return;

  case TYPE_POLYHEDRALSURFACE:
    writeInner<PolyhedralSurface, Polygon>(g.as<PolyhedralSurface>(), wkbOrder);
    return;

  default:
    std::ostringstream oss;
    oss << "WkbWriter: type '" << g.geometryType() << "' is not supported";
    BOOST_THROW_EXCEPTION(std::runtime_error(oss.str()));
  }
}

void
WkbWriter::write(const Geometry &g, const srid_t &srid,
                 boost::endian::order wkbOrder)
{

  _useSrid = true;
  _isEWKB  = true;
  _srid    = srid;

  write(g, wkbOrder);
}

void
WkbWriter::write(const Geometry &g, boost::endian::order wkbOrder)
{
  writeRec(g, wkbOrder);
}

void
WkbWriter::writeCoordinate(const Point &g, boost::endian::order wkbOrder)
{
  // x
  toByte(CGAL::to_double(g.x()), wkbOrder);

  // y
  toByte(CGAL::to_double(g.y()), wkbOrder);

  // z
  if (g.is3D()) {
    toByte(CGAL::to_double(g.z()), wkbOrder);
  }

  // m coordinate
  if (g.isMeasured()) {
    toByte(CGAL::to_double(g.m()), wkbOrder);
  }
}

void
WkbWriter::writeGeometryType(const Geometry &g, boost::endian::order wkbOrder)
{

  if (_isEWKB) {

    uint32_t ewkbtype = g.geometryTypeId();

    if (g.is3D()) {
      ewkbtype |= wkbZ;
    }
    if (g.isMeasured()) {
      ewkbtype |= wkbM;
    }
    if (_useSrid) {
      ewkbtype |= wkbSRID;
    }
    toByte<uint32_t>(ewkbtype, wkbOrder);

    if (_useSrid) {
      toByte<uint32_t>(_srid, wkbOrder);

      // once srid defined and written, we write geometry type without it
      // so "disable" it.
      _useSrid = false;
    }
  } else {
    toByte(static_cast<uint32_t>(g.geometryTypeId() +
                                 static_cast<int>(g.is3D()) * COORDINATE_XYZ +
                                 static_cast<int>(g.isMeasured()) *
                                     COORDINATE_XYM),
           wkbOrder);
  }
}

void
WkbWriter::writeInner(const Point &g, boost::endian::order wkbOrder)
{
  // Endianness
  toStream(std::array<std::byte, 1>{static_cast<std::byte>(wkbOrder)});

  // WkbType
  writeGeometryType(g, wkbOrder);

  if (g.isEmpty()) {
    toByte(std::numeric_limits<double>::quiet_NaN(), wkbOrder);
    toByte(std::numeric_limits<double>::quiet_NaN(), wkbOrder);
    if (g.is3D()) {
      toByte(std::numeric_limits<double>::quiet_NaN(), wkbOrder);
    }
    if (g.isMeasured()) {
      toByte(std::numeric_limits<double>::quiet_NaN(), wkbOrder);
    }
  } else {
    writeCoordinate(g, wkbOrder);
  }
}

void
WkbWriter::writeInner(const LineString &g, boost::endian::order wkbOrder)
{
  // Endianness
  toStream(std::array<std::byte, 1>{static_cast<std::byte>(wkbOrder)});

  // WkbType
  writeGeometryType(g, wkbOrder);

  writeInnerRing(g, wkbOrder);
}
void
WkbWriter::writeInnerRing(const LineString &g, boost::endian::order wkbOrder)
{
  toByte(static_cast<uint32_t>(g.numPoints()), wkbOrder);

  for (size_t i = 0; i < g.numPoints(); i++) {
    writeCoordinate(g.pointN(i), wkbOrder);
  }
}

void
WkbWriter::writeInner(const Polygon &g, boost::endian::order wkbOrder)
{
  // Endianness
  toStream(std::array<std::byte, 1>{static_cast<std::byte>(wkbOrder)});

  // WkbType
  writeGeometryType(g, wkbOrder);

  toByte(static_cast<uint32_t>(g.numRings()), wkbOrder);

  writeInnerRing(g.exteriorRing(), wkbOrder);
  for (size_t i = 0; i < g.numInteriorRings(); i++) {
    writeInnerRing(g.interiorRingN(i), wkbOrder);
  }
}

void
WkbWriter::writeInner(const Triangle &g, boost::endian::order wkbOrder)
{
  // Endianness
  toStream(std::array<std::byte, 1>{static_cast<std::byte>(wkbOrder)});

  // WkbType
  writeGeometryType(g, wkbOrder);

  if (!g.isEmpty()) {
    // One Ring
    toByte(static_cast<uint32_t>(1), wkbOrder);
    // 4 points
    toByte(static_cast<uint32_t>(4), wkbOrder);
    for (int i = 0; i < 4; i++) {
      writeCoordinate(g.vertex(i), wkbOrder);
    }
  }
}

void
WkbWriter::writeInner(const GeometryCollection &g,
                      boost::endian::order      wkbOrder)
{
  // Endianness
  toStream(std::array<std::byte, 1>{static_cast<std::byte>(wkbOrder)});

  // WkbType
  writeGeometryType(g, wkbOrder);

  // Number of Geometries
  toByte(static_cast<uint32_t>(g.numGeometries()), wkbOrder);

  for (size_t i = 0; i < g.numGeometries(); i++) {
    writeRec(g.geometryN(i), wkbOrder);
  }
}

template <typename M, typename G>
void
WkbWriter::writeInner(const M &g, boost::endian::order wkbOrder)
{
  // Endianness
  toStream(std::array<std::byte, 1>{static_cast<std::byte>(wkbOrder)});

  // WkbType
  writeGeometryType(g, wkbOrder);

  // Number of Geometries
  toByte(static_cast<uint32_t>(g.numGeometries()), wkbOrder);

  for (size_t i = 0; i < g.numGeometries(); i++) {
    writeInner(g.geometryN(i).template as<G>(), wkbOrder);
  }
}

} // namespace SFCGAL::detail::io
