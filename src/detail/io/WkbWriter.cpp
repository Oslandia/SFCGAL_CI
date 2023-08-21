// Copyright (c) 2023-2023, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <array>
#include <boost/endian/conversion.hpp> // don't use bit, since it requires c++20
#include <cstddef>
#include <iostream>

#include <SFCGAL/detail/io/WkbWriter.h>

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>

#include <boost/exception/all.hpp>
#include <exception>

namespace SFCGAL::detail::io {

static auto
toByte(const double x, boost::endian::order byteOrder)
    -> std::array<std::byte, 8>
{
  double y = x;
  if (boost::endian::order::native == boost::endian::order::big &&
      byteOrder == boost::endian::order::little) {
    boost::endian::endian_reverse_inplace(y);
  } else if (boost::endian::order::native == boost::endian::order::little &&
             byteOrder == boost::endian::order::big) {
    boost::endian::endian_reverse_inplace(y);
  }
  return *reinterpret_cast<std::array<std::byte, 8> *>(&y);
}

static auto
toByte(const uint32_t x, boost::endian::order byteOrder)
    -> std::array<std::byte, 4>
{
  uint32_t y = x;
  if (boost::endian::order::native == boost::endian::order::big &&
      byteOrder == boost::endian::order::little) {
    boost::endian::endian_reverse_inplace(y);
  } else if (boost::endian::order::native == boost::endian::order::little &&
             byteOrder == boost::endian::order::big) {
    boost::endian::endian_reverse_inplace(y);
  }
  return *reinterpret_cast<std::array<std::byte, 4> *>(&y);
}

auto
WkbWriter::toString(const bool asHex) -> std::string
{
  std::string       prefix{asHex ? "\\x" : ""};
  std::stringstream ss;
  for (const std::byte &byteVal : _wkb) {
    ss << prefix << std::hex << std::setw(2) << std::setfill('0')
       << static_cast<int>(byteVal);
  }
  return ss.str();
}

///
///
///
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
    writeInner(g.as<MultiPoint>(), wkbOrder);
    return;

  case TYPE_MULTILINESTRING:
    writeInner(g.as<MultiLineString>(), wkbOrder);
    return;

  case TYPE_MULTIPOLYGON:
    writeInner(g.as<MultiPolygon>(), wkbOrder);
    return;

  case TYPE_TRIANGLE:
    writeInner(g.as<Triangle>(), wkbOrder);
    return;

  case TYPE_TRIANGULATEDSURFACE:
    writeInner(g.as<TriangulatedSurface>(), wkbOrder);
    return;

  case TYPE_POLYHEDRALSURFACE:
    writeInner(g.as<PolyhedralSurface>(), wkbOrder);
    return;

  default:
    std::ostringstream oss;
    oss << "WkbWriter : '" << g.geometryType() << "' is not supported";
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

///
///
///
void
WkbWriter::write(const Geometry &g, boost::endian::order wkbOrder)
{
  switch (g.geometryTypeId()) {
  case TYPE_POINT:
    write(g.as<Point>(), wkbOrder);
    return;

  case TYPE_LINESTRING:
    write(g.as<LineString>(), wkbOrder);
    return;

  case TYPE_POLYGON:
    write(g.as<Polygon>(), wkbOrder);
    return;

  case TYPE_GEOMETRYCOLLECTION:
    write(g.as<GeometryCollection>(), wkbOrder);
    return;

  case TYPE_MULTIPOINT:
    write(g.as<MultiPoint>(), wkbOrder);
    return;

  case TYPE_MULTILINESTRING:
    write(g.as<MultiLineString>(), wkbOrder);
    return;

  case TYPE_MULTIPOLYGON:
    write(g.as<MultiPolygon>(), wkbOrder);
    return;

  case TYPE_TRIANGLE:
    write(g.as<Triangle>(), wkbOrder);
    return;

  case TYPE_TRIANGULATEDSURFACE:
    write(g.as<TriangulatedSurface>(), wkbOrder);
    return;

  case TYPE_POLYHEDRALSURFACE:
    write(g.as<PolyhedralSurface>(), wkbOrder);
    return;

  default:
    std::ostringstream oss;
    oss << "WkbWriter : '" << g.geometryType() << "' is not supported";
    BOOST_THROW_EXCEPTION(std::runtime_error(oss.str()));
  }
}

///
///
///
void
WkbWriter::writeCoordinate(const Point &g, boost::endian::order wkbOrder)
{
  std::array<std::byte, 8> coord;
  // x
  coord = toByte(CGAL::to_double(g.x()), wkbOrder);
  _wkb.insert(_wkb.end(), coord.begin(), coord.end());
  // y
  coord = toByte(CGAL::to_double(g.y()), wkbOrder);
  _wkb.insert(_wkb.end(), coord.begin(), coord.end());
  // z
  if (g.is3D()) {
    coord = toByte(CGAL::to_double(g.z()), wkbOrder);
    _wkb.insert(_wkb.end(), coord.begin(), coord.end());
  }

  // m coordinate
  if (g.isMeasured()) {
    coord = toByte(CGAL::to_double(g.m()), wkbOrder);
    _wkb.insert(_wkb.end(), coord.begin(), coord.end());
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
    const std::array<std::byte, 4> wkbType{toByte(ewkbtype, wkbOrder)};
    _wkb.insert(_wkb.end(), wkbType.begin(), wkbType.end());

    if (_useSrid) {
      const std::array<std::byte, 4> sridByte{toByte(_srid, wkbOrder)};
      _wkb.insert(_wkb.end(), sridByte.begin(), sridByte.end());

      _useSrid = false;
    }
  } else {
    const std::array<std::byte, 4> wkbType{toByte(
        static_cast<uint32_t>(
            g.geometryTypeId() + static_cast<int>(g.is3D()) * COORDINATE_XYZ +
            static_cast<int>(g.isMeasured()) * COORDINATE_XYM),
        wkbOrder)};
    _wkb.insert(_wkb.end(), wkbType.begin(), wkbType.end());
  }
}

///
///
///
void
WkbWriter::write(const Point &g, boost::endian::order wkbOrder)
{
  _wkb.clear();
  writeInner(g, wkbOrder);
}

///
///
///
void
WkbWriter::writeInner(const Point &g, boost::endian::order wkbOrder)
{
  // Endianness
  _wkb.push_back(static_cast<std::byte>(wkbOrder));

  // WkbType
  writeGeometryType(g, wkbOrder);

  if (g.isEmpty()) {
    const std::array<std::byte, 8> nanByte =
        toByte(std::numeric_limits<double>::quiet_NaN(), wkbOrder);
    _wkb.insert(_wkb.end(), nanByte.begin(), nanByte.end());
    _wkb.insert(_wkb.end(), nanByte.begin(), nanByte.end());
    if (g.is3D()) {
      _wkb.insert(_wkb.end(), nanByte.begin(), nanByte.end());
    }
    if (g.isMeasured()) {
      _wkb.insert(_wkb.end(), nanByte.begin(), nanByte.end());
    }
  } else {
    writeCoordinate(g, wkbOrder);
  }
}

///
///
///
void
WkbWriter::write(const LineString &g, boost::endian::order wkbOrder)
{
  _wkb.clear();
  writeInner(g, wkbOrder);
}

///
///
///
void
WkbWriter::writeInner(const LineString &g, boost::endian::order wkbOrder)
{
  // Endianness
  _wkb.push_back(static_cast<std::byte>(wkbOrder));

  // WkbType
  writeGeometryType(g, wkbOrder);

  writeInnerRing(g, wkbOrder);
}
void
WkbWriter::writeInnerRing(const LineString &g, boost::endian::order wkbOrder)
{
  const std::array<std::byte, 4> nbPoints{
      toByte(static_cast<uint32_t>(g.numPoints()), wkbOrder)};
  _wkb.insert(_wkb.end(), nbPoints.begin(), nbPoints.end());

  for (size_t i = 0; i < g.numPoints(); i++) {
    writeCoordinate(g.pointN(i), wkbOrder);
  }
}

///
///
///
void
WkbWriter::write(const Polygon &g, boost::endian::order wkbOrder)
{
  _wkb.clear();
  writeInner(g, wkbOrder);
}

///
///
///
void
WkbWriter::writeInner(const Polygon &g, boost::endian::order wkbOrder)
{
  // Endianness
  _wkb.push_back(static_cast<std::byte>(wkbOrder));

  // WkbType
  writeGeometryType(g, wkbOrder);

  const std::array<std::byte, 4> nbRings{
      toByte(static_cast<uint32_t>(g.numRings()), wkbOrder)};
  _wkb.insert(_wkb.end(), nbRings.begin(), nbRings.end());

  writeInnerRing(g.exteriorRing(), wkbOrder);
  for (size_t i = 0; i < g.numInteriorRings(); i++) {
    writeInnerRing(g.interiorRingN(i), wkbOrder);
  }
}

///
///
///
void
WkbWriter::write(const GeometryCollection &g, boost::endian::order wkbOrder)
{
  _wkb.clear();
  writeInner(g, wkbOrder);
}

void
WkbWriter::writeInner(const GeometryCollection &g,
                      boost::endian::order      wkbOrder)
{
  // Endianness
  _wkb.push_back(static_cast<std::byte>(wkbOrder));

  // WkbType
  writeGeometryType(g, wkbOrder);

  // Number of Geometries
  const std::array<std::byte, 4> numGeometries{
      toByte(static_cast<uint32_t>(g.numGeometries()), wkbOrder)};
  _wkb.insert(_wkb.end(), numGeometries.begin(), numGeometries.end());

  for (size_t i = 0; i < g.numGeometries(); i++) {
    writeRec(g.geometryN(i), wkbOrder);
  }
}
///
///
///
void
WkbWriter::write(const MultiPoint &g, boost::endian::order wkbOrder)
{
  _wkb.clear();
  writeInner(g, wkbOrder);
}

void
WkbWriter::writeInner(const MultiPoint &g, boost::endian::order wkbOrder)
{
  // Endianness
  _wkb.push_back(static_cast<std::byte>(wkbOrder));

  // WkbType
  writeGeometryType(g, wkbOrder);

  // Number of Geometries
  const std::array<std::byte, 4> numGeometries{
      toByte(static_cast<uint32_t>(g.numGeometries()), wkbOrder)};
  _wkb.insert(_wkb.end(), numGeometries.begin(), numGeometries.end());

  for (size_t i = 0; i < g.numGeometries(); i++) {
    writeInner(g.geometryN(i).as<Point>(), wkbOrder);
  }
}
///
///
///
void
WkbWriter::write(const MultiLineString &g, boost::endian::order wkbOrder)
{
  _wkb.clear();
  writeInner(g, wkbOrder);
}

void
WkbWriter::writeInner(const MultiLineString &g, boost::endian::order wkbOrder)
{
  // Endianness
  _wkb.push_back(static_cast<std::byte>(wkbOrder));

  // WkbType
  writeGeometryType(g, wkbOrder);

  // Number of Geometries
  const std::array<std::byte, 4> numGeometries{
      toByte(static_cast<uint32_t>(g.numGeometries()), wkbOrder)};
  _wkb.insert(_wkb.end(), numGeometries.begin(), numGeometries.end());

  for (size_t i = 0; i < g.numGeometries(); i++) {
    writeInner(g.geometryN(i).as<LineString>(), wkbOrder);
  }
}
///
///
///
void
WkbWriter::write(const MultiPolygon &g, boost::endian::order wkbOrder)
{
  _wkb.clear();
  writeInner(g, wkbOrder);
}

void
WkbWriter::writeInner(const MultiPolygon &g, boost::endian::order wkbOrder)
{
  // Endianness
  _wkb.push_back(static_cast<std::byte>(wkbOrder));

  // WkbType
  writeGeometryType(g, wkbOrder);

  // Number of Geometries
  const std::array<std::byte, 4> numGeometries{
      toByte(static_cast<uint32_t>(g.numGeometries()), wkbOrder)};
  _wkb.insert(_wkb.end(), numGeometries.begin(), numGeometries.end());

  for (size_t i = 0; i < g.numGeometries(); i++) {
    writeInner(g.geometryN(i).as<Polygon>(), wkbOrder);
  }
}

///
///
///
void
WkbWriter::write(const Triangle &g, boost::endian::order wkbOrder)
{
  _wkb.clear();
  writeInner(g, wkbOrder);
}

///
///
///
void
WkbWriter::writeInner(const Triangle &g, boost::endian::order wkbOrder)
{
  // Endianness
  _wkb.push_back(static_cast<std::byte>(wkbOrder));

  // WkbType
  writeGeometryType(g, wkbOrder);

  // 4 Points
  if (!g.isEmpty()) {
    // One Ring
    const std::array<std::byte, 4> oneRing =
        toByte(static_cast<uint32_t>(1), wkbOrder);
    _wkb.insert(_wkb.end(), oneRing.begin(), oneRing.end());
    // 4 points
    const std::array<std::byte, 4> fourPoints =
        toByte(static_cast<uint32_t>(4), wkbOrder);
    _wkb.insert(_wkb.end(), fourPoints.begin(), fourPoints.end());
    for (size_t i = 0; i < 4; i++) {
      writeCoordinate(g.vertex(i), wkbOrder);
    }
  }
}

///
///
///
void
WkbWriter::write(const TriangulatedSurface &g, boost::endian::order wkbOrder)
{
  _wkb.clear();
  writeInner(g, wkbOrder);
}

void
WkbWriter::writeInner(const TriangulatedSurface &g,
                      boost::endian::order       wkbOrder)
{
  // Endianness
  _wkb.push_back(static_cast<std::byte>(wkbOrder));

  // WkbType
  writeGeometryType(g, wkbOrder);

  // Number of Geometries
  const std::array<std::byte, 4> numGeometries{
      toByte(static_cast<uint32_t>(g.numGeometries()), wkbOrder)};
  _wkb.insert(_wkb.end(), numGeometries.begin(), numGeometries.end());

  for (size_t i = 0; i < g.numGeometries(); i++) {
    writeInner(g.geometryN(i).as<Triangle>(), wkbOrder);
  }
}
///
///
///
void
WkbWriter::write(const PolyhedralSurface &g, boost::endian::order wkbOrder)
{
  _wkb.clear();
  writeInner(g, wkbOrder);
}

///
///
///
void
WkbWriter::writeInner(const PolyhedralSurface &g, boost::endian::order wkbOrder)
{

  // Endianness
  _wkb.push_back(static_cast<std::byte>(wkbOrder));

  // WkbType
  writeGeometryType(g, wkbOrder);

  // Number of Geometries
  const std::array<std::byte, 4> numGeometries{
      toByte(static_cast<uint32_t>(g.numGeometries()), wkbOrder)};
  _wkb.insert(_wkb.end(), numGeometries.begin(), numGeometries.end());

  for (size_t i = 0; i < g.numGeometries(); i++) {
    writeInner(g.geometryN(i).as<Polygon>(), wkbOrder);
  }
}

} // namespace SFCGAL::detail::io
