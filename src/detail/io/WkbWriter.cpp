// Copyright (c) 2023-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <array>
#include <boost/endian/conversion.hpp> // don't use bit, since it requires c++20
#include <cstddef>

#include "SFCGAL/detail/io/WkbWriter.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"

#include <boost/exception/all.hpp>
#include <exception>

namespace SFCGAL::detail::io {

void
WkbWriter::writeRec(const Geometry &geometry, boost::endian::order wkbOrder)
{
  switch (geometry.geometryTypeId()) {
  case TYPE_POINT:
    writeInner(geometry.as<Point>(), wkbOrder);
    return;

  case TYPE_LINESTRING:
    writeInner(geometry.as<LineString>(), wkbOrder);
    return;

  case TYPE_POLYGON:
    writeInner(geometry.as<Polygon>(), wkbOrder);
    return;

  case TYPE_GEOMETRYCOLLECTION:
    writeInner(geometry.as<GeometryCollection>(), wkbOrder);
    return;

  case TYPE_MULTIPOINT:
    writeInner<MultiPoint, Point>(geometry.as<MultiPoint>(), wkbOrder);
    return;

  case TYPE_MULTILINESTRING:
    writeInner<MultiLineString, LineString>(geometry.as<MultiLineString>(),
                                            wkbOrder);
    return;

  case TYPE_MULTIPOLYGON:
    writeInner<MultiPolygon, Polygon>(geometry.as<MultiPolygon>(), wkbOrder);
    return;

  case TYPE_TRIANGLE:
    writeInner(geometry.as<Triangle>(), wkbOrder);
    return;

  case TYPE_TRIANGULATEDSURFACE:
    writeInner(geometry.as<TriangulatedSurface>(), wkbOrder);
    return;

  case TYPE_POLYHEDRALSURFACE:
    writeInner(geometry.as<PolyhedralSurface>(), wkbOrder);
    return;

  case TYPE_NURBSCURVE:
    writeInner(geometry.as<NURBSCurve>(), wkbOrder);
    return;

  case TYPE_SOLID:
    writeInner<Solid, PolyhedralSurface>(geometry.as<Solid>(), wkbOrder);
    return;

  default:
    std::ostringstream oss;
    oss << "WkbWriter: type '" << geometry.geometryType()
        << "' is not supported";
    BOOST_THROW_EXCEPTION(InappropriateGeometryException(oss.str()));
  }
}

void
WkbWriter::write(const Geometry &geometry, const srid_t &srid,
                 boost::endian::order wkbOrder)
{

  _useSrid = true;
  _isEWKB  = true;
  _srid    = srid;

  write(geometry, wkbOrder);
}

void
WkbWriter::write(const Geometry &geometry, boost::endian::order wkbOrder)
{
  writeRec(geometry, wkbOrder);
}

void
WkbWriter::writeCoordinate(const Point &geometry, boost::endian::order wkbOrder)
{
  // x
  toByte(CGAL::to_double(geometry.x()), wkbOrder);

  // y
  toByte(CGAL::to_double(geometry.y()), wkbOrder);

  // z
  if (geometry.is3D()) {
    toByte(CGAL::to_double(geometry.z()), wkbOrder);
  }

  // m coordinate
  if (geometry.isMeasured()) {
    toByte(CGAL::to_double(geometry.m()), wkbOrder);
  }
}

void
WkbWriter::writeGeometryType(const Geometry      &geometry,
                             boost::endian::order wkbOrder)
{

  if (_isEWKB) {

    uint32_t ewkbtype = geometry.geometryTypeId();

    if (geometry.is3D()) {
      ewkbtype |= wkbZ;
    }
    if (geometry.isMeasured()) {
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
    toByte(static_cast<uint32_t>(
               geometry.geometryTypeId() +
               static_cast<int>(geometry.is3D()) * COORDINATE_XYZ +
               static_cast<int>(geometry.isMeasured()) * COORDINATE_XYM),
           wkbOrder);
  }
}

void
WkbWriter::writeInner(const Point &geometry, boost::endian::order wkbOrder)
{
  // Endianness
  toStream(std::array<std::byte, 1>{static_cast<std::byte>(wkbOrder)});

  // WkbType
  writeGeometryType(geometry, wkbOrder);

  if (geometry.isEmpty()) {
    toByte(std::numeric_limits<double>::quiet_NaN(), wkbOrder);
    toByte(std::numeric_limits<double>::quiet_NaN(), wkbOrder);
    if (geometry.is3D()) {
      toByte(std::numeric_limits<double>::quiet_NaN(), wkbOrder);
    }
    if (geometry.isMeasured()) {
      toByte(std::numeric_limits<double>::quiet_NaN(), wkbOrder);
    }
  } else {
    writeCoordinate(geometry, wkbOrder);
  }
}

void
WkbWriter::writeInner(const LineString &geometry, boost::endian::order wkbOrder)
{
  // Endianness
  toStream(std::array<std::byte, 1>{static_cast<std::byte>(wkbOrder)});

  // WkbType
  writeGeometryType(geometry, wkbOrder);

  writeInnerRing(geometry, wkbOrder);
}
void
WkbWriter::writeInnerRing(const LineString    &geometry,
                          boost::endian::order wkbOrder)
{
  toByte(static_cast<uint32_t>(geometry.numPoints()), wkbOrder);

  for (size_t i = 0; i < geometry.numPoints(); i++) {
    writeCoordinate(geometry.pointN(i), wkbOrder);
  }
}

void
WkbWriter::writeInner(const Polygon &geometry, boost::endian::order wkbOrder)
{
  // Endianness
  toStream(std::array<std::byte, 1>{static_cast<std::byte>(wkbOrder)});

  // WkbType
  writeGeometryType(geometry, wkbOrder);

  toByte(static_cast<uint32_t>(geometry.numRings()), wkbOrder);

  writeInnerRing(geometry.exteriorRing(), wkbOrder);
  for (size_t i = 0; i < geometry.numInteriorRings(); i++) {
    writeInnerRing(geometry.interiorRingN(i), wkbOrder);
  }
}

void
WkbWriter::writeInner(const Triangle &geometry, boost::endian::order wkbOrder)
{
  // Endianness
  toStream(std::array<std::byte, 1>{static_cast<std::byte>(wkbOrder)});

  // WkbType
  writeGeometryType(geometry, wkbOrder);

  if (!geometry.isEmpty()) {
    // One Ring
    toByte(static_cast<uint32_t>(1), wkbOrder);
    // 4 points
    toByte(static_cast<uint32_t>(4), wkbOrder);
    for (int i = 0; i < 4; i++) {
      writeCoordinate(geometry.vertex(i), wkbOrder);
    }
  } else {
    // Empty triangle has 0 rings
    toByte(static_cast<uint32_t>(0), wkbOrder);
  }
}

void
WkbWriter::writeInner(const GeometryCollection &geometry,
                      boost::endian::order      wkbOrder)
{
  // Endianness
  toStream(std::array<std::byte, 1>{static_cast<std::byte>(wkbOrder)});

  // WkbType
  writeGeometryType(geometry, wkbOrder);

  // Number of Geometries
  toByte(static_cast<uint32_t>(geometry.numGeometries()), wkbOrder);

  for (size_t i = 0; i < geometry.numGeometries(); i++) {
    writeRec(geometry.geometryN(i), wkbOrder);
  }
}

void
WkbWriter::writeInner(const PolyhedralSurface &polyhedralSurface,
                      boost::endian::order     wkbOrder)
{
  // Endianness
  toStream(std::array<std::byte, 1>{static_cast<std::byte>(wkbOrder)});

  // WkbType
  writeGeometryType(polyhedralSurface, wkbOrder);

  // Number of Polygons
  toByte(static_cast<uint32_t>(polyhedralSurface.numPatches()), wkbOrder);

  for (size_t i = 0; i < polyhedralSurface.numPatches(); i++) {
    writeRec(polyhedralSurface.patchN(i), wkbOrder);
  }
}

void
WkbWriter::writeInner(const TriangulatedSurface &triangulatedSurface,
                      boost::endian::order       wkbOrder)
{
  // Endianness
  toStream(std::array<std::byte, 1>{static_cast<std::byte>(wkbOrder)});

  // WkbType
  writeGeometryType(triangulatedSurface, wkbOrder);

  // Number of Triangles
  toByte(static_cast<uint32_t>(triangulatedSurface.numPatches()), wkbOrder);

  for (size_t i = 0; i < triangulatedSurface.numPatches(); i++) {
    writeRec(triangulatedSurface.patchN(i), wkbOrder);
  }
}

void
WkbWriter::writeInner(const NURBSCurve &geometry, boost::endian::order wkbOrder)
{
  // Endianness
  toStream(std::array<std::byte, 1>{static_cast<std::byte>(wkbOrder)});

  // WkbType
  writeGeometryType(geometry, wkbOrder);

  // Degree (WKB standard format)
  toByte(static_cast<uint32_t>(geometry.degree()), wkbOrder);

  // Number of control points
  toByte(static_cast<uint32_t>(geometry.numControlPoints()), wkbOrder);

  // Write each control point with ISO WKB standard structure: [byte
  // order][coordinates][flag][weight?]
  for (size_t i = 0; i < geometry.numControlPoints(); i++) {
    auto weight          = geometry.weight(i);
    bool hasCustomWeight = (std::abs(CGAL::to_double(weight) - 1.0) > 1e-10);

    // Write byte order for this control point (required by ISO/IEC
    // 13249-3:2016)
    toStream(std::array<std::byte, 1>{static_cast<std::byte>(wkbOrder)});

    // Write coordinates
    writeCoordinate(geometry.controlPointN(i), wkbOrder);

    // Write weight flag (0 = default weight, 1 = custom weight follows)
    toByte(static_cast<uint8_t>(hasCustomWeight ? 1 : 0), wkbOrder);

    // Write custom weight if present
    if (hasCustomWeight) {
      toByte(CGAL::to_double(weight), wkbOrder);
    }
  }

  // Knot vector
  const auto &knots = geometry.knotVector();
  toByte(static_cast<uint32_t>(knots.size()), wkbOrder);

  for (const auto &knot : knots) {
    toByte(CGAL::to_double(knot), wkbOrder);
  }
}

template <typename M, typename G>
void
WkbWriter::writeInner(const M &geometry, boost::endian::order wkbOrder)
{
  // Endianness
  toStream(std::array<std::byte, 1>{static_cast<std::byte>(wkbOrder)});

  // WkbType
  writeGeometryType(geometry, wkbOrder);

  // Number of Geometries/Shells
  if constexpr (std::is_same_v<M, Solid>) {
    toByte(static_cast<uint32_t>(geometry.numShells()), wkbOrder);
    for (size_t i = 0; i < geometry.numShells(); i++) {
      writeRec(geometry.shellN(i).template as<G>(), wkbOrder);
    }
  } else {
    toByte(static_cast<uint32_t>(geometry.numGeometries()), wkbOrder);

    for (size_t i = 0; i < geometry.numGeometries(); i++) {
      writeRec(geometry.geometryN(i).template as<G>(), wkbOrder);
    }
  }
}

} // namespace SFCGAL::detail::io
