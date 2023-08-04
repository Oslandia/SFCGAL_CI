// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef _SFCGAL_IO_WKBREADER_H_
#define _SFCGAL_IO_WKBREADER_H_

#include <sstream>

#include <SFCGAL/config.h>

#include <SFCGAL/Geometry.h>
#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/PreparedGeometry.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>

#include <SFCGAL/detail/tools/InputStreamReader.h>

#include <boost/endian/conversion.hpp> // don't use bit, since it requires c++20
#include <sys/endian.h>
//
namespace SFCGAL {
namespace detail {
namespace io {

/**
 * read WKB geometry
 *
 */
class SFCGAL_API WkbReader {
public:
  /**
   * read WKB from input stream
   */
  explicit WkbReader(const std::string &wkbHexString) : _wkbData(wkbHexString)
  {
  }

  auto
  readWkb() -> std::unique_ptr<SFCGAL::Geometry>
  {
    // wkbOrder
    std::byte wkbOrder{read<std::byte>()};
    _swapEndian =
        boost::endian::order::native == boost::endian::order(wkbOrder);

    return readGeometry();
  }
  // Méthode pour lire une géométrie à partir du WKB
  auto
  readGeometry() -> std::unique_ptr<SFCGAL::Geometry>
  {
    GeometryType geometryType = readGeometryType();
    return readGeometryData(geometryType);
  }

private:
  template <typename T>
  auto
  read() -> T
  {

    const size_t      nbElements = 2;
    const size_t      sizeType   = sizeof(T);
    const std::string s    = _wkbData.substr(_index, nbElements * sizeType);
    const int         base = 16;
    union {
      std::array<std::byte, sizeType> byteArray;
      T                               d;
    };

    for (size_t i = 0; i < sizeType; i++) {
      size_t      chunkPos = nbElements * i;
      std::string byteStr  = s.substr(chunkPos, nbElements);
      byteArray[i] = static_cast<std::byte>(std::stoi(byteStr, nullptr, base));
    }

    _index += sizeType * nbElements;
    return d;
  }

  auto
  readGeometryType() -> GeometryType
  {
    uint32_t geometryType = read<uint32_t>();
    if (geometryType >= COORDINATE_XYZM) {
      _is3D       = true;
      _isMeasured = true;
      geometryType -= COORDINATE_XYZM;
    } else if (geometryType >= COORDINATE_XYM) {
      _is3D       = false;
      _isMeasured = true;
      geometryType -= COORDINATE_XYM;
    } else if (geometryType >= COORDINATE_XYZ) {
      _is3D       = true;
      _isMeasured = false;
      geometryType -= COORDINATE_XYZ;
    }
    return static_cast<GeometryType>(geometryType);
  }

  auto
  readGeometryData(GeometryType geometryType) -> std::unique_ptr<Geometry>
  {
    switch (geometryType) {
    case TYPE_POINT:
      return std::unique_ptr<SFCGAL::Geometry>{readInnerPoint().clone()};

    case TYPE_LINESTRING:
      return std::unique_ptr<SFCGAL::Geometry>(readInnerLineString().clone());

    case TYPE_POLYGON:
      return std::unique_ptr<SFCGAL::Geometry>(readInnerPolygon().clone());

    case TYPE_GEOMETRYCOLLECTION:
      return std::unique_ptr<SFCGAL::Geometry>(
          readInnerGeometryCollection().clone());

    case TYPE_MULTIPOINT:
      return std::unique_ptr<SFCGAL::Geometry>(readInnerMultiPoint().clone());

    case TYPE_MULTILINESTRING:
      return std::unique_ptr<SFCGAL::Geometry>(
          readInnerMultiLineString().clone());

    case TYPE_MULTIPOLYGON:
      return std::unique_ptr<SFCGAL::Geometry>(readInnerMultiPolygon().clone());

    case TYPE_TRIANGLE:
      return std::unique_ptr<SFCGAL::Geometry>(readInnerTriangle().clone());

    case TYPE_TRIANGULATEDSURFACE:
      return std::unique_ptr<SFCGAL::Geometry>(
          readInnerTriangulatedSurface().clone());

    case TYPE_POLYHEDRALSURFACE:
      return std::unique_ptr<SFCGAL::Geometry>(
          readInnerPolyhedralSurface().clone());

    default:
      std::ostringstream oss;
      oss << "WkbWriter : '" << geometryType << "' is not supported";
      std::cerr << oss.str() << std::endl;

      return std::unique_ptr<SFCGAL::Geometry>();
      // BOOST_THROW_EXCEPTION(std::runtime_error(oss.str()));
    }
  }
  /**
   * read an SRID, if present
   *
   */
  srid_t
  readSRID();

  /**
   * Read Point content from wkb
   */
  auto
  readInnerPoint() -> Point;

  /**
   * Read LineString content from wkb
   */
  auto
  readInnerLineString() -> LineString;

  /**
   * Read Polygon content from wkb
   */
  auto
  readInnerPolygon() -> Polygon;

  /**
   * Read Triangle content from wkb
   */
  auto
  readInnerTriangle() -> Triangle;

  /**
   * Read MultiPoint content from wkb
   */
  auto
  readInnerMultiPoint() -> MultiPoint;

  /**
   * Read MultiLineString content from wkb
   */
  auto
  readInnerMultiLineString() -> MultiLineString;

  /**
   * Read MultiPolygon content from wkb
   */
  auto
  readInnerMultiPolygon() -> MultiPolygon;

  /**
   * Read GeometryCollection content from wkb
   */
  auto
  readInnerGeometryCollection() -> GeometryCollection;

  /**
   * Read TriangulatedSurface content from wkb
   */
  auto
  readInnerTriangulatedSurface() -> TriangulatedSurface;

  /**
   * Read PolyhedralSurface content from wkb
   */
  auto
  readInnerPolyhedralSurface() -> PolyhedralSurface;

  /**
   * actually reading 3D ?
   */
  bool _is3D = false;
  /**
   * actually reading Measured ?
   */
  bool _isMeasured = false;

  /**
   * wkb data
   */
  std::string _wkbData;

  bool _swapEndian = false;

  size_t _index = 0;
};

} // namespace io
} // namespace detail
} // namespace SFCGAL

#endif
