// Copyright (c) 2023-2023, Oslandia.
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

//
namespace SFCGAL::detail::io {

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
  readWkb() -> void
  {
    // wkbOrder
    std::byte wkbOrder{read<std::byte>()};
    _swapEndian =
        boost::endian::order::native == boost::endian::order(wkbOrder);

    _geometry = readGeometry();
  }

  auto
  geometry() -> std::unique_ptr<SFCGAL::Geometry>
  {
    return std::move(_geometry);
  }

  auto
  preparedGeometry() -> std::unique_ptr<SFCGAL::PreparedGeometry>
  {
    return std::make_unique<SFCGAL::PreparedGeometry>(std::move(_geometry),
                                                      _srid);
  }

  [[nodiscard]] auto
  srid() const -> srid_t
  {
    return _srid;
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

  // Méthode pour lire une géométrie à partir du WKB
  auto
  readGeometry() -> std::unique_ptr<SFCGAL::Geometry>
  {
    GeometryType geometryType = readGeometryType();
    return readGeometryData(geometryType);
  }

  auto
  readGeometryType() -> GeometryType
  {
    auto geometryType = read<uint32_t>();

    if (_isEWKB || ((geometryType & wkbSRID) == wkbSRID)) {
      if (!_isEWKB) {
        readSRID();
        _isEWKB = true;
      }

      if ((geometryType & wkbZ) == wkbZ) {
        _is3D = true;
      }
      if ((geometryType & wkbM) == wkbM) {
        _isMeasured = true;
      }
      geometryType &= 0x0FFFFFFF;
    } else {
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
      return std::unique_ptr<SFCGAL::Geometry>(
          readInnerMultiGeometries<MultiPoint, Point>().clone());

    case TYPE_MULTILINESTRING:
      return std::unique_ptr<SFCGAL::Geometry>(
          readInnerMultiGeometries<MultiLineString, LineString>().clone());

    case TYPE_MULTIPOLYGON:
      return std::unique_ptr<SFCGAL::Geometry>(
          readInnerMultiGeometries<MultiPolygon, Polygon>().clone());

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
  auto
  readSRID() -> void
  {
    _srid = read<uint32_t>();
  }

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
   * Read MultiGeometries (MultiPoint, MultiLineString, MultiPolygon content
   * from Wkb
   */
  template <typename M, typename G>
  auto
  readInnerMultiGeometries() -> M;

  /**
   * Read Triangle content from wkb
   */
  auto
  readInnerTriangle() -> Triangle;

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

  srid_t _srid = 0;

  bool                              _isEWKB = false;
  std::unique_ptr<SFCGAL::Geometry> _geometry;
};

} // namespace SFCGAL::detail::io

#endif
