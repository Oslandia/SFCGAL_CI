// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/detail/io/WkbReader.h>

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
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>
#include <boost/exception/all.hpp>
#include <exception>

namespace SFCGAL::detail::io {

/**
 * Read Point content from wkb
 */
auto
WkbReader::readInnerPoint() -> Point
{
  double x{read<double>()};
  double y{read<double>()};

  if (!(std::isfinite(x) && std::isfinite(y))) {
    return {};
  }
  if (_is3D && _isMeasured) {
    double z{read<double>()};
    double m{read<double>()};

    if (!(std::isfinite(z) && std::isfinite(m))) {
      return {};
    }
    return SFCGAL::Point(x, y, z, m);
  }
  if (_is3D) {
    double z{read<double>()};
    if (!(std::isfinite(z))) {
      return {};
    }
    return SFCGAL::Point(x, y, z);
  }
  if (_isMeasured) {

    double m{read<double>()};
    if (!(std::isfinite(m))) {
      return SFCGAL::Point();
    }
    SFCGAL::Point result{x, y};
    result.setM(m);
    return result;
  }
  return SFCGAL::Point{x, y};
}

/**
 * Read LineString content from wkb
 */
auto
WkbReader::readInnerLineString() -> LineString
{
  SFCGAL::LineString result;
  try {
    const uint32_t numPoints{read<uint32_t>()};
    for (uint32_t i = 0; i < numPoints; ++i) {
      result.addPoint(readInnerPoint());
    }
  } catch (std::exception &e) {
    std::cerr << e.what();
    return {};
  }
  return result;
}

/**
 * Read Polygon content from wkb
 */
auto
WkbReader::readInnerPolygon() -> Polygon
{
  SFCGAL::Polygon result;
  try {
    const uint32_t numRings{read<uint32_t>()};
    for (uint32_t i = 0; i < numRings; ++i) {
      SFCGAL::LineString ls{readInnerLineString()};

      if (i == 0) {
        result.setExteriorRing(ls);
      } else {
        result.addInteriorRing(ls);
      }
    }
  } catch (std::exception &e) {
    std::cerr << e.what();
    return {};
  }
  return result;
}

/**
 * Read Triangle content from wkb
 */
auto
WkbReader::readInnerTriangle() -> Triangle
{
  if (_index >= _wkbData.size()) {
    return {};
  }
  try {
    SFCGAL::LineString geom{readInnerPolygon().exteriorRing()};

    return SFCGAL::Triangle{geom.pointN(0), geom.pointN(1), geom.pointN(2)};
  } catch (std::exception &e) {
    std::cerr << e.what();
  }
  return {};
}

/**
 * Read MultiPoint content from wkb
 */
auto
WkbReader::readInnerMultiPoint() -> MultiPoint
{
  SFCGAL::MultiPoint result;
  try {
    const uint32_t numGeoms{read<uint32_t>()};
    for (uint32_t i = 0; i < numGeoms; ++i) {
      readWkb();
      SFCGAL::Point geom{_geometry->as<SFCGAL::Point>()};
      result.addGeometry(geom);
    }
  } catch (std::exception &e) {
    std::cerr << e.what();
    return {};
  }
  return result;
}

/**
 * Read MultiLineString content from wkb
 */
auto
WkbReader::readInnerMultiLineString() -> MultiLineString
{
  SFCGAL::MultiLineString result;
  try {
    const uint32_t numGeoms{read<uint32_t>()};
    for (uint32_t i = 0; i < numGeoms; ++i) {
      readWkb();
      SFCGAL::LineString geom{_geometry->as<SFCGAL::LineString>()};
      result.addGeometry(geom);
    }
  } catch (std::exception &e) {
    std::cerr << e.what();
    return {};
  }
  return result;
}

/**
 * Read MultiPolygon content from wkb
 */
auto
WkbReader::readInnerMultiPolygon() -> MultiPolygon
{
  SFCGAL::MultiPolygon result;
  try {
    const uint32_t numGeoms{read<uint32_t>()};
    for (uint32_t i = 0; i < numGeoms; ++i) {
      readWkb();
      SFCGAL::Polygon geom{_geometry->as<SFCGAL::Polygon>()};
      result.addGeometry(geom);
    }
  } catch (std::exception &e) {
    std::cerr << e.what();
    return {};
  }
  return result;
}

/**
 * Read GeometryCollection content from wkb
 */
auto
WkbReader::readInnerGeometryCollection() -> GeometryCollection
{
  SFCGAL::GeometryCollection result;
  try {
    const uint32_t numGeoms{read<uint32_t>()};
    for (uint32_t i = 0; i < numGeoms; ++i) {
      readWkb();
      result.addGeometry(_geometry.release());
    }
  } catch (std::exception &e) {
    std::cerr << e.what();
    return {};
  }
  return result;
}

/**
 * Read TriangulatedSurface content from wkb
 */
auto
WkbReader::readInnerTriangulatedSurface() -> TriangulatedSurface
{
  SFCGAL::TriangulatedSurface result;
  try {
    const uint32_t numGeoms{read<uint32_t>()};
    for (uint32_t i = 0; i < numGeoms; ++i) {
      readWkb();
      SFCGAL::Triangle geom{_geometry->as<SFCGAL::Triangle>()};
      result.addTriangle(geom);
    }
  } catch (std::exception &e) {
    std::cerr << e.what();
    return {};
  }
  return result;
}

/**
 * Read PolyhedralSurface content from wkb
 */
auto
WkbReader::readInnerPolyhedralSurface() -> PolyhedralSurface
{
  std::vector<Polygon> geoms;
  try {
    const uint32_t numGeoms{read<uint32_t>()};
    for (uint32_t i = 0; i < numGeoms; ++i) {
      readWkb();
      geoms.push_back(_geometry->as<SFCGAL::Polygon>());
    }
  } catch (std::exception &e) {
    std::cerr << e.what();
    return {};
  }
  SFCGAL::PolyhedralSurface result{geoms};
  return result;
}

} // namespace SFCGAL::detail::io
