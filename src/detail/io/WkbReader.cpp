// Copyright (c) 2023-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/io/WkbReader.h"

#include "SFCGAL/config.h"

#include "SFCGAL/Geometry.h"
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

/**
 * Read Point content from wkb
 */
auto
WkbReader::readInnerPoint() -> Point
{
  double const x{read<double>()};
  double const y{read<double>()};

  if (!(std::isfinite(x) && std::isfinite(y))) {
    return {};
  }
  if (_is3D && _isMeasured) {
    double const z{read<double>()};
    double const m{read<double>()};

    if (!(std::isfinite(z) && std::isfinite(m))) {
      return {};
    }
    return SFCGAL::Point{x, y, z, m};
  }
  if (_is3D) {
    double const z{read<double>()};
    if (!(std::isfinite(z))) {
      return {};
    }
    return SFCGAL::Point{x, y, z};
  }
  if (_isMeasured) {

    double const m{read<double>()};
    if (!(std::isfinite(m))) {
      return {};
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
      SFCGAL::LineString const ls{readInnerLineString()};

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
  try {
    SFCGAL::Polygon poly{readInnerPolygon()};
    if (poly.isEmpty()) {
      return {};
    }

    SFCGAL::LineString geom{poly.exteriorRing()};
    if (geom.isEmpty()) {
      return {};
    }

    return SFCGAL::Triangle{geom.pointN(0), geom.pointN(1), geom.pointN(2)};
  } catch (std::exception &e) {
    std::cerr << e.what();
  }
  return {};
}

/**
 * Read MultiGeometries content from wkb
 */
template <typename M, typename G>
auto
WkbReader::readInnerMultiGeometries() -> M
{
  M result;
  try {
    const uint32_t numGeoms{read<uint32_t>()};
    for (uint32_t i = 0; i < numGeoms; ++i) {
      readWkb();
      G geom{_geometry->template as<G>()};
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
      if (_geometry != nullptr) {
        result.addGeometry(_geometry.release());
      }
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
      if (_geometry != nullptr) {
        SFCGAL::Triangle const geom{_geometry->as<SFCGAL::Triangle>()};
        result.addPatch(geom);
      }
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
      if (_geometry != nullptr) {
        geoms.push_back(_geometry->as<SFCGAL::Polygon>());
      }
    }
  } catch (std::exception &e) {
    std::cerr << e.what();
    return {};
  }
  SFCGAL::PolyhedralSurface const result{geoms};
  return result;
}

auto
WkbReader::readInnerNURBSCurve() -> NURBSCurve
{
  try {
    // Read degree (ISO/SQL-MM format)
    uint32_t degree = read<uint32_t>();

    // Read number of control points
    uint32_t numControlPoints = read<uint32_t>();

    std::vector<Point>          controlPoints;
    std::vector<NURBSCurve::FT> weights;

    controlPoints.reserve(numControlPoints);
    weights.reserve(numControlPoints);

    // Read each control point with ISO WKB standard structure: [byte
    // order][coordinates][flag][weight?]
    for (uint32_t i = 0; i < numControlPoints; i++) {
      // Read byte order for this control point (required by ISO/IEC
      // 13249-3:2016)
      std::byte pointByteOrder = read<std::byte>();
      // Store current swap state and temporarily switch if needed
      bool savedSwapEndian = _swapEndian;
      _swapEndian =
          boost::endian::order::native == boost::endian::order(pointByteOrder);

      // Read coordinates
      double x = read<double>();
      double y = read<double>();

      Point controlPoint;
      if (_is3D && _isMeasured) {
        double z = read<double>();
        double m = read<double>();
        controlPoint =
            Point(NURBSCurve::FT(x), NURBSCurve::FT(y), NURBSCurve::FT(z));
        controlPoint.setM(m);
      } else if (_is3D) {
        double z = read<double>();
        controlPoint =
            Point(NURBSCurve::FT(x), NURBSCurve::FT(y), NURBSCurve::FT(z));
      } else if (_isMeasured) {
        double m     = read<double>();
        controlPoint = Point(NURBSCurve::FT(x), NURBSCurve::FT(y));
        controlPoint.setM(m);
      } else {
        controlPoint = Point(NURBSCurve::FT(x), NURBSCurve::FT(y));
      }

      controlPoints.push_back(controlPoint);

      // Read weight flag
      uint8_t hasWeight = read<uint8_t>();

      // Read weight if flag indicates custom weight
      if (hasWeight != 0) {
        double weight = read<double>();
        weights.push_back(NURBSCurve::FT(weight));
      } else {
        // Default weight is 1.0
        weights.push_back(NURBSCurve::FT(1.0));
      }

      // Restore original swap state
      _swapEndian = savedSwapEndian;
    }

    // Read knot vector
    uint32_t                      numKnots = read<uint32_t>();
    std::vector<NURBSCurve::Knot> knots;
    knots.reserve(numKnots);

    for (uint32_t i = 0; i < numKnots; i++) {
      double knot = read<double>();
      knots.push_back(NURBSCurve::FT(knot));
    }

    // Create the NURBS curve
    NURBSCurve curve(controlPoints, weights, degree, knots);

    return curve;

  } catch (std::exception &e) {
    std::cerr << "WkbReader::readInnerNURBSCurve error: " << e.what()
              << std::endl;
    return NURBSCurve{};
  }
}

} // namespace SFCGAL::detail::io
