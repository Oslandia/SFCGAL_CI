// Copyright (c) 2023-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/io/WkbReader.h"

#include "SFCGAL/config.h"

#include "SFCGAL/Exception.h"
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
    BOOST_THROW_EXCEPTION(
        Exception((boost::format("WkbReader error: %s") % e.what()).str()));
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
    BOOST_THROW_EXCEPTION(
        Exception((boost::format("WkbReader error: %s") % e.what()).str()));
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
    BOOST_THROW_EXCEPTION(
        Exception((boost::format("WkbReader error: %s") % e.what()).str()));
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
    BOOST_THROW_EXCEPTION(
        Exception((boost::format("WkbReader error: %s") % e.what()).str()));
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
    BOOST_THROW_EXCEPTION(
        Exception((boost::format("WkbReader error: %s") % e.what()).str()));
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
    BOOST_THROW_EXCEPTION(
        Exception((boost::format("WkbReader error: %s") % e.what()).str()));
  }
  SFCGAL::PolyhedralSurface const result{geoms};
  return result;
}

auto
WkbReader::readInnerNURBSCurve() -> NURBSCurve
{
  try {
    // Read degree (ISO/SQL-MM format)
    auto degree = read<uint32_t>();

    // Validate degree
    if (degree == 0 || degree > 100) {
      BOOST_THROW_EXCEPTION(Exception(
          (boost::format("Invalid NURBS degree: %d (must be > 0 and <= 100)") %
           degree)
              .str()));
    }

    // Read number of control points
    auto numControlPoints = read<uint32_t>();

    // Validate control points count
    if (numControlPoints < degree + 1) {
      BOOST_THROW_EXCEPTION(
          Exception((boost::format("Invalid control points count: %d (must be "
                                   ">= degree + 1 = %d)") %
                     numControlPoints % (degree + 1))
                        .str()));
    }

    std::vector<Point>          controlPoints;
    std::vector<NURBSCurve::FT> weights;

    controlPoints.reserve(numControlPoints);
    weights.reserve(numControlPoints);

    // Read each control point with ISO WKB standard structure: [byte
    // order][coordinates][flag][weight?]
    for (uint32_t i = 0; i < numControlPoints; i++) {
      // Read byte order for this control point (required by ISO/IEC
      // 13249-3:2016)
      auto pointByteOrder = read<std::byte>();

      // Store current swap state and temporarily switch if needed
      // Use RAII to ensure swap state is always restored
      bool savedSwapEndian = _swapEndian;
      _swapEndian =
          boost::endian::order::native == boost::endian::order(pointByteOrder);

      // RAII guard to restore swap state in case of exception
      struct SwapEndianGuard {
        bool &swapEndian;
        bool  originalValue;
        SwapEndianGuard(bool &swapEndianRef, bool orig)
            : swapEndian(swapEndianRef), originalValue(orig)
        {
        }
        ~SwapEndianGuard() { swapEndian = originalValue; }
      } guard(_swapEndian, savedSwapEndian);

      // Validate coordinate values before using them
      auto x = read<double>();
      auto y = read<double>();

      if (!std::isfinite(x) || !std::isfinite(y)) {
        BOOST_THROW_EXCEPTION(Exception(
            (boost::format(
                 "Invalid coordinates at control point %d: x=%f, y=%f") %
             i % x % y)
                .str()));
      }

      Point controlPoint;
      if (_is3D && _isMeasured) {
        auto z = read<double>();
        auto m = read<double>();
        if (!std::isfinite(z) || !std::isfinite(m)) {
          BOOST_THROW_EXCEPTION(
              Exception((boost::format("Invalid z or m coordinate at control "
                                       "point %d: z=%f, m=%f") %
                         i % z % m)
                            .str()));
        }
        controlPoint =
            Point(NURBSCurve::FT(x), NURBSCurve::FT(y), NURBSCurve::FT(z));
        controlPoint.setM(m);
      } else if (_is3D) {
        auto z = read<double>();
        if (!std::isfinite(z)) {
          BOOST_THROW_EXCEPTION(Exception(
              (boost::format("Invalid z coordinate at control point %d: z=%f") %
               i % z)
                  .str()));
        }
        controlPoint =
            Point(NURBSCurve::FT(x), NURBSCurve::FT(y), NURBSCurve::FT(z));
      } else if (_isMeasured) {
        auto m = read<double>();
        if (!std::isfinite(m)) {
          BOOST_THROW_EXCEPTION(Exception(
              (boost::format("Invalid m coordinate at control point %d: m=%f") %
               i % m)
                  .str()));
        }
        controlPoint = Point(NURBSCurve::FT(x), NURBSCurve::FT(y));
        controlPoint.setM(m);
      } else {
        controlPoint = Point(NURBSCurve::FT(x), NURBSCurve::FT(y));
      }

      controlPoints.push_back(controlPoint);

      // Read weight flag
      auto hasWeight = read<uint8_t>();

      // Read weight if flag indicates custom weight
      if (hasWeight != 0) {
        auto weight = read<double>();
        if (!std::isfinite(weight) || weight <= 0.0) {
          BOOST_THROW_EXCEPTION(
              Exception((boost::format("Invalid weight at control point %d: %f "
                                       "(must be finite and positive)") %
                         i % weight)
                            .str()));
        }
        weights.emplace_back(weight);
      } else {
        // Default weight is 1.0
        weights.emplace_back(1.0);
      }

      // Swap state will be restored by the RAII guard destructor
    }

    // Read knot vector
    auto numKnots = read<uint32_t>();

    // Validate knot vector size
    uint32_t expectedKnots = numControlPoints + degree + 1;
    if (numKnots != expectedKnots) {
      BOOST_THROW_EXCEPTION(Exception(
          (boost::format(
               "Invalid knot vector size: %d (expected %d = %d + %d + 1)") %
           numKnots % expectedKnots % numControlPoints % degree)
              .str()));
    }

    std::vector<NURBSCurve::Knot> knots;
    knots.reserve(numKnots);

    for (uint32_t i = 0; i < numKnots; i++) {
      auto knot = read<double>();
      if (!std::isfinite(knot)) {
        BOOST_THROW_EXCEPTION(Exception(
            (boost::format(
                 "Invalid knot value at index %d: %f (must be finite)") %
             i % knot)
                .str()));
      }
      knots.emplace_back(knot);
    }

    // Validate that knots are non-decreasing
    for (uint32_t i = 1; i < numKnots; i++) {
      if (knots[i] < knots[i - 1]) {
        BOOST_THROW_EXCEPTION(
            Exception((boost::format("Knot vector not non-decreasing: knot[%d] "
                                     "= %f < knot[%d] = %f") %
                       i % CGAL::to_double(knots[i]) % (i - 1) %
                       CGAL::to_double(knots[i - 1]))
                          .str()));
      }
    }

    // Create the NURBS curve
    NURBSCurve curve(controlPoints, weights, degree, knots);

    return curve;

  } catch (std::exception &e) {
    BOOST_THROW_EXCEPTION(Exception(
        (boost::format("WkbReader::readInnerNURBSCurve error: %s") % e.what())
            .str()));
  }
}

} // namespace SFCGAL::detail::io
