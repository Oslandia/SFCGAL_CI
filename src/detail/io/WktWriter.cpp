// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/io/WktWriter.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
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

namespace impl {
auto
writeFT(std::ostream &outStream, const CGAL::Gmpq &fraction) -> std::ostream &
{
  outStream << fraction;
  return outStream;
}

#ifdef CGAL_USE_GMPXX
auto
writeFT(std::ostream &outStream, const mpq_class &fraction) -> std::ostream &
{
  outStream << fraction.get_num() << "/" << fraction.get_den();
  return outStream;
}
#endif
} // namespace impl

WktWriter::WktWriter(std::ostream &outStream) : _s(outStream) {}

void
WktWriter::writeRec(const Geometry &geometry)
{
  switch (geometry.geometryTypeId()) {
  case TYPE_POINT:
    write(geometry.as<Point>());
    return;

  case TYPE_LINESTRING:
    write(geometry.as<LineString>());
    return;

  case TYPE_POLYGON:
    write(geometry.as<Polygon>());
    return;

  case TYPE_GEOMETRYCOLLECTION:
    write(geometry.as<GeometryCollection>());
    return;

  case TYPE_MULTIPOINT:
    write(geometry.as<MultiPoint>());
    return;

  case TYPE_MULTILINESTRING:
    write(geometry.as<MultiLineString>());
    return;

  case TYPE_MULTIPOLYGON:
    write(geometry.as<MultiPolygon>());
    return;

  case TYPE_TRIANGLE:
    write(geometry.as<Triangle>());
    return;

  case TYPE_TRIANGULATEDSURFACE:
    write(geometry.as<TriangulatedSurface>());
    return;

  case TYPE_POLYHEDRALSURFACE:
    write(geometry.as<PolyhedralSurface>());
    return;

  case TYPE_SOLID:
    write(geometry.as<Solid>());
    return;

  case TYPE_MULTISOLID:
    write(geometry.as<MultiSolid>());
    return;

  case TYPE_NURBSCURVE:
    write(geometry.as<NURBSCurve>());
    return;
  }

  std::ostringstream oss;
  oss << "WktWriter : '" << geometry.geometryType() << "' is not supported";
  BOOST_THROW_EXCEPTION(InappropriateGeometryException(oss.str()));
}

void
WktWriter::write(const Geometry &geometry, bool exact)
{
  _exactWrite = exact;
  writeRec(geometry);
}

void
WktWriter::writeCoordinateType(const Geometry &geometry)
{
  if (geometry.is3D() && !geometry.isMeasured()) {
    _s << "Z ";
  } else if (!geometry.is3D() && geometry.isMeasured()) {
    _s << "M ";
  } else if (geometry.is3D() && geometry.isMeasured()) {
    _s << "ZM ";
  }
}

static auto
fixZeroNeg(double val, int precision) -> double
{
  if (std::abs(val) < std::pow(10, -precision)) {
    return 0;
  }
  return val;
}

static auto
fixZeroNegForWeights(double val, int precision) -> double
{
  // NURBS weights must be positive
  if (val <= 0) {
    // Convert non-positive weights to minimum positive value
    return std::pow(10, -precision);
  }

  // For positive weights, preserve small values instead of converting to 0
  if (val < std::pow(10, -precision)) {
    return std::pow(10, -precision);
  }

  return val;
}

void
WktWriter::writeCoordinate(const Point &point)
{
  if (_exactWrite) {
    impl::writeFT(_s, CGAL::exact(point.x())) << " ";
    impl::writeFT(_s, CGAL::exact(point.y()));

    if (point.is3D()) {
      _s << " ";
      impl::writeFT(_s, CGAL::exact(point.z()));
    }
  } else {
    _s << fixZeroNeg(CGAL::to_double(point.x()),
                     static_cast<int>(_s.precision()))
       << " "
       << fixZeroNeg(CGAL::to_double(point.y()),
                     static_cast<int>(_s.precision()));

    if (point.is3D()) {
      _s << " "
         << fixZeroNeg(CGAL::to_double(point.z()),
                       static_cast<int>(_s.precision()));
    }
  }

  // m coordinate
  if (point.isMeasured()) {
    _s << " "
       << fixZeroNeg(CGAL::to_double(point.m()),
                     static_cast<int>(_s.precision()));
  }
}

void
WktWriter::write(const Point &point)
{
  _s << "POINT ";
  writeCoordinateType(point);

  if (point.isEmpty()) {
    _s << "EMPTY";
    return;
  }

  writeInner(point);
}

void
WktWriter::writeInner(const Point &point)
{
  if (point.isEmpty()) {
    _s << "EMPTY";
    return;
  }

  _s << "(";
  writeCoordinate(point);
  _s << ")";
}

void
WktWriter::write(const LineString &lineString)
{
  _s << "LINESTRING ";
  writeCoordinateType(lineString);

  if (lineString.isEmpty()) {
    _s << "EMPTY";
    return;
  }

  writeInner(lineString);
}

void
WktWriter::writeInner(const LineString &lineString)
{
  _s << "(";

  for (size_t i = 0; i < lineString.numPoints(); i++) {
    if (i != 0) {
      _s << ",";
    }

    writeCoordinate(lineString.pointN(i));
  }

  _s << ")";
}

void
WktWriter::write(const Polygon &polygon)
{
  _s << "POLYGON ";
  writeCoordinateType(polygon);

  if (polygon.isEmpty()) {
    _s << "EMPTY";
    return;
  }

  writeInner(polygon);
}

void
WktWriter::writeInner(const Polygon &polygon)
{
  _s << "(";
  writeInner(polygon.exteriorRing());

  for (size_t i = 0; i < polygon.numInteriorRings(); i++) {
    _s << ",";
    writeInner(polygon.interiorRingN(i));
  }

  _s << ")";
}

void
WktWriter::write(const GeometryCollection &collection)
{
  _s << "GEOMETRYCOLLECTION ";
  writeCoordinateType(collection);

  if (collection.isEmpty()) {
    _s << "EMPTY";
    return;
  }

  _s << "(";

  for (size_t i = 0; i < collection.numGeometries(); i++) {
    if (i != 0) {
      _s << ",";
    }

    writeRec(collection.geometryN(i));
  }

  _s << ")";
}

void
WktWriter::write(const MultiPoint &multiPoint)
{
  _s << "MULTIPOINT ";
  writeCoordinateType(multiPoint);

  if (multiPoint.isEmpty()) {
    _s << "EMPTY";
    return;
  }

  _s << "(";

  for (size_t i = 0; i < multiPoint.numGeometries(); i++) {
    if (i != 0) {
      _s << ",";
    }

    writeInner(multiPoint.geometryN(i).as<Point>());
  }

  _s << ")";
}

void
WktWriter::write(const MultiLineString &multiLineString)
{
  _s << "MULTILINESTRING ";
  writeCoordinateType(multiLineString);

  if (multiLineString.isEmpty()) {
    _s << "EMPTY";
    return;
  }

  _s << "(";

  for (size_t i = 0; i < multiLineString.numGeometries(); i++) {
    if (i != 0) {
      _s << ",";
    }

    writeInner(multiLineString.geometryN(i).as<LineString>());
  }

  _s << ")";
}

void
WktWriter::write(const MultiPolygon &multiPolygon)
{
  _s << "MULTIPOLYGON ";
  writeCoordinateType(multiPolygon);

  if (multiPolygon.isEmpty()) {
    _s << "EMPTY";
    return;
  }

  _s << "(";

  for (size_t i = 0; i < multiPolygon.numGeometries(); i++) {
    if (i != 0) {
      _s << ",";
    }

    writeInner(multiPolygon.geometryN(i).as<Polygon>());
  }

  _s << ")";
}

void
WktWriter::write(const MultiSolid &multiSolid)
{
  _s << "MULTISOLID ";
  writeCoordinateType(multiSolid);

  if (multiSolid.isEmpty()) {
    _s << "EMPTY";
    return;
  }

  _s << "(";

  for (size_t i = 0; i < multiSolid.numGeometries(); i++) {
    if (i != 0) {
      _s << ",";
    }

    writeInner(multiSolid.geometryN(i).as<Solid>());
  }

  _s << ")";
}

void
WktWriter::write(const Triangle &triangle)
{
  _s << "TRIANGLE ";
  writeCoordinateType(triangle);

  if (triangle.isEmpty()) {
    _s << "EMPTY";
    return;
  }

  writeInner(triangle);
}

void
WktWriter::writeInner(const Triangle &triangle)
{
  _s << "(";
  _s << "(";

  // close triangle
  for (size_t i = 0; i < 4; i++) {
    if (i != 0) {
      _s << ",";
    }

    writeCoordinate(triangle.vertex(static_cast<int>(i)));
  }

  _s << ")";
  _s << ")";
}

void
WktWriter::write(const TriangulatedSurface &tin)
{
  _s << "TIN ";
  writeCoordinateType(tin);

  if (tin.isEmpty()) {
    _s << "EMPTY";
    return;
  }

  _s << "("; // begin TIN

  for (size_t i = 0; i < tin.numPatches(); i++) {
    if (i != 0) {
      _s << ",";
    }

    writeInner(tin.patchN(i));
  }

  _s << ")"; // end TIN
}

void
WktWriter::write(const PolyhedralSurface &surface)
{
  _s << "POLYHEDRALSURFACE ";
  writeCoordinateType(surface);

  if (surface.isEmpty()) {
    _s << "EMPTY";
    return;
  }

  writeInner(surface);
}

void
WktWriter::writeInner(const PolyhedralSurface &surface)
{
  _s << "("; // begin POLYHEDRALSURFACE

  for (size_t i = 0; i < surface.numPatches(); i++) {
    if (i != 0) {
      _s << ",";
    }

    writeInner(surface.patchN(i));
  }

  _s << ")"; // end POLYHEDRALSURFACE
}

void
WktWriter::write(const Solid &solid)
{
  _s << "SOLID ";
  writeCoordinateType(solid);

  if (solid.isEmpty()) {
    _s << "EMPTY";
    return;
  }

  writeInner(solid);
}

void
WktWriter::writeInner(const Solid &solid)
{
  _s << "("; // begin SOLID
  writeInner(solid.exteriorShell());

  for (size_t i = 0; i < solid.numInteriorShells(); i++) {
    _s << ",";
    writeInner(solid.interiorShellN(i));
  }

  _s << ")"; // end SOLID
}

void
WktWriter::write(const NURBSCurve &nurbsCurve)
{
  _s << "NURBSCURVE ";
  writeCoordinateType(nurbsCurve);

  if (nurbsCurve.isEmpty()) {
    _s << "EMPTY";
    return;
  }

  writeInner(nurbsCurve);
}

void
WktWriter::writeInner(const NURBSCurve &nurbsCurve)
{
  _s << "(";

  // NEW FORMAT: Always write degree first (ISO compliant)
  _s << nurbsCurve.degree();

  // Write control points
  _s << ",(";
  for (size_t i = 0; i < nurbsCurve.numControlPoints(); i++) {
    if (i != 0) {
      _s << ",";
    }
    writeCoordinate(nurbsCurve.controlPointN(i));
  }
  _s << ")";

  // Always write weights if they exist to ensure round-trip consistency
  if (!nurbsCurve.weights().empty()) {
    _s << ",";
    writeWeights(nurbsCurve.weights());
  }

  // Always write knot vector if available to ensure round-trip consistency
  const auto &knots = nurbsCurve.knotVector();
  if (!knots.empty()) {
    _s << ",";
    writeKnots(knots);
  }

  _s << ")";
}

void
WktWriter::writeWeights(const std::vector<Kernel::FT> &weights)
{
  _s << "(";
  for (size_t i = 0; i < weights.size(); i++) {
    if (i != 0) {
      _s << ",";
    }

    // Always ensure weight is positive - critical for NURBS validity
    // Check if weight is zero or negative using CGAL comparison
    if (weights[i] <= Kernel::FT(0)) {
      // Throw exception for invalid weights instead of silent coercion
      std::ostringstream oss;
      oss << "Invalid NURBS weight at index " << i << ": "
          << CGAL::to_double(weights[i]) << ". Weights must be positive.";
      BOOST_THROW_EXCEPTION(Exception(oss.str()));
    } else {
      // Original weight is valid, use it
      if (_exactWrite) {
        impl::writeFT(_s, CGAL::exact(weights[i]));
      } else {
        double weight_val = CGAL::to_double(weights[i]);
        _s << fixZeroNegForWeights(weight_val,
                                   static_cast<int>(_s.precision()));
      }
    }
  }
  _s << ")";
}

void
WktWriter::writeKnots(const std::vector<Kernel::FT> &knots)
{
  _s << "(";
  for (size_t i = 0; i < knots.size(); i++) {
    if (i != 0) {
      _s << ",";
    }
    if (_exactWrite) {
      impl::writeFT(_s, CGAL::exact(knots[i]));
    } else {
      _s << fixZeroNeg(CGAL::to_double(knots[i]),
                       static_cast<int>(_s.precision()));
    }
  }
  _s << ")";
}

} // namespace SFCGAL::detail::io
