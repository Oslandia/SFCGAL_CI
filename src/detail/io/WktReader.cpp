// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/io/WktReader.h"

#include <memory>

#include "SFCGAL/BSplineCurve.h"
#include "SFCGAL/BezierCurve.h"
#include "SFCGAL/GeometryCollection.h"
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

#include "SFCGAL/Exception.h"

namespace SFCGAL::detail::io {

WktReader::WktReader(std::istream &s) : _reader(s) {}

auto
WktReader::readSRID() -> srid_t
{
  srid_t srid = 0;

  if (_reader.imatch("SRID=")) {
    _reader.read(srid);

    if (!_reader.match(";")) {
      BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
    }
  }

  return srid;
}

auto
WktReader::readGeometry() -> Geometry *
{
  GeometryType const geometryType = readGeometryType();
  _is3D                           = _reader.imatch("Z");
  _isMeasured                     = _reader.imatch("M");

  switch (geometryType) {
  case TYPE_POINT: {
    std::unique_ptr<Point> g(new Point());
    readInnerPoint(*g);
    return g.release();
  }

  case TYPE_LINESTRING: {
    std::unique_ptr<LineString> g(new LineString());
    readInnerLineString(*g);
    return g.release();
  }

  case TYPE_TRIANGLE: {
    std::unique_ptr<Triangle> g(new Triangle());
    readInnerTriangle(*g);
    return g.release();
  }

  case TYPE_POLYGON: {
    std::unique_ptr<Polygon> g(new Polygon());
    readInnerPolygon(*g);
    return g.release();
  }

  case TYPE_MULTIPOINT: {
    std::unique_ptr<MultiPoint> g(new MultiPoint());
    readInnerMultiPoint(*g);
    return g.release();
  }

  case TYPE_MULTILINESTRING: {
    std::unique_ptr<MultiLineString> g(new MultiLineString());
    readInnerMultiLineString(*g);
    return g.release();
  }

  case TYPE_MULTIPOLYGON: {
    std::unique_ptr<MultiPolygon> g(new MultiPolygon());
    readInnerMultiPolygon(*g);
    return g.release();
  }

  case TYPE_GEOMETRYCOLLECTION: {
    std::unique_ptr<GeometryCollection> g(new GeometryCollection());
    readInnerGeometryCollection(*g);
    return g.release();
  }

  case TYPE_TRIANGULATEDSURFACE: {
    std::unique_ptr<TriangulatedSurface> g(new TriangulatedSurface());
    readInnerTriangulatedSurface(*g);
    return g.release();
  }

  case TYPE_POLYHEDRALSURFACE: {
    std::unique_ptr<PolyhedralSurface> g(new PolyhedralSurface());
    readInnerPolyhedralSurface(*g);
    return g.release();
  }

  case TYPE_SOLID: {
    std::unique_ptr<Solid> g(new Solid());
    readInnerSolid(*g);
    return g.release();
  }

  case TYPE_MULTISOLID: {
    std::unique_ptr<MultiSolid> g(new MultiSolid());
    readInnerMultiSolid(*g);
    return g.release();
  }

  case TYPE_BEZIERCURVE: {
    std::unique_ptr<BezierCurve> g(new BezierCurve());
    readInnerBezierCurve(*g);
    return g.release();
  }

  case TYPE_BSPLINECURVE: {
    std::unique_ptr<BSplineCurve> g(new BSplineCurve());
    readInnerBSplineCurve(*g);
    return g.release();
  }

  case TYPE_NURBSCURVE: {
    std::unique_ptr<NURBSCurve> g(new NURBSCurve());
    readInnerNURBSCurve(*g);
    return g.release();
  }
  }

  BOOST_THROW_EXCEPTION(WktParseException("unexpected geometry"));
}

auto
WktReader::readGeometryType() -> GeometryType
{
  if (_reader.imatch("POINT")) {
    return TYPE_POINT;
  }
  if (_reader.imatch("LINESTRING")) {
    return TYPE_LINESTRING;
  }
  if (_reader.imatch("POLYGON")) {
    return TYPE_POLYGON;
  }
  if (_reader.imatch("TRIANGLE")) {
    // not official
    return TYPE_TRIANGLE;
  }
  if (_reader.imatch("MULTIPOINT")) {
    return TYPE_MULTIPOINT;
  }
  if (_reader.imatch("MULTILINESTRING")) {
    return TYPE_MULTILINESTRING;
  }
  if (_reader.imatch("MULTIPOLYGON")) {
    return TYPE_MULTIPOLYGON;
  }
  if (_reader.imatch("GEOMETRYCOLLECTION")) {
    return TYPE_GEOMETRYCOLLECTION;
  }
  if (_reader.imatch("TIN")) {
    return TYPE_TRIANGULATEDSURFACE;
  }
  if (_reader.imatch("POLYHEDRALSURFACE")) {
    return TYPE_POLYHEDRALSURFACE;
  }
  if (_reader.imatch("SOLID")) {
    // not official
    return TYPE_SOLID;
  }
  if (_reader.imatch("MULTISOLID")) {
    // not official
    return TYPE_MULTISOLID;
  }
  if (_reader.imatch("BEZIERCURVE")) {
    // not official
    return TYPE_BEZIERCURVE;
  }
  if (_reader.imatch("BSPLINECURVE")) {
    // not official
    return TYPE_BSPLINECURVE;
  }

  if (_reader.imatch("NURBSCURVE")) {
    return TYPE_NURBSCURVE;
  }

  std::ostringstream oss;
  oss << "can't parse WKT geometry type (" << _reader.context() << ")";
  BOOST_THROW_EXCEPTION(WktParseException(oss.str()));
}

void
WktReader::readInnerPoint(Point &g)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  readPointCoordinate(g);

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }
}

void
WktReader::readInnerLineString(LineString &g)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  while (!_reader.eof()) {

    std::unique_ptr<Point> p(new Point());

    if (readPointCoordinate(*p)) {
      g.addPoint(p.release());
    } else {
      BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
    }

    if (_reader.match(',')) {
      continue;
    }

    break;
  }

  if (g.numPoints() < 2U) {
    BOOST_THROW_EXCEPTION(WktParseException(
        "WKT parse error, LineString should have at least 2 points"));
  }

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }
}

void
WktReader::readInnerPolygon(Polygon &g)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  for (int i = 0; !_reader.eof(); i++) {
    if (i == 0) {
      readInnerLineString(g.exteriorRing());
    } else {
      std::unique_ptr<LineString> interiorRing(new LineString);
      readInnerLineString(*interiorRing);
      g.addRing(interiorRing.release());
    }

    // break if not followed by another ring
    if (!_reader.match(',')) {
      break;
    }
  }

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }
}

void
WktReader::readInnerTriangle(Triangle &g)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  // 4 points to read
  std::vector<Point> points;

  while (!_reader.eof()) {
    points.emplace_back();
    readPointCoordinate(points.back());

    if (!_reader.match(",")) {
      break;
    }
  }

  if (points.size() != 4) {
    BOOST_THROW_EXCEPTION(
        WktParseException("WKT parse error, expected 4 points for triangle"));
  }

  if (points.back() != points.front()) {
    BOOST_THROW_EXCEPTION(
        WktParseException("WKT parse error, first point different of the last "
                          "point for triangle"));
  }

  g = Triangle(points[0], points[1], points[2]);

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }
}

void
WktReader::readInnerMultiPoint(MultiPoint &g)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  while (!_reader.eof()) {
    std::unique_ptr<Point> p(new Point());

    if (!_reader.imatch("EMPTY")) {
      // optional open/close parenthesis
      bool parenthesisOpen = false;

      if (_reader.match('(')) {
        parenthesisOpen = true;
      }

      readPointCoordinate(*p);

      if (parenthesisOpen && !_reader.match(')')) {
        BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
      }
    }

    if (!p->isEmpty()) {
      g.addGeometry(p.release());
    }

    // break if not followed by another points
    if (!_reader.match(',')) {
      break;
    }
  }

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }
}

void
WktReader::readInnerMultiLineString(MultiLineString &g)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  while (!_reader.eof()) {

    std::unique_ptr<LineString> lineString(new LineString());
    readInnerLineString(*lineString);
    if (!lineString->isEmpty()) {
      g.addGeometry(lineString.release());
    }

    // break if not followed by another points
    if (!_reader.match(',')) {
      break;
    }
  }

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }
}

void
WktReader::readInnerMultiPolygon(MultiPolygon &g)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  while (!_reader.eof()) {

    std::unique_ptr<Polygon> polygon(new Polygon());
    readInnerPolygon(*polygon);
    if (!polygon->isEmpty()) {
      g.addGeometry(polygon.release());
    }

    // break if not followed by another points
    if (!_reader.match(',')) {
      break;
    }
  }

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }
}

void
WktReader::readInnerGeometryCollection(GeometryCollection &g)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  while (!_reader.eof()) {

    // read a full wkt geometry ex : POINT (2.0 6.0)
    Geometry *gg = readGeometry();
    if (!gg->isEmpty()) {
      g.addGeometry(gg);
    }

    // break if not followed by another points
    if (!_reader.match(',')) {
      break;
    }
  }

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }
}

void
WktReader::readInnerTriangulatedSurface(TriangulatedSurface &g)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  while (!_reader.eof()) {
    std::unique_ptr<Triangle> triangle(new Triangle());
    readInnerTriangle(*triangle);
    g.addPatch(triangle.release());

    if (!_reader.match(',')) {
      break;
    }
  }

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }
}

void
WktReader::readInnerPolyhedralSurface(PolyhedralSurface &g)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  while (!_reader.eof()) {
    std::unique_ptr<Polygon> polygon(new Polygon());
    readInnerPolygon(*polygon);
    g.addPatch(polygon.release());

    // break if not followed by another points
    if (!_reader.match(',')) {
      break;
    }
  }

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }
}

void
WktReader::readInnerSolid(Solid &g)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  // solid begin
  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  for (int i = 0; !_reader.eof(); i++) {
    if (i == 0) {
      readInnerPolyhedralSurface(g.exteriorShell());
    } else {
      std::unique_ptr<PolyhedralSurface> shell(new PolyhedralSurface);
      readInnerPolyhedralSurface(*shell);
      g.addInteriorShell(shell.release());
    }

    // break if not followed by another points
    if (!_reader.match(',')) {
      break;
    }
  }

  // solid end
  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }
}

void
WktReader::readInnerMultiSolid(MultiSolid &g)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  while (!_reader.eof()) {

    std::unique_ptr<Solid> solid(new Solid());
    readInnerSolid(*solid);
    if (!solid->isEmpty()) {
      g.addGeometry(solid.release());
    }

    // break if not followed by another points
    if (!_reader.match(',')) {
      break;
    }
  }

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }
}

void
WktReader::readInnerBezierCurve(BezierCurve &g)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  while (!_reader.eof()) {
    Point point;
    readPointCoordinate(point);
    g.addControlPoint(point);

    // Break if not followed by comma
    if (!_reader.match(',')) {
      break;
    }
  }

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }
}

void
WktReader::readInnerBSplineCurve(BSplineCurve &g)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  while (!_reader.eof()) {
    Point point;
    readPointCoordinate(point);
    g.addControlPoint(point);

    if (!_reader.match(',')) {
      break;
    }
  }

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  if (!_reader.match(',')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  int degree;
  if (!_reader.read(degree)) {
    BOOST_THROW_EXCEPTION(
        WktParseException("Expected degree value for BSplineCurve"));
  }

  if (degree < 0) {
    BOOST_THROW_EXCEPTION(WktParseException("Degree cannot be negative"));
  }

  g.setDegree(static_cast<unsigned int>(degree));

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }
}

void
WktReader::readInnerNURBSCurve(NURBSCurve &g)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  // Read control points
  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  std::vector<Point> controlPoints;
  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  while (!_reader.eof()) {
    Point point;
    readPointCoordinate(point);
    controlPoints.push_back(point);

    if (!_reader.match(',')) {
      break;
    }
  }

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  // Parse optional elements after control points
  std::vector<double> weights;
  std::vector<double> knots;
  unsigned int        degree =
      controlPoints.empty()
                 ? 0
                 : static_cast<unsigned int>(controlPoints.size() - 1);

  if (_reader.match(',')) {
    _reader.begin(); // Save position

    if (_reader.match('(')) {
      _reader.rollback(); // Restore position

      // Read first vector (could be weights or knots)
      std::vector<double> firstVector = readWeightsVector();

      if (_reader.match(',')) {
        _reader.begin(); // Save position

        if (_reader.match('(')) {
          _reader.rollback(); // Restore position
          // First was weights, second is knots
          weights = firstVector;
          knots   = readKnotsVector();

          if (_reader.match(',')) {
            degree = readDegree();
          }
        } else {
          _reader.rollback(); // Restore position
          // First was weights, next is degree
          weights = firstVector;
          degree  = readDegree();
        }
      } else {
        // Only one vector provided - assume it's weights
        weights = firstVector;
      }
    } else {
      _reader.rollback(); // Restore position
      // Direct degree
      degree = readDegree();
    }
  }

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  // Validate knot vector if provided
  if (!knots.empty() &&
      !validateKnotVector(knots, controlPoints.size(), degree)) {
    BOOST_THROW_EXCEPTION(
        WktParseException("Invalid knot vector for NURBS curve"));
  }

  // Construct the NURBS curve
  if (!knots.empty()) {
    if (weights.empty()) {
      weights = std::vector<double>(controlPoints.size(), 1.0);
    }
    g = NURBSCurve(controlPoints, degree, weights, knots);
  } else if (!weights.empty()) {
    g = NURBSCurve(controlPoints, degree, weights);
  } else {
    g = NURBSCurve(controlPoints, degree);
  }
}

std::vector<double>
WktReader::readWeightsVector()
{
  std::vector<double> weights;

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  while (!_reader.eof()) {
    double weight;
    if (!_reader.read(weight)) {
      BOOST_THROW_EXCEPTION(WktParseException("Expected weight value"));
    }

    if (weight <= 0.0) {
      BOOST_THROW_EXCEPTION(
          WktParseException("NURBS weights must be positive"));
    }

    weights.push_back(weight);

    if (!_reader.match(',')) {
      break;
    }
  }

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  return weights;
}

std::vector<double>
WktReader::readKnotsVector()
{
  std::vector<double> knots;

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  while (!_reader.eof()) {
    double knot;
    if (!_reader.read(knot)) {
      BOOST_THROW_EXCEPTION(WktParseException("Expected knot value"));
    }

    knots.push_back(knot);

    if (!_reader.match(',')) {
      break;
    }
  }

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  return knots;
}

unsigned int
WktReader::readDegree()
{
  unsigned int degree;
  if (!_reader.read(degree)) {
    BOOST_THROW_EXCEPTION(WktParseException("Expected degree value"));
  }
  return degree;
}

auto
WktReader::readPointCoordinate(Point &p) -> bool
{
  std::vector<Kernel::Exact_kernel::FT> coordinates;
  Kernel::Exact_kernel::FT              d;

  if (_reader.imatch("EMPTY")) {
    p = Point();
    return false;
  }

  while (_reader.read(d)) {
    coordinates.push_back(d);
  }

  if (coordinates.size() < 2) {
    BOOST_THROW_EXCEPTION(WktParseException(
        (boost::format("WKT parse error, Coordinate dimension < 2 (%s)") %
         _reader.context())
            .str()));
  }

  if (coordinates.size() > 4) {
    BOOST_THROW_EXCEPTION(
        WktParseException("WKT parse error, Coordinate dimension > 4"));
  }

  if (_isMeasured && _is3D) {
    // XYZM
    if (coordinates.size() != 4) {
      BOOST_THROW_EXCEPTION(WktParseException("bad coordinate dimension"));
    }

    p = Point(coordinates[0], coordinates[1], coordinates[2]);
    p.setM(CGAL::to_double(coordinates[3]));
  } else if (_isMeasured && !_is3D) {
    // XYM
    if (coordinates.size() != 3) {
      BOOST_THROW_EXCEPTION(WktParseException(
          "bad coordinate dimension (expecting XYM coordinates)"));
    }

    p = Point(coordinates[0], coordinates[1]);
    p.setM(CGAL::to_double(coordinates[2]));
  } else if (coordinates.size() == 3) {
    // XYZ
    p = Point(coordinates[0], coordinates[1], coordinates[2]);
  } else {
    // XY
    p = Point(coordinates[0], coordinates[1]);
  }

  return true;
}

auto
WktReader::parseErrorMessage() -> std::string
{
  std::ostringstream oss;
  oss << "WKT parse error (" << _reader.context() << ")";
  return oss.str();
}

/// Helper functions for knot vector

std::vector<double>
WktReader::generateUniformKnots(size_t numControlPoints, unsigned int degree)
{
  if (numControlPoints == 0) {
    return {};
  }

  size_t              knotCount = numControlPoints + degree + 1;
  std::vector<double> knots(knotCount);

  // Clamp knots at start
  for (size_t i = 0; i <= degree; ++i) {
    knots[i] = 0.0;
  }

  // Internal knots uniformly spaced
  size_t segments = numControlPoints - degree;
  for (size_t i = degree + 1; i < knotCount - degree - 1; ++i) {
    knots[i] = static_cast<double>(i - degree) / static_cast<double>(segments);
  }

  // Clamp knots at end
  for (size_t i = knotCount - degree - 1; i < knotCount; ++i) {
    knots[i] = 1.0;
  }

  return knots;
}

bool
WktReader::validateKnotVector(const std::vector<double> &knots,
                              size_t numControlPoints, unsigned int degree)
{
  if (numControlPoints == 0) {
    return knots.empty();
  }

  size_t expectedSize = numControlPoints + degree + 1;
  if (knots.size() != expectedSize) {
    return false;
  }

  // Check non-decreasing order
  for (size_t i = 1; i < knots.size(); ++i) {
    if (knots[i] < knots[i - 1]) {
      return false;
    }
  }

  return true;
}

} // namespace SFCGAL::detail::io
