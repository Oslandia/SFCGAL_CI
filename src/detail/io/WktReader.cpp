// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/io/WktReader.h"

#include <memory>

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

WktReader::WktReader(std::istream &inputStream) : _reader(inputStream) {}

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
WktReader::readGeometry() -> std::unique_ptr<Geometry>
{
  // Check recursion depth to prevent stack overflow (CWE-674)
  if (_recursionDepth >= MAX_RECURSION_DEPTH) {
    BOOST_THROW_EXCEPTION(
        WktParseException("WktReader: maximum recursion depth exceeded"));
  }

  // RAII guard to ensure depth is decremented on exit
  struct RecursionGuard {
    int &depth;
    RecursionGuard(int &d) : depth(d) { ++depth; }
    ~RecursionGuard() { --depth; }
  } guard(_recursionDepth);

  GeometryType const geometryType = readGeometryType();
  _is3D                           = _reader.imatch("Z");
  _isMeasured                     = _reader.imatch("M");

  switch (geometryType) {
  case TYPE_POINT: {
    auto geom = std::make_unique<Point>();
    readInnerPoint(*geom);
    return geom;
  }

  case TYPE_LINESTRING: {
    auto geom = std::make_unique<LineString>();
    readInnerLineString(*geom);
    return geom;
  }

  case TYPE_TRIANGLE: {
    auto geom = std::make_unique<Triangle>();
    readInnerTriangle(*geom);
    return geom;
  }

  case TYPE_POLYGON: {
    auto geom = std::make_unique<Polygon>();
    readInnerPolygon(*geom);
    return geom;
  }

  case TYPE_MULTIPOINT: {
    auto geom = std::make_unique<MultiPoint>();
    readInnerMultiPoint(*geom);
    return geom;
  }

  case TYPE_MULTILINESTRING: {
    auto geom = std::make_unique<MultiLineString>();
    readInnerMultiLineString(*geom);
    return geom;
  }

  case TYPE_MULTIPOLYGON: {
    auto geom = std::make_unique<MultiPolygon>();
    readInnerMultiPolygon(*geom);
    return geom;
  }

  case TYPE_GEOMETRYCOLLECTION: {
    auto geom = std::make_unique<GeometryCollection>();
    readInnerGeometryCollection(*geom);
    return geom;
  }

  case TYPE_TRIANGULATEDSURFACE: {
    auto geom = std::make_unique<TriangulatedSurface>();
    readInnerTriangulatedSurface(*geom);
    return geom;
  }

  case TYPE_POLYHEDRALSURFACE: {
    auto geom = std::make_unique<PolyhedralSurface>();
    readInnerPolyhedralSurface(*geom);
    return geom;
  }

  case TYPE_SOLID: {
    auto geom = std::make_unique<Solid>();
    readInnerSolid(*geom);
    return geom;
  }

  case TYPE_MULTISOLID: {
    auto geom = std::make_unique<MultiSolid>();
    readInnerMultiSolid(*geom);
    return geom;
  }

  case TYPE_NURBSCURVE: {
    auto geom = std::make_unique<NURBSCurve>();
    readInnerNURBSCurve(*geom);
    return geom;
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
  if (_reader.imatch("NURBSCURVE")) {
    return TYPE_NURBSCURVE;
  }

  std::ostringstream errorStream;
  errorStream << "can't parse WKT geometry type (" << _reader.context() << ")";
  BOOST_THROW_EXCEPTION(WktParseException(errorStream.str()));
}

void
WktReader::readInnerPoint(Point &point)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  readPointCoordinate(point);

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }
}

void
WktReader::readInnerLineString(LineString &lineString)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  while (!_reader.eof()) {

    std::unique_ptr<Point> point(new Point());

    if (readPointCoordinate(*point)) {
      lineString.addPoint(point.release());
    } else {
      BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
    }

    if (_reader.match(',')) {
      continue;
    }

    break;
  }

  if (lineString.numPoints() < 2U) {
    BOOST_THROW_EXCEPTION(WktParseException(
        "WKT parse error, LineString should have at least 2 points"));
  }

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }
}

void
WktReader::readInnerPolygon(Polygon &polygon)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  for (int i = 0; !_reader.eof(); i++) {
    if (i == 0) {
      readInnerLineString(polygon.exteriorRing());
    } else {
      std::unique_ptr<LineString> interiorRing(new LineString);
      readInnerLineString(*interiorRing);
      polygon.addRing(interiorRing.release());
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
WktReader::readInnerTriangle(Triangle &triangle)
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

  triangle = Triangle(points[0], points[1], points[2]);

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }
}

void
WktReader::readInnerMultiPoint(MultiPoint &multiPoint)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  while (!_reader.eof()) {
    std::unique_ptr<Point> point(new Point());

    if (!_reader.imatch("EMPTY")) {
      // optional open/close parenthesis
      bool parenthesisOpen = false;

      if (_reader.match('(')) {
        parenthesisOpen = true;
      }

      readPointCoordinate(*point);

      if (parenthesisOpen && !_reader.match(')')) {
        BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
      }
    }

    if (!point->isEmpty()) {
      multiPoint.addGeometry(point.release());
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
WktReader::readInnerMultiLineString(MultiLineString &multiLineString)
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
      multiLineString.addGeometry(lineString.release());
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
WktReader::readInnerMultiPolygon(MultiPolygon &multiPolygon)
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
      multiPolygon.addGeometry(polygon.release());
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
WktReader::readInnerGeometryCollection(GeometryCollection &collection)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  while (!_reader.eof()) {
    // Save current state to avoid transactional state leaks
    bool saved_is3D       = _is3D;
    bool saved_isMeasured = _isMeasured;

    // read a full wkt geometry ex : POINT (2.0 6.0)
    auto gg = readGeometry();

    // Check for null before dereferencing to prevent null pointer access
    if (gg && !gg->isEmpty()) {
      collection.addGeometry(gg.release());
    }

    // Restore state for next iteration
    _is3D       = saved_is3D;
    _isMeasured = saved_isMeasured;

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
WktReader::readInnerTriangulatedSurface(TriangulatedSurface &tin)
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
    tin.addPatch(triangle.release());

    if (!_reader.match(',')) {
      break;
    }
  }

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }
}

void
WktReader::readInnerPolyhedralSurface(PolyhedralSurface &surface)
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
    surface.addPatch(polygon.release());

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
WktReader::readInnerSolid(Solid &solid)
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
      readInnerPolyhedralSurface(solid.exteriorShell());
    } else {
      std::unique_ptr<PolyhedralSurface> shell(new PolyhedralSurface);
      readInnerPolyhedralSurface(*shell);
      solid.addInteriorShell(shell.release());
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
WktReader::readInnerMultiSolid(MultiSolid &multiSolid)
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
      multiSolid.addGeometry(solid.release());
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

/// Helper functions for knot vector validation

static auto
validateKnotVector(const std::vector<Kernel::FT> &knots,
                   size_t numControlPoints, unsigned int degree) -> bool
{
  if (numControlPoints == 0) {
    return knots.empty();
  }

  size_t expectedSize = numControlPoints + degree + 1;
  if (knots.size() != expectedSize) {
    return false;
  }

  // Check non-decreasing order
  for (size_t knotIdx = 1; knotIdx < knots.size(); ++knotIdx) {
    if (knots[knotIdx] < knots[knotIdx - 1]) {
      return false;
    }
  }

  return true;
}

// NOLINTBEGIN(readability-function-cognitive-complexity)
void
WktReader::readInnerNURBSCurve(NURBSCurve &geometry)
{
  if (_reader.imatch("EMPTY")) {
    return;
  }

  // NEW FORMAT: Read degree first (ISO compliant)
  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  unsigned int degree = readDegree();

  if (!_reader.match(',')) {
    BOOST_THROW_EXCEPTION(WktParseException("Expected comma after degree"));
  }

  // Read control points - always present after degree
  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException("Expected control points list"));
  }

  std::vector<Point> controlPoints;
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

  if (controlPoints.empty()) {
    BOOST_THROW_EXCEPTION(
        WktParseException("NURBS curve must have at least one control point"));
  }

  // Validate degree
  if (degree >= controlPoints.size()) {
    BOOST_THROW_EXCEPTION(WktParseException(
        "NURBS degree must be less than number of control points"));
  }

  // Parse optional elements after control points: weights and/or knots
  std::vector<Kernel::FT> weights;
  std::vector<Kernel::FT> knots;

  if (_reader.match(',')) {
    // Look ahead to determine what comes next
    _reader.begin(); // Save position

    if (_reader.match('(')) {
      _reader.rollback(); // Restore position

      // Read first vector (could be weights or knots)
      std::vector<Kernel::FT> firstVector = readWeightsVector();

      if (_reader.match(',')) {
        _reader.begin(); // Save position

        if (_reader.match('(')) {
          _reader.rollback(); // Restore position
          // First was weights, second is knots
          weights = firstVector;
          knots   = readKnotsVector();
        } else {
          _reader.rollback(); // Restore position
          // Only weights provided (knot can be a scalar)
          weights = firstVector;
        }
      } else {
        // Only one vector provided - could be weights or knots
        // Assume weights if size matches control points, otherwise knots
        if (firstVector.size() == controlPoints.size()) {
          weights = firstVector;
        } else {
          knots = firstVector;
        }
      }
    }
  }

  if (!_reader.match(')')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  // Validate weights if provided
  if (!weights.empty() && weights.size() != controlPoints.size()) {
    BOOST_THROW_EXCEPTION(WktParseException(
        "Number of weights must match number of control points"));
  }

  // Validate that all weights are positive
  if (!weights.empty()) {
    for (size_t i = 0; i < weights.size(); ++i) {
      if (weights[i] <= 0) {
        BOOST_THROW_EXCEPTION(WktParseException(
            "Weight " + std::to_string(i) +
            " is non-positive: " + std::to_string(CGAL::to_double(weights[i])) +
            ". All weights must be positive"));
      }
    }
  }

  // Validate knot vector if provided
  if (!knots.empty() &&
      !validateKnotVector(knots, controlPoints.size(), degree)) {
    BOOST_THROW_EXCEPTION(WktParseException(
        "Invalid knot vector: must be non-decreasing and have exactly " +
        std::to_string(controlPoints.size() + degree + 1) + " elements"));
  }

  // Construct the NURBS curve using appropriate constructor
  if (!knots.empty()) {
    // Full constructor with knots
    if (weights.empty()) {
      weights = std::vector<Kernel::FT>(controlPoints.size(), Kernel::FT(1));
    }
    geometry = NURBSCurve(controlPoints, weights, degree, knots);
  } else if (!weights.empty()) {
    // Constructor with weights but no knots
    geometry = NURBSCurve(controlPoints, weights, degree);
  } else {
    // Simple constructor with just control points and degree
    geometry = NURBSCurve(controlPoints, degree);
  }
}
// NOLINTEND(readability-function-cognitive-complexity)

auto
WktReader::readWeightsVector() -> std::vector<Kernel::FT>
{
  std::vector<Kernel::FT> weights;

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  while (!_reader.eof()) {
    Kernel::FT weight;
    if (!_reader.read(weight)) {
      BOOST_THROW_EXCEPTION(WktParseException("Expected weight value"));
    }

    if (weight <= Kernel::FT(0)) {
      // Instead of throwing an error, replace invalid weights with 1.0
      // This provides better robustness for WKB round-trip compatibility
      weight = Kernel::FT(1);
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

auto
WktReader::readKnotsVector() -> std::vector<Kernel::FT>
{
  std::vector<Kernel::FT> knots;

  if (!_reader.match('(')) {
    BOOST_THROW_EXCEPTION(WktParseException(parseErrorMessage()));
  }

  while (!_reader.eof()) {
    Kernel::FT knot;
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

auto
WktReader::readDegree() -> unsigned int
{
  unsigned int degree;
  if (!_reader.read(degree)) {
    BOOST_THROW_EXCEPTION(WktParseException("Expected degree value"));
  }
  return degree;
}

auto
WktReader::readPointCoordinate(Point &point) -> bool
{
  std::vector<Kernel::FT> coordinates;
  Kernel::FT              coordinate;

  if (_reader.imatch("EMPTY")) {
    point = Point();
    return false;
  }

  while (_reader.read(coordinate)) {
    coordinates.push_back(coordinate);
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

    point = Point(coordinates[0], coordinates[1], coordinates[2]);
    point.setM(CGAL::to_double(coordinates[3]));
  } else if (_isMeasured && !_is3D) {
    // XYM
    if (coordinates.size() != 3) {
      BOOST_THROW_EXCEPTION(WktParseException(
          "bad coordinate dimension (expecting XYM coordinates)"));
    }

    point = Point(coordinates[0], coordinates[1]);
    point.setM(CGAL::to_double(coordinates[2]));
  } else if (coordinates.size() == 3) {
    // XYZ
    point = Point(coordinates[0], coordinates[1], coordinates[2]);
  } else {
    // XY
    point = Point(coordinates[0], coordinates[1]);
  }

  return true;
}

auto
WktReader::parseErrorMessage() -> std::string
{
  std::ostringstream errorStream;
  errorStream << "WKT parse error (" << _reader.context() << ")";
  return errorStream.str();
}

} // namespace SFCGAL::detail::io
