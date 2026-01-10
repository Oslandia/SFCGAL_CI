// Copyright (c) 2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <memory>
#include <string>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/PreparedGeometry.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/io/geojson.h"

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

using namespace SFCGAL;
using namespace SFCGAL::io;

BOOST_AUTO_TEST_SUITE(SFCGAL_io_GeoJSONTest)

// ============================================================================
// READING TESTS
// ============================================================================

BOOST_AUTO_TEST_CASE(testReadPoint2D)
{
  std::string const json = R"({"type": "Point", "coordinates": [1.0, 2.0]})";
  std::unique_ptr<Geometry> g = readGeoJSON(json);

  BOOST_CHECK(g != nullptr);
  BOOST_CHECK(g->is<Point>());
  BOOST_CHECK(!g->isEmpty());

  const auto &pt = g->as<Point>();
  BOOST_CHECK_CLOSE(CGAL::to_double(pt.x()), 1.0, 0.0001);
  BOOST_CHECK_CLOSE(CGAL::to_double(pt.y()), 2.0, 0.0001);
  BOOST_CHECK(!pt.is3D());
}

BOOST_AUTO_TEST_CASE(testReadPoint3D)
{
  std::string const json =
      R"({"type": "Point", "coordinates": [1.0, 2.0, 3.0]})";
  std::unique_ptr<Geometry> g = readGeoJSON(json);

  BOOST_CHECK(g != nullptr);
  BOOST_CHECK(g->is<Point>());

  const auto &pt = g->as<Point>();
  BOOST_CHECK(pt.is3D());
  BOOST_CHECK_CLOSE(CGAL::to_double(pt.z()), 3.0, 0.0001);
}

BOOST_AUTO_TEST_CASE(testReadLineString)
{
  std::string const json =
      R"({"type": "LineString", "coordinates": [[0,0], [1,1], [2,0]]})";
  std::unique_ptr<Geometry> g = readGeoJSON(json);

  BOOST_CHECK(g->is<LineString>());
  const auto &lineString = g->as<LineString>();
  BOOST_CHECK_EQUAL(lineString.numPoints(), 3U);
}

BOOST_AUTO_TEST_CASE(testReadLineString3D)
{
  std::string const json =
      R"({"type": "LineString", "coordinates": [[0,0,0], [1,1,1], [2,0,2]]})";
  std::unique_ptr<Geometry> g = readGeoJSON(json);

  BOOST_CHECK(g->is<LineString>());
  const auto &lineString = g->as<LineString>();
  BOOST_CHECK_EQUAL(lineString.numPoints(), 3U);
  BOOST_CHECK(lineString.is3D());
}

BOOST_AUTO_TEST_CASE(testReadPolygon)
{
  std::string const         json = R"({
        "type": "Polygon",
        "coordinates": [
            [[0,0], [10,0], [10,10], [0,10], [0,0]]
        ]
    })";
  std::unique_ptr<Geometry> g    = readGeoJSON(json);

  BOOST_CHECK(g->is<Polygon>());
  const auto &poly = g->as<Polygon>();
  BOOST_CHECK_EQUAL(poly.exteriorRing().numPoints(), 5U);
  BOOST_CHECK_EQUAL(poly.numInteriorRings(), 0U);
}

BOOST_AUTO_TEST_CASE(testReadPolygonWithHole)
{
  std::string const         json = R"({
        "type": "Polygon",
        "coordinates": [
            [[0,0], [10,0], [10,10], [0,10], [0,0]],
            [[2,2], [8,2], [8,8], [2,8], [2,2]]
        ]
    })";
  std::unique_ptr<Geometry> g    = readGeoJSON(json);

  BOOST_CHECK(g->is<Polygon>());
  const auto &poly = g->as<Polygon>();
  BOOST_CHECK_EQUAL(poly.numInteriorRings(), 1U);
}

BOOST_AUTO_TEST_CASE(testReadMultiPoint)
{
  std::string const json =
      R"({"type": "MultiPoint", "coordinates": [[0,0], [1,1], [2,2]]})";
  std::unique_ptr<Geometry> g = readGeoJSON(json);

  BOOST_CHECK(g->is<MultiPoint>());
  const auto &mp = g->as<MultiPoint>();
  BOOST_CHECK_EQUAL(mp.numGeometries(), 3U);
}

BOOST_AUTO_TEST_CASE(testReadMultiLineString)
{
  std::string const         json = R"({
        "type": "MultiLineString",
        "coordinates": [
            [[0,0], [1,1]],
            [[2,2], [3,3]]
        ]
    })";
  std::unique_ptr<Geometry> g    = readGeoJSON(json);

  BOOST_CHECK(g->is<MultiLineString>());
  const auto &mls = g->as<MultiLineString>();
  BOOST_CHECK_EQUAL(mls.numGeometries(), 2U);
}

BOOST_AUTO_TEST_CASE(testReadMultiPolygon)
{
  std::string const         json = R"({
        "type": "MultiPolygon",
        "coordinates": [
            [[[0,0], [1,0], [1,1], [0,1], [0,0]]],
            [[[2,2], [3,2], [3,3], [2,3], [2,2]]]
        ]
    })";
  std::unique_ptr<Geometry> g    = readGeoJSON(json);

  BOOST_CHECK(g->is<MultiPolygon>());
  const auto &mp = g->as<MultiPolygon>();
  BOOST_CHECK_EQUAL(mp.numGeometries(), 2U);
}

BOOST_AUTO_TEST_CASE(testReadGeometryCollection)
{
  std::string const         json = R"({
        "type": "GeometryCollection",
        "geometries": [
            {"type": "Point", "coordinates": [0, 0]},
            {"type": "LineString", "coordinates": [[0,0], [1,1]]}
        ]
    })";
  std::unique_ptr<Geometry> g    = readGeoJSON(json);

  BOOST_CHECK(g->is<GeometryCollection>());
  const auto &gc = g->as<GeometryCollection>();
  BOOST_CHECK_EQUAL(gc.numGeometries(), 2U);
}

BOOST_AUTO_TEST_CASE(testReadFeature)
{
  std::string const         json = R"({
        "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [1, 2]},
        "properties": {"name": "test"}
    })";
  std::unique_ptr<Geometry> g    = readGeoJSON(json);

  BOOST_CHECK(g->is<Point>());
}

BOOST_AUTO_TEST_CASE(testReadFeatureCollection)
{
  std::string const         json = R"({
        "type": "FeatureCollection",
        "features": [
            {"type": "Feature", "geometry": {"type": "Point", "coordinates": [0, 0]}, "properties": {}},
            {"type": "Feature", "geometry": {"type": "Point", "coordinates": [1, 1]}, "properties": {}}
        ]
    })";
  std::unique_ptr<Geometry> g    = readGeoJSON(json);

  BOOST_CHECK(g->is<GeometryCollection>());
  const auto &gc = g->as<GeometryCollection>();
  BOOST_CHECK_EQUAL(gc.numGeometries(), 2U);
}

BOOST_AUTO_TEST_CASE(testReadTriangleExtension)
{
  std::string const         json = R"({
        "type": "Triangle",
        "coordinates": [[[0,0,0], [1,0,0], [0.5,1,0], [0,0,0]]]
    })";
  std::unique_ptr<Geometry> g    = readGeoJSON(json);

  BOOST_CHECK(g->is<Triangle>());
  const auto &tri = g->as<Triangle>();
  BOOST_CHECK(tri.is3D());
}

BOOST_AUTO_TEST_CASE(testReadTINExtension)
{
  std::string const         json = R"({
        "type": "TIN",
        "coordinates": [
            [[[0,0,0], [1,0,0], [0.5,1,0], [0,0,0]]],
            [[[1,0,0], [2,0,0], [1.5,1,0], [1,0,0]]]
        ]
    })";
  std::unique_ptr<Geometry> g    = readGeoJSON(json);

  BOOST_CHECK(g->is<TriangulatedSurface>());
  const auto &tin = g->as<TriangulatedSurface>();
  BOOST_CHECK_EQUAL(tin.numTriangles(), 2U);
}

BOOST_AUTO_TEST_CASE(testReadPolyhedralSurfaceExtension)
{
  std::string const         json = R"({
        "type": "PolyhedralSurface",
        "coordinates": [
            [[[0,0,0], [1,0,0], [1,1,0], [0,1,0], [0,0,0]]],
            [[[0,0,1], [1,0,1], [1,1,1], [0,1,1], [0,0,1]]]
        ]
    })";
  std::unique_ptr<Geometry> g    = readGeoJSON(json);

  BOOST_CHECK(g->is<PolyhedralSurface>());
  const auto &polyhedralSurface = g->as<PolyhedralSurface>();
  BOOST_CHECK_EQUAL(polyhedralSurface.numPolygons(), 2U);
}

BOOST_AUTO_TEST_CASE(testReadPreparedWithCRS)
{
  std::string const                 json     = R"({
        "type": "Point",
        "coordinates": [1, 2],
        "crs": {
            "type": "name",
            "properties": {"name": "urn:ogc:def:crs:EPSG::4326"}
        }
    })";
  std::unique_ptr<PreparedGeometry> prepared = readGeoJSONPrepared(json);

  BOOST_CHECK_EQUAL(prepared->SRID(), 4326U);
}

BOOST_AUTO_TEST_CASE(testReadPreparedWithSimpleCRS)
{
  std::string const                 json     = R"({
        "type": "Point",
        "coordinates": [1, 2],
        "crs": {
            "type": "name",
            "properties": {"name": "EPSG:2154"}
        }
    })";
  std::unique_ptr<PreparedGeometry> prepared = readGeoJSONPrepared(json);

  BOOST_CHECK_EQUAL(prepared->SRID(), 2154U);
}

// ============================================================================
// WRITING TESTS
// ============================================================================

BOOST_AUTO_TEST_CASE(testWritePoint2D)
{
  Point const pt(1.5, 2.5);
  std::string json = writeGeoJSON(pt);

  auto parsed = nlohmann::json::parse(json);
  BOOST_CHECK_EQUAL(parsed["type"], "Point");
  BOOST_CHECK_CLOSE(parsed["coordinates"][0].get<double>(), 1.5, 0.0001);
  BOOST_CHECK_CLOSE(parsed["coordinates"][1].get<double>(), 2.5, 0.0001);
  BOOST_CHECK_EQUAL(parsed["coordinates"].size(), 2U);
}

BOOST_AUTO_TEST_CASE(testWritePoint3D)
{
  Point const pt(1.0, 2.0, 3.0);
  std::string json = writeGeoJSON(pt);

  auto parsed = nlohmann::json::parse(json);
  BOOST_CHECK_EQUAL(parsed["coordinates"].size(), 3U);
  BOOST_CHECK_CLOSE(parsed["coordinates"][2].get<double>(), 3.0, 0.0001);
}

BOOST_AUTO_TEST_CASE(testWriteLineString)
{
  LineString lineString;
  lineString.addPoint(Point(0, 0));
  lineString.addPoint(Point(1, 1));
  lineString.addPoint(Point(2, 0));

  std::string json   = writeGeoJSON(lineString);
  auto        parsed = nlohmann::json::parse(json);

  BOOST_CHECK_EQUAL(parsed["type"], "LineString");
  BOOST_CHECK_EQUAL(parsed["coordinates"].size(), 3U);
}

BOOST_AUTO_TEST_CASE(testWritePolygon)
{
  LineString ring;
  ring.addPoint(Point(0, 0));
  ring.addPoint(Point(1, 0));
  ring.addPoint(Point(1, 1));
  ring.addPoint(Point(0, 1));
  ring.addPoint(Point(0, 0));
  Polygon const poly(ring);

  std::string json   = writeGeoJSON(poly);
  auto        parsed = nlohmann::json::parse(json);

  BOOST_CHECK_EQUAL(parsed["type"], "Polygon");
  BOOST_CHECK_EQUAL(parsed["coordinates"].size(), 1U);    // 1 ring
  BOOST_CHECK_EQUAL(parsed["coordinates"][0].size(), 5U); // 5 points
}

BOOST_AUTO_TEST_CASE(testWriteTriangleStrict)
{
  Triangle const tri(Point(0, 0, 0), Point(1, 0, 0), Point(0.5, 1, 0));

  GeoJSONOptions opts;
  opts.strict      = true;
  std::string json = writeGeoJSON(tri, opts);

  auto parsed = nlohmann::json::parse(json);
  BOOST_CHECK_EQUAL(parsed["type"], "Polygon");
  BOOST_CHECK_EQUAL(parsed["properties"]["__sfcgal_original_type"], "Triangle");
}

BOOST_AUTO_TEST_CASE(testWriteTriangleExtension)
{
  Triangle const tri(Point(0, 0, 0), Point(1, 0, 0), Point(0.5, 1, 0));

  GeoJSONOptions opts;
  opts.strict      = false;
  std::string json = writeGeoJSON(tri, opts);

  auto parsed = nlohmann::json::parse(json);
  BOOST_CHECK_EQUAL(parsed["type"], "Triangle");
}

BOOST_AUTO_TEST_CASE(testWriteTINStrict)
{
  TriangulatedSurface tin;
  tin.addTriangle(Triangle(Point(0, 0, 0), Point(1, 0, 0), Point(0.5, 1, 0)));

  GeoJSONOptions opts;
  opts.strict      = true;
  std::string json = writeGeoJSON(tin, opts);

  auto parsed = nlohmann::json::parse(json);
  BOOST_CHECK_EQUAL(parsed["type"], "MultiPolygon");
  BOOST_CHECK_EQUAL(parsed["properties"]["__sfcgal_original_type"],
                    "TriangulatedSurface");
}

BOOST_AUTO_TEST_CASE(testWriteTINExtension)
{
  TriangulatedSurface tin;
  tin.addTriangle(Triangle(Point(0, 0, 0), Point(1, 0, 0), Point(0.5, 1, 0)));

  GeoJSONOptions opts;
  opts.strict      = false;
  std::string json = writeGeoJSON(tin, opts);

  auto parsed = nlohmann::json::parse(json);
  BOOST_CHECK_EQUAL(parsed["type"], "TIN");
}

BOOST_AUTO_TEST_CASE(testWritePrecision)
{
  Point const pt(1.123456789, 2.987654321);

  GeoJSONOptions opts;
  opts.precision   = 2;
  std::string json = writeGeoJSON(pt, opts);

  auto parsed = nlohmann::json::parse(json);
  BOOST_CHECK_CLOSE(parsed["coordinates"][0].get<double>(), 1.12, 0.01);
  BOOST_CHECK_CLOSE(parsed["coordinates"][1].get<double>(), 2.99, 0.01);
}

BOOST_AUTO_TEST_CASE(testWritePreparedWithSRID)
{
  auto             pt = std::make_unique<Point>(1.0, 2.0);
  PreparedGeometry prepared(std::move(pt), 4326);

  std::string json   = writeGeoJSON(prepared);
  auto        parsed = nlohmann::json::parse(json);

  BOOST_CHECK(parsed.contains("crs"));
  BOOST_CHECK_EQUAL(parsed["crs"]["properties"]["name"],
                    "urn:ogc:def:crs:EPSG::4326");
}

BOOST_AUTO_TEST_CASE(testWriteWithBbox)
{
  LineString lineString;
  lineString.addPoint(Point(0, 0));
  lineString.addPoint(Point(10, 10));

  GeoJSONOptions opts;
  opts.includeBbox = true;
  std::string json = writeGeoJSON(lineString, opts);

  auto parsed = nlohmann::json::parse(json);
  BOOST_CHECK(parsed.contains("bbox"));
  BOOST_CHECK_EQUAL(parsed["bbox"].size(), 4U); // 2D bbox
  BOOST_CHECK_CLOSE(parsed["bbox"][0].get<double>(), 0.0, 0.0001);
  BOOST_CHECK_CLOSE(parsed["bbox"][1].get<double>(), 0.0, 0.0001);
  BOOST_CHECK_CLOSE(parsed["bbox"][2].get<double>(), 10.0, 0.0001);
  BOOST_CHECK_CLOSE(parsed["bbox"][3].get<double>(), 10.0, 0.0001);
}

BOOST_AUTO_TEST_CASE(testWriteGeometryCollection)
{
  GeometryCollection gc;
  gc.addGeometry(new Point(0, 0));
  gc.addGeometry(new Point(1, 1));

  std::string json   = writeGeoJSON(gc);
  auto        parsed = nlohmann::json::parse(json);

  BOOST_CHECK_EQUAL(parsed["type"], "GeometryCollection");
  BOOST_CHECK_EQUAL(parsed["geometries"].size(), 2U);
}

// ============================================================================
// ROUND-TRIP TESTS
// ============================================================================

BOOST_AUTO_TEST_CASE(testRoundTripPoint2D)
{
  Point const original(1.5, 2.5);
  std::string json   = writeGeoJSON(original);
  auto        parsed = readGeoJSON(json);

  BOOST_CHECK(parsed->is<Point>());
  const auto &result = parsed->as<Point>();
  BOOST_CHECK_CLOSE(CGAL::to_double(result.x()), 1.5, 0.0001);
  BOOST_CHECK_CLOSE(CGAL::to_double(result.y()), 2.5, 0.0001);
}

BOOST_AUTO_TEST_CASE(testRoundTripPoint3D)
{
  Point const original(1.5, 2.5, 3.5);
  std::string json   = writeGeoJSON(original);
  auto        parsed = readGeoJSON(json);

  const auto &result = parsed->as<Point>();
  BOOST_CHECK(result.is3D());
  BOOST_CHECK_CLOSE(CGAL::to_double(result.x()), 1.5, 0.0001);
  BOOST_CHECK_CLOSE(CGAL::to_double(result.y()), 2.5, 0.0001);
  BOOST_CHECK_CLOSE(CGAL::to_double(result.z()), 3.5, 0.0001);
}

BOOST_AUTO_TEST_CASE(testRoundTripPolygonWithHole)
{
  LineString exterior;
  exterior.addPoint(Point(0, 0));
  exterior.addPoint(Point(10, 0));
  exterior.addPoint(Point(10, 10));
  exterior.addPoint(Point(0, 10));
  exterior.addPoint(Point(0, 0));

  LineString interior;
  interior.addPoint(Point(2, 2));
  interior.addPoint(Point(8, 2));
  interior.addPoint(Point(8, 8));
  interior.addPoint(Point(2, 8));
  interior.addPoint(Point(2, 2));

  Polygon original(exterior);
  original.addInteriorRing(interior);

  std::string json   = writeGeoJSON(original);
  auto        parsed = readGeoJSON(json);

  BOOST_CHECK(parsed->is<Polygon>());
  const auto &result = parsed->as<Polygon>();
  BOOST_CHECK_EQUAL(result.numInteriorRings(), 1U);
}

BOOST_AUTO_TEST_CASE(testRoundTripTINExtension)
{
  TriangulatedSurface original;
  original.addTriangle(
      Triangle(Point(0, 0, 0), Point(1, 0, 0), Point(0.5, 1, 0)));
  original.addTriangle(
      Triangle(Point(1, 0, 0), Point(2, 0, 0), Point(1.5, 1, 0)));

  GeoJSONOptions opts;
  opts.strict = false; // Use extension mode for round-trip

  std::string json   = writeGeoJSON(original, opts);
  auto        parsed = readGeoJSON(json);

  BOOST_CHECK(parsed->is<TriangulatedSurface>());
  const auto &result = parsed->as<TriangulatedSurface>();
  BOOST_CHECK_EQUAL(result.numTriangles(), 2U);
}

// ============================================================================
// HELPER FUNCTION TESTS
// ============================================================================

BOOST_AUTO_TEST_CASE(testIsNativeGeoJSONType)
{
  BOOST_CHECK(isNativeGeoJSONType(TYPE_POINT));
  BOOST_CHECK(isNativeGeoJSONType(TYPE_LINESTRING));
  BOOST_CHECK(isNativeGeoJSONType(TYPE_POLYGON));
  BOOST_CHECK(isNativeGeoJSONType(TYPE_MULTIPOINT));
  BOOST_CHECK(isNativeGeoJSONType(TYPE_MULTILINESTRING));
  BOOST_CHECK(isNativeGeoJSONType(TYPE_MULTIPOLYGON));
  BOOST_CHECK(isNativeGeoJSONType(TYPE_GEOMETRYCOLLECTION));

  BOOST_CHECK(!isNativeGeoJSONType(TYPE_TRIANGLE));
  BOOST_CHECK(!isNativeGeoJSONType(TYPE_TRIANGULATEDSURFACE));
  BOOST_CHECK(!isNativeGeoJSONType(TYPE_POLYHEDRALSURFACE));
  BOOST_CHECK(!isNativeGeoJSONType(TYPE_SOLID));
  BOOST_CHECK(!isNativeGeoJSONType(TYPE_MULTISOLID));
}

BOOST_AUTO_TEST_CASE(testGeoJSONTypeName)
{
  Point const pt(0, 0);
  BOOST_CHECK_EQUAL(geoJSONTypeName(pt, true), "Point");
  BOOST_CHECK_EQUAL(geoJSONTypeName(pt, false), "Point");

  Triangle const tri(Point(0, 0), Point(1, 0), Point(0.5, 1));
  BOOST_CHECK_EQUAL(geoJSONTypeName(tri, true), "Polygon");
  BOOST_CHECK_EQUAL(geoJSONTypeName(tri, false), "Triangle");

  TriangulatedSurface tin;
  BOOST_CHECK_EQUAL(geoJSONTypeName(tin, true), "MultiPolygon");
  BOOST_CHECK_EQUAL(geoJSONTypeName(tin, false), "TriangulatedSurface");
}

// ============================================================================
// ERROR HANDLING TESTS
// ============================================================================

BOOST_AUTO_TEST_CASE(testInvalidJSON)
{
  BOOST_CHECK_THROW(readGeoJSON(std::string("not valid json")), std::exception);
}

BOOST_AUTO_TEST_CASE(testMissingType)
{
  BOOST_CHECK_THROW(readGeoJSON(std::string(R"({"coordinates": [1, 2]})")),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(testMissingCoordinates)
{
  BOOST_CHECK_THROW(readGeoJSON(std::string(R"({"type": "Point"})")),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(testUnknownType)
{
  BOOST_CHECK_THROW(
      readGeoJSON(std::string(R"({"type": "Unknown", "coordinates": []})")),
      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(testFeatureWithoutGeometry)
{
  BOOST_CHECK_THROW(
      readGeoJSON(std::string(R"({"type": "Feature", "properties": {}})")),
      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(testFeatureWithNullGeometry)
{
  BOOST_CHECK_THROW(
      readGeoJSON(std::string(
          R"({"type": "Feature", "geometry": null, "properties": {}})")),
      std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
