// Copyright (c) 2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/io/geojson.h"

#include "SFCGAL/Envelope.h"
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
#include "SFCGAL/PreparedGeometry.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"

#include <cmath>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace SFCGAL::io {

namespace {

// ============================================================================
// READING HELPERS
// ============================================================================

auto
parsePoint(const nlohmann::json &coords) -> std::unique_ptr<Point>
{
  if (!coords.is_array() || coords.size() < 2) {
    throw std::runtime_error(
        "Invalid coordinate: expected array of at least 2 numbers");
  }

  double x = coords[0].get<double>();
  double y = coords[1].get<double>();

  // Validate coordinate values to prevent NaN or infinity
  if (std::isnan(x) || std::isnan(y) || std::isinf(x) || std::isinf(y)) {
    throw std::runtime_error(
        "Invalid coordinate: NaN or infinity values not allowed");
  }

  if (coords.size() >= 3 && !coords[2].is_null()) {
    double z = coords[2].get<double>();
    if (std::isnan(z) || std::isinf(z)) {
      throw std::runtime_error(
          "Invalid coordinate: NaN or infinity values not allowed");
    }
    return std::make_unique<Point>(x, y, z);
  }
  return std::make_unique<Point>(x, y);
}

auto
parseLineString(const nlohmann::json &coords) -> std::unique_ptr<LineString>
{
  auto lineString = std::make_unique<LineString>();
  for (const auto &coord : coords) {
    lineString->addPoint(*parsePoint(coord));
  }
  return lineString;
}

auto
parseLinearRing(const nlohmann::json &coords) -> LineString
{
  LineString ring;
  for (const auto &coord : coords) {
    ring.addPoint(*parsePoint(coord));
  }
  return ring;
}

auto
parsePolygon(const nlohmann::json &coords) -> std::unique_ptr<Polygon>
{
  if (coords.empty()) {
    return std::make_unique<Polygon>();
  }

  // Exterior ring
  auto exterior = parseLinearRing(coords[0]);
  auto poly     = std::make_unique<Polygon>(exterior);

  // Interior rings (holes)
  for (size_t i = 1; i < coords.size(); ++i) {
    auto interior = parseLinearRing(coords[i]);
    poly->addInteriorRing(interior);
  }

  return poly;
}

auto
parseTriangle(const nlohmann::json &coords) -> std::unique_ptr<Triangle>
{
  // Triangle in GeoJSON extension mode: [[p0, p1, p2, p0]]
  if (coords.empty() || coords[0].size() < 3) {
    throw std::runtime_error("Invalid Triangle: expected ring with at least 3 "
                             "points");
  }

  const auto &ring = coords[0];
  auto        pt0  = parsePoint(ring[0]);
  auto        pt1  = parsePoint(ring[1]);
  auto        pt2  = parsePoint(ring[2]);

  return std::make_unique<Triangle>(*pt0, *pt1, *pt2);
}

auto
parseMultiPoint(const nlohmann::json &coords) -> std::unique_ptr<MultiPoint>
{
  auto multiPoint = std::make_unique<MultiPoint>();
  for (const auto &coord : coords) {
    multiPoint->addGeometry(parsePoint(coord).release());
  }
  return multiPoint;
}

auto
parseMultiLineString(const nlohmann::json &coords)
    -> std::unique_ptr<MultiLineString>
{
  auto multiLineString = std::make_unique<MultiLineString>();
  for (const auto &lineCoords : coords) {
    multiLineString->addGeometry(parseLineString(lineCoords).release());
  }
  return multiLineString;
}

auto
parseMultiPolygon(const nlohmann::json &coords) -> std::unique_ptr<MultiPolygon>
{
  auto multiPolygon = std::make_unique<MultiPolygon>();
  for (const auto &polyCoords : coords) {
    multiPolygon->addGeometry(parsePolygon(polyCoords).release());
  }
  return multiPolygon;
}

auto
parseTIN(const nlohmann::json &coords) -> std::unique_ptr<TriangulatedSurface>
{
  auto tin = std::make_unique<TriangulatedSurface>();
  for (const auto &triCoords : coords) {
    // Each element is a polygon-like structure for a triangle
    if (!triCoords.empty() && triCoords[0].size() >= 3) {
      const auto &ring = triCoords[0];
      auto        pt0  = parsePoint(ring[0]);
      auto        pt1  = parsePoint(ring[1]);
      auto        pt2  = parsePoint(ring[2]);
      tin->addTriangle(Triangle(*pt0, *pt1, *pt2));
    }
  }
  return tin;
}

auto
parsePolyhedralSurface(const nlohmann::json &coords)
    -> std::unique_ptr<PolyhedralSurface>
{
  auto polyhedralSurface = std::make_unique<PolyhedralSurface>();
  for (const auto &polyCoords : coords) {
    polyhedralSurface->addPolygon(*parsePolygon(polyCoords));
  }
  return polyhedralSurface;
}

auto
parseSolid(const nlohmann::json &coords) -> std::unique_ptr<Solid>
{
  // Solid in extension mode: array of shells (each shell is a MultiPolygon)
  // First shell is exterior, rest are interior
  if (coords.empty()) {
    return std::make_unique<Solid>();
  }

  // Parse exterior shell
  auto exteriorShell = parsePolyhedralSurface(coords[0]);
  auto solid         = std::make_unique<Solid>(*exteriorShell);

  // Parse interior shells
  for (size_t i = 1; i < coords.size(); ++i) {
    auto interiorShell = parsePolyhedralSurface(coords[i]);
    solid->addInteriorShell(*interiorShell);
  }

  return solid;
}

// Forward declaration for recursive parsing
auto
parseGeometryObject(const nlohmann::json &json) -> std::unique_ptr<Geometry>;

auto
parseGeometryCollection(const nlohmann::json &geometries)
    -> std::unique_ptr<GeometryCollection>
{
  auto geometryCollection = std::make_unique<GeometryCollection>();
  for (const auto &geomJson : geometries) {
    geometryCollection->addGeometry(parseGeometryObject(geomJson).release());
  }
  return geometryCollection;
}

auto
parseGeometryObject(const nlohmann::json &json) -> std::unique_ptr<Geometry>
{
  if (!json.contains("type")) {
    throw std::runtime_error("GeoJSON geometry must have a 'type' member");
  }

  std::string type = json["type"].get<std::string>();

  // Handle GeometryCollection separately (has "geometries" not "coordinates")
  if (type == "GeometryCollection") {
    if (!json.contains("geometries")) {
      throw std::runtime_error(
          "GeometryCollection must have 'geometries' member");
    }
    return parseGeometryCollection(json["geometries"]);
  }

  // All other types have "coordinates"
  if (!json.contains("coordinates")) {
    throw std::runtime_error("GeoJSON geometry must have 'coordinates' member");
  }

  const auto &coords = json["coordinates"];

  // Standard GeoJSON types
  if (type == "Point") {
    return parsePoint(coords);
  }
  if (type == "LineString") {
    return parseLineString(coords);
  }
  if (type == "Polygon") {
    return parsePolygon(coords);
  }
  if (type == "MultiPoint") {
    return parseMultiPoint(coords);
  }
  if (type == "MultiLineString") {
    return parseMultiLineString(coords);
  }
  if (type == "MultiPolygon") {
    return parseMultiPolygon(coords);
  }
  // SFCGAL extension types
  if (type == "Triangle") {
    return parseTriangle(coords);
  }
  if (type == "TIN" || type == "TriangulatedSurface") {
    return parseTIN(coords);
  }
  if (type == "PolyhedralSurface") {
    return parsePolyhedralSurface(coords);
  }
  if (type == "Solid") {
    return parseSolid(coords);
  }

  throw std::runtime_error("Unknown geometry type: " + type);
}

auto
extractSRIDFromCRS(const nlohmann::json &json) -> srid_t
{
  if (!json.contains("crs") || json["crs"].is_null()) {
    return 0;
  }

  const auto &crs = json["crs"];
  if (!crs.contains("properties") || !crs["properties"].contains("name")) {
    return 0;
  }

  std::string name = crs["properties"]["name"].get<std::string>();

  // Parse EPSG code from URN or simple format
  // Formats: "urn:ogc:def:crs:EPSG::4326", "EPSG:4326", "epsg:4326"
  std::regex const epsgRegex(R"((?:urn:ogc:def:crs:)?EPSG::?(\d+))",
                             std::regex::icase);
  std::smatch      match;
  if (std::regex_search(name, match, epsgRegex)) {
    return static_cast<srid_t>(std::stoi(match[1].str()));
  }

  return 0;
}

// ============================================================================
// WRITING HELPERS
// ============================================================================

auto
roundCoordinate(double value, int precision) -> double
{
  if (precision < 0) {
    return value;
  }
  double const factor = std::pow(10.0, precision);
  return std::round(value * factor) / factor;
}

auto
coordinateToJson(const Point &point, int precision) -> nlohmann::json
{
  nlohmann::json coord = nlohmann::json::array();

  coord.push_back(roundCoordinate(CGAL::to_double(point.x()), precision));
  coord.push_back(roundCoordinate(CGAL::to_double(point.y()), precision));
  if (point.is3D()) {
    coord.push_back(roundCoordinate(CGAL::to_double(point.z()), precision));
  }

  return coord;
}

auto
lineStringToCoords(const LineString &lineString, int precision)
    -> nlohmann::json
{
  nlohmann::json coords = nlohmann::json::array();
  if (lineString.isEmpty()) {
    return coords; // Return empty array for empty LineString
  }
  for (size_t i = 0; i < lineString.numPoints(); ++i) {
    coords.push_back(coordinateToJson(lineString.pointN(i), precision));
  }
  return coords;
}

auto
polygonToCoords(const Polygon &poly, int precision) -> nlohmann::json
{
  nlohmann::json coords = nlohmann::json::array();

  if (poly.isEmpty()) {
    return coords; // Return empty array for empty Polygon
  }

  // Exterior ring
  coords.push_back(lineStringToCoords(poly.exteriorRing(), precision));

  // Interior rings
  for (size_t i = 0; i < poly.numInteriorRings(); ++i) {
    coords.push_back(lineStringToCoords(poly.interiorRingN(i), precision));
  }

  return coords;
}

auto
triangleToCoords(const Triangle &tri, int precision) -> nlohmann::json
{
  nlohmann::json coords = nlohmann::json::array();
  if (tri.isEmpty()) {
    return coords; // Return empty array for empty Triangle
  }

  // Triangle as polygon: [[p0, p1, p2, p0]]
  nlohmann::json ring = nlohmann::json::array();
  ring.push_back(coordinateToJson(tri.vertex(0), precision));
  ring.push_back(coordinateToJson(tri.vertex(1), precision));
  ring.push_back(coordinateToJson(tri.vertex(2), precision));
  ring.push_back(coordinateToJson(tri.vertex(0), precision)); // Close ring

  coords.push_back(ring);
  return coords;
}

// Forward declaration
auto
geometryToJsonObject(const Geometry &geom, const GeoJSONOptions &options)
    -> nlohmann::json;

auto
geometryCollectionToJson(const GeometryCollection &geometryCollection,
                         const GeoJSONOptions     &options) -> nlohmann::json
{
  nlohmann::json geometries = nlohmann::json::array();
  for (size_t i = 0; i < geometryCollection.numGeometries(); ++i) {
    geometries.push_back(
        geometryToJsonObject(geometryCollection.geometryN(i), options));
  }
  return geometries;
}

auto
multiPointToCoords(const MultiPoint &multiPoint, int precision)
    -> nlohmann::json
{
  nlohmann::json coords = nlohmann::json::array();
  if (multiPoint.isEmpty()) {
    return coords; // Return empty array for empty MultiPoint
  }
  for (size_t i = 0; i < multiPoint.numGeometries(); ++i) {
    coords.push_back(coordinateToJson(multiPoint.pointN(i), precision));
  }
  return coords;
}

auto
multiLineStringToCoords(const MultiLineString &multiLineString, int precision)
    -> nlohmann::json
{
  nlohmann::json coords = nlohmann::json::array();
  if (multiLineString.isEmpty()) {
    return coords; // Return empty array for empty MultiLineString
  }
  for (size_t i = 0; i < multiLineString.numGeometries(); ++i) {
    coords.push_back(
        lineStringToCoords(multiLineString.lineStringN(i), precision));
  }
  return coords;
}

auto
multiPolygonToCoords(const MultiPolygon &multiPolygon, int precision)
    -> nlohmann::json
{
  nlohmann::json coords = nlohmann::json::array();
  if (multiPolygon.isEmpty()) {
    return coords; // Return empty array for empty MultiPolygon
  }
  for (size_t i = 0; i < multiPolygon.numGeometries(); ++i) {
    coords.push_back(polygonToCoords(multiPolygon.polygonN(i), precision));
  }
  return coords;
}

auto
tinToCoords(const TriangulatedSurface &tin, int precision) -> nlohmann::json
{
  nlohmann::json coords = nlohmann::json::array();
  if (tin.isEmpty()) {
    return coords; // Return empty array for empty TriangulatedSurface
  }
  for (size_t i = 0; i < tin.numTriangles(); ++i) {
    coords.push_back(triangleToCoords(tin.triangleN(i), precision));
  }
  return coords;
}

auto
polyhedralSurfaceToCoords(const PolyhedralSurface &polyhedralSurface,
                          int                      precision) -> nlohmann::json
{
  nlohmann::json coords = nlohmann::json::array();
  if (polyhedralSurface.isEmpty()) {
    return coords; // Return empty array for empty PolyhedralSurface
  }
  for (size_t i = 0; i < polyhedralSurface.numPolygons(); ++i) {
    coords.push_back(polygonToCoords(polyhedralSurface.polygonN(i), precision));
  }
  return coords;
}

auto
solidToCoords(const Solid &solid, int precision) -> nlohmann::json
{
  nlohmann::json coords = nlohmann::json::array();
  if (solid.isEmpty()) {
    return coords; // Return empty array for empty Solid
  }

  // Exterior shell
  coords.push_back(polyhedralSurfaceToCoords(solid.exteriorShell(), precision));

  // Interior shells
  for (size_t i = 0; i < solid.numInteriorShells(); ++i) {
    coords.push_back(
        polyhedralSurfaceToCoords(solid.interiorShellN(i), precision));
  }

  return coords;
}

auto
computeBbox(const Geometry &geom) -> nlohmann::json
{
  if (geom.isEmpty()) {
    return nullptr;
  }

  auto envelope = geom.envelope();
  if (envelope.isEmpty()) {
    return nullptr;
  }

  nlohmann::json bbox = nlohmann::json::array();
  bbox.push_back(CGAL::to_double(envelope.xMin()));
  bbox.push_back(CGAL::to_double(envelope.yMin()));
  if (envelope.is3D()) {
    bbox.push_back(CGAL::to_double(envelope.zMin()));
  }
  bbox.push_back(CGAL::to_double(envelope.xMax()));
  bbox.push_back(CGAL::to_double(envelope.yMax()));
  if (envelope.is3D()) {
    bbox.push_back(CGAL::to_double(envelope.zMax()));
  }

  return bbox;
}

// NOLINTBEGIN(readability-function-cognitive-complexity)
auto
geometryToJsonObject(const Geometry &geom, const GeoJSONOptions &options)
    -> nlohmann::json
{
  nlohmann::json result;
  int const      precision = options.precision;
  bool const     strict    = options.strict;

  switch (geom.geometryTypeId()) {
  case TYPE_POINT: {
    const auto &point = geom.as<Point>();
    result["type"]    = "Point";
    if (point.isEmpty()) {
      result["coordinates"] = nlohmann::json::array();
    } else {
      result["coordinates"] = coordinateToJson(point, precision);
    }
    break;
  }
  case TYPE_LINESTRING: {
    result["type"] = "LineString";
    result["coordinates"] =
        lineStringToCoords(geom.as<LineString>(), precision);
    break;
  }
  case TYPE_POLYGON: {
    result["type"]        = "Polygon";
    result["coordinates"] = polygonToCoords(geom.as<Polygon>(), precision);
    break;
  }
  case TYPE_TRIANGLE: {
    const auto &tri = geom.as<Triangle>();
    if (strict) {
      result["type"]        = "Polygon";
      result["coordinates"] = triangleToCoords(tri, precision);
      result["properties"]["__sfcgal_original_type"] = "Triangle";
    } else {
      result["type"]        = "Triangle";
      result["coordinates"] = triangleToCoords(tri, precision);
    }
    break;
  }
  case TYPE_MULTIPOINT: {
    result["type"] = "MultiPoint";
    result["coordinates"] =
        multiPointToCoords(geom.as<MultiPoint>(), precision);
    break;
  }
  case TYPE_MULTILINESTRING: {
    result["type"] = "MultiLineString";
    result["coordinates"] =
        multiLineStringToCoords(geom.as<MultiLineString>(), precision);
    break;
  }
  case TYPE_MULTIPOLYGON: {
    result["type"] = "MultiPolygon";
    result["coordinates"] =
        multiPolygonToCoords(geom.as<MultiPolygon>(), precision);
    break;
  }
  case TYPE_GEOMETRYCOLLECTION: {
    result["type"] = "GeometryCollection";
    result["geometries"] =
        geometryCollectionToJson(geom.as<GeometryCollection>(), options);
    break;
  }
  case TYPE_TRIANGULATEDSURFACE: {
    const auto &tin = geom.as<TriangulatedSurface>();
    if (strict) {
      result["type"]        = "MultiPolygon";
      result["coordinates"] = tinToCoords(tin, precision);
      result["properties"]["__sfcgal_original_type"] = "TriangulatedSurface";
    } else {
      result["type"]        = "TIN";
      result["coordinates"] = tinToCoords(tin, precision);
    }
    break;
  }
  case TYPE_POLYHEDRALSURFACE: {
    const auto &polyhedralSurface = geom.as<PolyhedralSurface>();
    if (strict) {
      result["type"] = "MultiPolygon";
      result["coordinates"] =
          polyhedralSurfaceToCoords(polyhedralSurface, precision);
      result["properties"]["__sfcgal_original_type"] = "PolyhedralSurface";
    } else {
      result["type"] = "PolyhedralSurface";
      result["coordinates"] =
          polyhedralSurfaceToCoords(polyhedralSurface, precision);
    }
    break;
  }
  case TYPE_SOLID: {
    const auto &solid = geom.as<Solid>();
    if (strict) {
      // In strict mode, export only exterior shell as MultiPolygon
      result["type"] = "MultiPolygon";
      result["coordinates"] =
          polyhedralSurfaceToCoords(solid.exteriorShell(), precision);
      result["properties"]["__sfcgal_original_type"] = "Solid";
    } else {
      result["type"]        = "Solid";
      result["coordinates"] = solidToCoords(solid, precision);
    }
    break;
  }
  case TYPE_MULTISOLID: {
    const auto &multiSolid = geom.as<MultiSolid>();
    if (strict) {
      // Convert to GeometryCollection of MultiPolygons
      result["type"]            = "GeometryCollection";
      nlohmann::json geometries = nlohmann::json::array();
      for (size_t i = 0; i < multiSolid.numGeometries(); ++i) {
        nlohmann::json solidJson;
        solidJson["type"]        = "MultiPolygon";
        solidJson["coordinates"] = polyhedralSurfaceToCoords(
            multiSolid.solidN(i).exteriorShell(), precision);
        geometries.push_back(solidJson);
      }
      result["geometries"]                           = geometries;
      result["properties"]["__sfcgal_original_type"] = "MultiSolid";
    } else {
      result["type"]        = "MultiSolid";
      nlohmann::json solids = nlohmann::json::array();
      for (size_t i = 0; i < multiSolid.numGeometries(); ++i) {
        solids.push_back(solidToCoords(multiSolid.solidN(i), precision));
      }
      result["coordinates"] = solids;
    }
    break;
  }
  case TYPE_NURBSCURVE: {
    // NURBS curves are converted to LineString
    const auto &nurbsCurve    = geom.as<NURBSCurve>();
    auto        lineStringPtr = nurbsCurve.toLineString();
    result["type"]            = "LineString";
    result["coordinates"]     = lineStringToCoords(*lineStringPtr, precision);
    result["properties"]["__sfcgal_original_type"] = "NURBSCurve";
    break;
  }
  default:
    throw std::runtime_error("Unsupported geometry type for GeoJSON export");
  }

  if (options.includeBbox) {
    auto bbox = computeBbox(geom);
    if (!bbox.is_null()) {
      result["bbox"] = bbox;
    }
  }

  return result;
}
// NOLINTEND(readability-function-cognitive-complexity)

} // anonymous namespace

// ============================================================================
// PUBLIC API IMPLEMENTATION
// ============================================================================

auto
readGeoJSON(const std::string &json) -> std::unique_ptr<Geometry>
{
  nlohmann::json parsed = nlohmann::json::parse(json);
  return readGeoJSON(parsed);
}

auto
readGeoJSON(const char *str, size_t len) -> std::unique_ptr<Geometry>
{
  std::string    json_str(str, len);
  nlohmann::json parsed = nlohmann::json::parse(json_str);
  return readGeoJSON(parsed);
}

auto
readGeoJSON(const nlohmann::json &json) -> std::unique_ptr<Geometry>
{
  if (!json.contains("type")) {
    throw std::runtime_error("Invalid GeoJSON: missing 'type' member");
  }

  std::string type = json["type"].get<std::string>();

  // Handle Feature - extract geometry
  if (type == "Feature") {
    if (!json.contains("geometry") || json["geometry"].is_null()) {
      throw std::runtime_error("Feature has no geometry");
    }
    return parseGeometryObject(json["geometry"]);
  }

  // Handle FeatureCollection - combine geometries
  if (type == "FeatureCollection") {
    if (!json.contains("features")) {
      throw std::runtime_error("FeatureCollection has no features");
    }
    auto geometryCollection = std::make_unique<GeometryCollection>();
    for (const auto &feature : json["features"]) {
      if (feature.contains("geometry") && !feature["geometry"].is_null()) {
        geometryCollection->addGeometry(
            parseGeometryObject(feature["geometry"]).release());
      }
    }
    return geometryCollection;
  }

  // Handle direct geometry object
  return parseGeometryObject(json);
}

auto
readGeoJSONPrepared(const std::string &json)
    -> std::unique_ptr<PreparedGeometry>
{
  nlohmann::json parsed = nlohmann::json::parse(json);
  return readGeoJSONPrepared(parsed);
}

auto
readGeoJSONPrepared(const char *str, size_t len)
    -> std::unique_ptr<PreparedGeometry>
{
  std::string    json_str(str, len);
  nlohmann::json parsed = nlohmann::json::parse(json_str);
  return readGeoJSONPrepared(parsed);
}

auto
readGeoJSONPrepared(const nlohmann::json &json)
    -> std::unique_ptr<PreparedGeometry>
{
  auto   geometry = readGeoJSON(json);
  srid_t srid     = extractSRIDFromCRS(json);
  return std::make_unique<PreparedGeometry>(std::move(geometry), srid);
}

auto
writeGeoJSON(const Geometry &geometry, const GeoJSONOptions &options)
    -> std::string
{
  nlohmann::json json = writeGeoJSONObject(geometry, options);
  return json.dump();
}

auto
writeGeoJSONObject(const Geometry &geometry, const GeoJSONOptions &options)
    -> nlohmann::json
{
  return geometryToJsonObject(geometry, options);
}

auto
writeGeoJSON(const PreparedGeometry &prepared, const GeoJSONOptions &options)
    -> std::string
{
  nlohmann::json json = writeGeoJSONObject(prepared.geometry(), options);

  // Add CRS if SRID is set
  if (prepared.SRID() != 0) {
    json["crs"] = {{"type", "name"},
                   {"properties",
                    {{"name", "urn:ogc:def:crs:EPSG::" +
                                  std::to_string(prepared.SRID())}}}};
  }

  return json.dump();
}

auto
isNativeGeoJSONType(GeometryType geomType) -> bool
{
  switch (geomType) {
  case TYPE_POINT:
  case TYPE_LINESTRING:
  case TYPE_POLYGON:
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_GEOMETRYCOLLECTION:
    return true;
  default:
    return false;
  }
}

auto
geoJSONTypeName(const Geometry &geometry, bool strict) -> std::string
{
  if (!strict) {
    return geometry.geometryType();
  }

  switch (geometry.geometryTypeId()) {
  case TYPE_POINT:
    return "Point";
  case TYPE_LINESTRING:
  case TYPE_NURBSCURVE:
    return "LineString";
  case TYPE_POLYGON:
  case TYPE_TRIANGLE:
    return "Polygon";
  case TYPE_MULTIPOINT:
    return "MultiPoint";
  case TYPE_MULTILINESTRING:
    return "MultiLineString";
  case TYPE_MULTIPOLYGON:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
  case TYPE_SOLID:
    return "MultiPolygon";
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_MULTISOLID:
  default:
    return "GeometryCollection";
  }
}

} // namespace SFCGAL::io
