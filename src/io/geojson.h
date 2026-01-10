// Copyright (c) 2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_IO_GEOJSON_H_
#define SFCGAL_IO_GEOJSON_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/config.h"

#include <memory>
#include <string>

#include <nlohmann/json.hpp>

namespace SFCGAL {
class PreparedGeometry;
} // namespace SFCGAL

namespace SFCGAL::io {

/**
 * @brief Options for GeoJSON serialization
 */
struct SFCGAL_API GeoJSONOptions {
  /**
   * If true, output strictly RFC 7946 compliant GeoJSON.
   * Non-standard types (TIN, Solid, etc.) are converted to standard types.
   * If false, use SFCGAL type names as extensions (allows round-trip).
   * Default: true
   */
  bool strict = true;

  /**
   * Number of decimal places for coordinates.
   * -1 means full precision (no rounding).
   *  0 means round to integers.
   *  N means round to N decimal places.
   * Default: -1
   */
  int precision = -1;

  /**
   * Include bounding box in output.
   * Default: false
   */
  bool includeBbox = false;
};

/**
 * @brief Read a Geometry from a GeoJSON string.
 * Supports Geometry, Feature, and FeatureCollection.
 * For Feature/FeatureCollection, only the geometry part is extracted.
 *
 * @param json GeoJSON string
 * @return Parsed geometry
 * @throws std::runtime_error on parse error
 */
SFCGAL_API auto
readGeoJSON(const std::string &json) -> std::unique_ptr<Geometry>;

/**
 * @brief Read a Geometry from a GeoJSON char array.
 * Supports Geometry, Feature, and FeatureCollection.
 * For Feature/FeatureCollection, only the geometry part is extracted.
 *
 * @param str The GeoJSON character array
 * @param len The length of the character array
 * @return Parsed geometry
 * @throws std::runtime_error on parse error
 */
SFCGAL_API auto
readGeoJSON(const char *str, size_t len) -> std::unique_ptr<Geometry>;

/**
 * @brief Read a Geometry from a parsed JSON object.
 *
 * @param json nlohmann::json object
 * @return Parsed geometry
 * @throws std::runtime_error on invalid geometry
 */
SFCGAL_API auto
readGeoJSON(const nlohmann::json &json) -> std::unique_ptr<Geometry>;

/**
 * @brief Read a PreparedGeometry (Geometry + SRID) from GeoJSON.
 * If CRS is specified in GeoJSON, extracts the EPSG code as SRID.
 *
 * @param json GeoJSON string
 * @return PreparedGeometry with geometry and SRID
 */
SFCGAL_API auto
readGeoJSONPrepared(const std::string &json)
    -> std::unique_ptr<PreparedGeometry>;

/**
 * @brief Read a PreparedGeometry (Geometry + SRID) from a GeoJSON char array.
 * If CRS is specified in GeoJSON, extracts the EPSG code as SRID.
 *
 * @param str The GeoJSON character array
 * @param len The length of the character array
 * @return PreparedGeometry with geometry and SRID
 */
SFCGAL_API auto
readGeoJSONPrepared(const char *str, size_t len)
    -> std::unique_ptr<PreparedGeometry>;

/**
 * @brief Read a PreparedGeometry from a parsed JSON object.
 *
 * @param json nlohmann::json object
 * @return PreparedGeometry with geometry and SRID
 */
SFCGAL_API auto
readGeoJSONPrepared(const nlohmann::json &json)
    -> std::unique_ptr<PreparedGeometry>;

/**
 * @brief Write a Geometry to GeoJSON string.
 *
 * @param geometry Geometry to serialize
 * @param options Serialization options
 * @return GeoJSON string
 */
SFCGAL_API auto
writeGeoJSON(const Geometry &geometry, const GeoJSONOptions &options = {})
    -> std::string;

/**
 * @brief Write a Geometry to a JSON object.
 *
 * @param geometry Geometry to serialize
 * @param options Serialization options
 * @return nlohmann::json object
 */
SFCGAL_API auto
writeGeoJSONObject(const Geometry &geometry, const GeoJSONOptions &options = {})
    -> nlohmann::json;

/**
 * @brief Write a PreparedGeometry to GeoJSON string.
 * If SRID is non-zero, includes a CRS member (non-standard but useful).
 *
 * @param prepared PreparedGeometry to serialize
 * @param options Serialization options
 * @return GeoJSON string
 */
SFCGAL_API auto
writeGeoJSON(const PreparedGeometry &prepared,
             const GeoJSONOptions   &options = {}) -> std::string;

/**
 * @brief Check if a geometry type is natively supported by GeoJSON (RFC 7946).
 *
 * @param geomType Geometry type ID
 * @return true if type is Point, LineString, Polygon, Multi*, or
 * GeometryCollection
 */
SFCGAL_API auto
isNativeGeoJSONType(GeometryType geomType) -> bool;

/**
 * @brief Get the GeoJSON type name for a geometry.
 * In strict mode, returns the RFC 7946 equivalent type.
 * In extension mode, returns the SFCGAL type name.
 *
 * @param geometry Geometry
 * @param strict Use strict RFC 7946 types
 * @return Type name string
 */
SFCGAL_API auto
geoJSONTypeName(const Geometry &geometry, bool strict = true) -> std::string;

} // namespace SFCGAL::io

#endif // SFCGAL_IO_GEOJSON_H_
