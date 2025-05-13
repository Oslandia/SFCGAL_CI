/**
 * io.h - Function declarations for input/output operations
 */

#ifndef IO_H
#define IO_H

#include "sfcgal_c.h"
#include <stdbool.h>

/**
 * Output format enumeration
 */
typedef enum {
    FORMAT_WKT,      // Well-Known Text
    FORMAT_WKB,      // Well-Known Binary (hex)
    FORMAT_TXT,      // Plain text
    FORMAT_GEOJSON   // GeoJSON format
} OutputFormat;

/**
 * Load a geometry from a source
 * 
 * @param source Source string (WKT, WKB, file path, "stdin", "stdin.wkb")
 * @param geometry Pointer to geometry pointer that will be set
 * @return true on success, false on failure
 */
bool load_geometry(const char* source, sfcgal_geometry_t** geometry);

/**
 * Output a geometry in the specified format
 * 
 * @param geometry Geometry to output
 * @param format Output format
 * @param precision Precision for floating-point coordinates
 * @return true on success, false on failure
 */
bool output_geometry(const sfcgal_geometry_t* geometry, OutputFormat format, int precision);

/**
 * Convert string to output format
 * 
 * @param format_str Format string ("wkt", "wkb", "txt", "geojson")
 * @param format Pointer to store the output format
 * @return true if format was successfully parsed, false otherwise
 */
bool parse_output_format(const char* format_str, OutputFormat* format);

#endif /* IO_H */
