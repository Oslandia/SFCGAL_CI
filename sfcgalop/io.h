/**
 * io.h - Function declarations for input/output operations
 */

#ifndef IO_H
#define IO_H

#include "sfcgal_c.h"

/**
 * Output format enumeration
 */
typedef enum {
    FORMAT_WKT,
    FORMAT_WKB,
    FORMAT_TXT,
    FORMAT_GEOJSON
} OutputFormat;

/**
 * Load a geometry from a source
 * 
 * @param source Source string (WKT, WKB, file path, "stdin", "stdin.wkb")
 * @param geometry Pointer to geometry pointer that will be set
 * @return 0 on success, non-zero on failure
 */
int load_geometry(const char* source, sfcgal_geometry_t** geometry);

/**
 * Output a geometry in the specified format
 * 
 * @param geometry Geometry to output
 * @param format Output format
 * @param precision Precision for floating-point coordinates
 */
void output_geometry(const sfcgal_geometry_t* geometry, OutputFormat format, int precision);

#endif /* IO_H */
