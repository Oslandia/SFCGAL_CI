/**
 * io.c - Implementation of input/output operations
 */

#include "io.h"
#include "util.h"
#include <ctype.h>
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INITIAL_BUFFER_SIZE 8192

// WKT Geometry type prefixes
static const char *const WKT_PREFIXES[] = {"POINT",
                                           "LINESTRING",
                                           "POLYGON",
                                           "MULTIPOINT",
                                           "MULTILINESTRING",
                                           "MULTIPOLYGON",
                                           "GEOMETRYCOLLECTION",
                                           "TRIANGLE",
                                           "POLYHEDRALSURFACE",
                                           "TIN",
                                           "SOLID",
                                           "MULTISOLID"};

static const size_t WKT_PREFIXES_COUNT =
    sizeof(WKT_PREFIXES) / sizeof(WKT_PREFIXES[0]);

/**
 * Remove UTF-8 BOM from the beginning of a buffer
 *
 * @param buffer Buffer to process
 * @param size Pointer to buffer size (will be updated if BOM is removed)
 */
static void
remove_utf8_bom(char **buffer, size_t *size)
{
  if (!buffer || !*buffer || !size || *size < 3) {
    return;
  }

  // Check for UTF-8 BOM (0xEF 0xBB 0xBF)
  unsigned char *data = (unsigned char *)*buffer;
  if (data[0] == 0xEF && data[1] == 0xBB && data[2] == 0xBF) {
    // Move the content 3 bytes to the left
    memmove(*buffer, *buffer + 3, *size - 3);
    *size -= 3;
    // Null-terminate
    (*buffer)[*size] = '\0';
  }
}

/**
 * Trim whitespace and control characters from both ends of a string
 *
 * @param str String to trim (modified in place)
 * @return Pointer to the trimmed string (same as input)
 */
static char *
trim_string(char *str)
{
  if (!str) {
    return str;
  }

  size_t len = strlen(str);
  if (len == 0) {
    return str;
  }

  // Trim trailing whitespace and control characters
  char *end = str + len - 1;
  while (end >= str &&
         (isspace((unsigned char)*end) || iscntrl((unsigned char)*end))) {
    *end = '\0';
    end--;
  }

  // Trim leading whitespace and control characters
  char *start = str;
  while (*start &&
         (isspace((unsigned char)*start) || iscntrl((unsigned char)*start))) {
    start++;
  }

  // Move string to the beginning if needed
  if (start != str) {
    size_t new_len = strlen(start);
    memmove(str, start, new_len + 1);
  }

  return str;
}

/**
 * Check if a string starts with any WKT geometry type
 *
 * @param str String to check
 * @return true if it's a WKT string, false otherwise
 */
static bool
is_wkt(const char *str)
{
  if (!str) {
    return false;
  }

  // Skip leading whitespace and control characters
  while (*str &&
         (isspace((unsigned char)*str) || iscntrl((unsigned char)*str))) {
    str++;
  }

  // Check each geometry type prefix
  for (size_t i = 0; i < WKT_PREFIXES_COUNT; i++) {
    if (strncasecmp(str, WKT_PREFIXES[i], strlen(WKT_PREFIXES[i])) == 0) {
      return true;
    }
  }

  return false;
}

/**
 * Check if a string is a hexadecimal representation of WKB
 *
 * @param str String to check
 * @return true if it's a hex WKB string, false otherwise
 */
static bool
is_hex_wkb(const char *str)
{
  if (!str) {
    return false;
  }

  // Skip leading whitespace
  while (*str && isspace((unsigned char)*str)) {
    str++;
  }

  // Empty string is not valid hex WKB
  if (*str == '\0') {
    return false;
  }

  // Check if string contains only hex characters and whitespace
  bool has_hex = false;
  while (*str) {
    if (isxdigit((unsigned char)*str)) {
      has_hex = true;
    } else if (!isspace((unsigned char)*str)) {
      return false;
    }
    str++;
  }

  return has_hex;
}

/**
 * Remove whitespace from a string in-place
 *
 * @param str String to clean
 */
static void
remove_whitespace(char *str)
{
  if (!str) {
    return;
  }

  size_t i = 0, j = 0;
  while (str[j]) {
    if (!isspace((unsigned char)str[j])) {
      str[i++] = str[j];
    }
    j++;
  }
  str[i] = '\0';
}

/**
 * Load a geometry from a WKT string
 *
 * @param wkt WKT string
 * @param geometry Pointer to store the loaded geometry
 * @return true on success, false on failure
 */
static bool
load_from_wkt(const char *wkt, sfcgal_geometry_t **geometry)
{
  if (!wkt || !geometry) {
    return false;
  }

  *geometry = sfcgal_io_read_wkt(wkt, strlen(wkt));
  return (*geometry != NULL);
}

/**
 * Load a geometry from a WKB hex string
 *
 * @param hex_wkb WKB hex string
 * @param geometry Pointer to store the loaded geometry
 * @return true on success, false on failure
 */
static bool
load_from_hex_wkb(const char *hex_wkb, sfcgal_geometry_t **geometry)
{
  if (!hex_wkb || !geometry) {
    return false;
  }

  // Copy the string so we can modify it
  char *cleaned_hex = safe_strdup(hex_wkb);
  if (!cleaned_hex) {
    return false;
  }

  // Remove whitespace
  remove_whitespace(cleaned_hex);

  // Load the geometry
  *geometry = sfcgal_io_read_wkb(cleaned_hex, strlen(cleaned_hex));

  free(cleaned_hex);
  return (*geometry != NULL);
}

/**
 * Safely grow a buffer to accommodate more data
 *
 * @param buffer Pointer to buffer pointer
 * @param current_size Current buffer size
 * @param required_size Required buffer size
 * @return true if buffer was successfully grown, false on failure
 */
static bool
grow_buffer(void **buffer, size_t *current_size, size_t required_size)
{
  if (!buffer || !current_size || *current_size >= required_size) {
    return false;
  }

  // Calculate new size (double the current size)
  size_t new_size = *current_size * 2;

  // Ensure new size is at least as large as required
  if (new_size < required_size) {
    new_size = required_size;
  }

  // Reallocate buffer
  void *new_buffer = realloc(*buffer, new_size);
  if (!new_buffer) {
    return false;
  }

  *buffer       = new_buffer;
  *current_size = new_size;
  return true;
}

/**
 * Load contents from a file into a buffer
 *
 * @param filename Path to file
 * @param buffer Pointer to store the buffer
 * @param size Pointer to store the buffer size
 * @return true on success, false on failure
 */
static bool
load_file_contents(const char *filename, void **buffer, size_t *size)
{
  if (!filename || !buffer || !size) {
    return false;
  }

  FILE *file = fopen(filename, "rb");
  if (!file) {
    fprintf(stderr, "Error opening file '%s': %s\n", filename, strerror(errno));
    return false;
  }

  // Get file size
  if (fseek(file, 0, SEEK_END) != 0) {
    fclose(file);
    return false;
  }

  long file_size = ftell(file);
  if (file_size < 0) {
    fclose(file);
    return false;
  }

  if (fseek(file, 0, SEEK_SET) != 0) {
    fclose(file);
    return false;
  }

  // Allocate buffer
  *buffer = malloc((size_t)file_size + 1); // +1 for null terminator
  if (!*buffer) {
    fclose(file);
    return false;
  }

  // Read file contents
  size_t bytes_read = fread(*buffer, 1, (size_t)file_size, file);
  fclose(file);

  if (bytes_read != (size_t)file_size) {
    free(*buffer);
    *buffer = NULL;
    return false;
  }

  // Null-terminate buffer
  ((char *)*buffer)[bytes_read] = '\0';
  *size                         = bytes_read;

  return true;
}

/**
 * Load a geometry from a file containing WKT
 *
 * @param filename Path to file
 * @param geometry Pointer to store the loaded geometry
 * @return true on success, false on failure
 */
static bool
load_from_file_wkt(const char *filename, sfcgal_geometry_t **geometry)
{
  if (!filename || !geometry) {
    return false;
  }

  char  *buffer = NULL;
  size_t size   = 0;

  if (!load_file_contents(filename, (void **)&buffer, &size)) {
    return false;
  }

  // Remove UTF-8 BOM if present
  remove_utf8_bom(&buffer, &size);

  // Clean up line breaks and trim
  for (size_t i = 0; i < size; i++) {
    if (buffer[i] == '\r' || buffer[i] == '\n') {
      buffer[i] = ' ';
    }
  }

  // Trim whitespace from both ends
  trim_string(buffer);

  // Load geometry from WKT
  bool success = is_wkt(buffer) && load_from_wkt(buffer, geometry);

  free(buffer);
  return success;
}

/**
 * Load a geometry from a file containing WKB (binary)
 *
 * @param filename Path to file
 * @param geometry Pointer to store the loaded geometry
 * @return true on success, false on failure
 */
static bool
load_from_file_wkb(const char *filename, sfcgal_geometry_t **geometry)
{
  if (!filename || !geometry) {
    return false;
  }

  void  *buffer = NULL;
  size_t size   = 0;

  if (!load_file_contents(filename, &buffer, &size)) {
    return false;
  }

  // Load geometry from WKB
  *geometry = sfcgal_io_read_wkb((const char *)buffer, size);

  free(buffer);
  return (*geometry != NULL);
}

/**
 * Load content from stdin into a buffer
 *
 * @param buffer Pointer to store the buffer
 * @param size Pointer to store the buffer size
 * @param binary Whether to read in binary mode
 * @return true on success, false on failure
 */
static bool
load_stdin_contents(void **buffer, size_t *size, bool binary)
{
  if (!buffer || !size) {
    return false;
  }

  size_t buffer_size = INITIAL_BUFFER_SIZE;
  *buffer            = malloc(buffer_size);
  if (!*buffer) {
    return false;
  }

  size_t bytes_read = 0;
  int    c;

  // Read input until EOF
  while ((c = getchar()) != EOF) {
    // Check if buffer needs to grow
    if (bytes_read >= buffer_size - 1) {
      if (!grow_buffer(buffer, &buffer_size, bytes_read + 2)) {
        free(*buffer);
        *buffer = NULL;
        return false;
      }
    }

    ((char *)*buffer)[bytes_read++] = (char)c;
  }

  if (!binary) {
    // Null-terminate for text mode
    ((char *)*buffer)[bytes_read] = '\0';
  }

  *size = bytes_read;
  return true;
}

/**
 * Load a geometry from stdin (WKT)
 *
 * @param geometry Pointer to store the loaded geometry
 * @return true on success, false on failure
 */
static bool
load_from_stdin_wkt(sfcgal_geometry_t **geometry)
{
  if (!geometry) {
    return false;
  }

  char  *buffer = NULL;
  size_t size   = 0;

  if (!load_stdin_contents((void **)&buffer, &size, false)) {
    return false;
  }

  // Remove UTF-8 BOM if present
  remove_utf8_bom(&buffer, &size);

  // Trim whitespace and control characters from both ends
  trim_string(buffer);

  // Load geometry from WKT
  bool success = load_from_wkt(buffer, geometry);

  free(buffer);
  return success;
}

/**
 * Load a geometry from stdin (WKB binary)
 *
 * @param geometry Pointer to store the loaded geometry
 * @return true on success, false on failure
 */
static bool
load_from_stdin_wkb(sfcgal_geometry_t **geometry)
{
  if (!geometry) {
    return false;
  }

  void  *buffer = NULL;
  size_t size   = 0;

  if (!load_stdin_contents(&buffer, &size, true)) {
    return false;
  }

  // Load geometry from WKB
  *geometry = sfcgal_io_read_wkb((const char *)buffer, size);

  free(buffer);
  return (*geometry != NULL);
}

/**
 * Load a geometry from a source
 */
bool
load_geometry(const char *source, sfcgal_geometry_t **geometry)
{
  if (!source || !geometry) {
    return false;
  }

  // Initialize output parameter
  *geometry = NULL;

  // Check if source is a WKT string
  if (is_wkt(source)) {
    return load_from_wkt(source, geometry);
  }

  // Check if source is a hex WKB string
  if (is_hex_wkb(source)) {
    return load_from_hex_wkb(source, geometry);
  }

  // Check if source is "stdin"
  if (strcmp(source, "stdin") == 0) {
    return load_from_stdin_wkt(geometry);
  }

  // Check if source is "stdin.wkb"
  if (strcmp(source, "stdin.wkb") == 0) {
    return load_from_stdin_wkb(geometry);
  }

  // Try to load from file (first as WKT, then as WKB)
  if (load_from_file_wkt(source, geometry)) {
    return true;
  }

  return load_from_file_wkb(source, geometry);
}

/**
 * Output a geometry as WKT
 *
 * @param geometry Geometry to output
 * @param precision Precision for floating-point coordinates
 * @return true on success, false on failure
 */
static bool
output_as_wkt(const sfcgal_geometry_t *geometry, int precision)
{
  if (!geometry) {
    return false;
  }

  char  *wkt = NULL;
  size_t len = 0;

  sfcgal_geometry_as_text_decim(geometry, precision, &wkt, &len);

  if (!wkt) {
    return false;
  }

  printf("%s\n", wkt);
  sfcgal_free_buffer(wkt);

  return true;
}

/**
 * Output a geometry as WKB (hex)
 *
 * @param geometry Geometry to output
 * @return true on success, false on failure
 */
static bool
output_as_wkb(const sfcgal_geometry_t *geometry)
{
  if (!geometry) {
    return false;
  }

  char  *wkb = NULL;
  size_t len = 0;

  sfcgal_geometry_as_hexwkb(geometry, &wkb, &len);

  if (!wkb) {
    return false;
  }

  printf("%s\n", wkb);
  sfcgal_free_buffer(wkb);

  return true;
}

/**
 * Output a geometry as plain text (SFCGAL native format with exact fractions)
 *
 * @param geometry Geometry to output
 * @param precision Precision for floating-point coordinates (ignored for SFCGAL
 * format)
 * @return true on success, false on failure
 */
static bool
output_as_txt(const sfcgal_geometry_t *geometry, int precision)
{
  if (!geometry) {
    return false;
  }

  // Ignore precision parameter for SFCGAL native format
  (void)precision;

  char  *txt = NULL;
  size_t len = 0;

  // Use SFCGAL's native format with exact fractions
  sfcgal_geometry_as_text(geometry, &txt, &len);

  if (!txt) {
    return false;
  }

  printf("%s\n", txt);
  sfcgal_free_buffer(txt);

  return true;
}

/**
 * Get GeoJSON type string for a WKT type
 *
 * @param wkt_type WKT geometry type string
 * @return GeoJSON type string
 */
static const char *
get_geojson_type(const char *wkt_type)
{
  if (!wkt_type) {
    return "Unknown";
  }

  struct {
    const char *wkt_prefix;
    const char *json_type;
  } type_map[] = {{"POINT", "Point"},
                  {"LINESTRING", "LineString"},
                  {"POLYGON", "Polygon"},
                  {"MULTIPOINT", "MultiPoint"},
                  {"MULTILINESTRING", "MultiLineString"},
                  {"MULTIPOLYGON", "MultiPolygon"},
                  {"GEOMETRYCOLLECTION", "GeometryCollection"},
                  {NULL, "Unknown"}};

  for (int i = 0; type_map[i].wkt_prefix != NULL; i++) {
    if (strncasecmp(wkt_type, type_map[i].wkt_prefix,
                    strlen(type_map[i].wkt_prefix)) == 0) {
      return type_map[i].json_type;
    }
  }

  return "Unknown";
}

/**
 * Output a geometry as GeoJSON
 *
 * @param geometry Geometry to output
 * @param precision Precision for floating-point coordinates
 * @return true on success, false on failure
 */
static bool
output_as_geojson(const sfcgal_geometry_t *geometry, int precision)
{
  if (!geometry) {
    return false;
  }

  char  *wkt = NULL;
  size_t len = 0;

  // Get WKT representation first
  sfcgal_geometry_as_text_decim(geometry, precision, &wkt, &len);

  if (!wkt) {
    return false;
  }

  // This is a basic GeoJSON conversion - for production a proper library should
  // be used
  const char *geojson_type = get_geojson_type(wkt);

  printf("{\n");
  printf("  \"type\": \"Feature\",\n");
  printf("  \"geometry\": {\n");
  printf("    \"type\": \"%s\",\n", geojson_type);

  // For a proper implementation, we would need to parse the WKT and convert it
  // to GeoJSON coordinates
  printf("    \"coordinates\": \"Detailed coordinates conversion not "
         "implemented\"\n");
  printf("  },\n");
  printf("  \"properties\": {}\n");
  printf("}\n");

  sfcgal_free_buffer(wkt);
  return true;
}

/**
 * Output a geometry in the specified format
 */
bool
output_geometry(const sfcgal_geometry_t *geometry, OutputFormat format,
                int precision)
{
  if (!geometry || precision < 0) {
    return false;
  }

  switch (format) {
  case FORMAT_WKT:
    return output_as_wkt(geometry, precision);
  case FORMAT_WKB:
    return output_as_wkb(geometry);
  case FORMAT_TXT:
    return output_as_txt(geometry, precision);
  case FORMAT_GEOJSON:
    return output_as_geojson(geometry, precision);
  default:
    fprintf(stderr, "Unknown output format\n");
    return false;
  }
}

/**
 * Convert string to output format
 */
bool
parse_output_format(const char *format_str, OutputFormat *format)
{
  if (!format_str || !format) {
    return false;
  }

  if (strcasecmp(format_str, "wkt") == 0) {
    *format = FORMAT_WKT;
  } else if (strcasecmp(format_str, "wkb") == 0) {
    *format = FORMAT_WKB;
  } else if (strcasecmp(format_str, "txt") == 0) {
    *format = FORMAT_TXT;
  } else if (strcasecmp(format_str, "geojson") == 0) {
    *format = FORMAT_GEOJSON;
  } else {
    return false;
  }

  return true;
}
