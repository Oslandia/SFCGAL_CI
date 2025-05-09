/**
 * io.c - Implementation of input/output operations
 */

#include "io.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define BUFFER_SIZE 8192

/**
 * Check if a string starts with "POINT", "LINESTRING", etc.
 */
static bool is_wkt(const char* str) {
    const char* wkt_types[] = {
        "POINT", "LINESTRING", "POLYGON", "MULTIPOINT", "MULTILINESTRING",
        "MULTIPOLYGON", "GEOMETRYCOLLECTION", "TRIANGLE", "POLYHEDRALSURFACE",
        "TIN", "SOLID", "MULTISOLID"
    };
    
    int n_types = sizeof(wkt_types) / sizeof(wkt_types[0]);
    while (*str && (isspace((unsigned char)*str) || !isprint((unsigned char)*str))) 
        str++;

    for (int i = 0; i < n_types; i++) {
        if (strncasecmp(str, wkt_types[i], strlen(wkt_types[i])) == 0) {
            return true;
        }
    }
    return false;
}

/**
 * Check if a string is a hexadecimal representation of WKB
 */
static bool is_hex_wkb(const char* str) {
    // Skip leading whitespace
    while (isspace((unsigned char)*str)) str++;
    
    // Check if string contains only hex characters
    while (*str) {
        if (!isxdigit((unsigned char)*str) && !isspace((unsigned char)*str)) {
            return false;
        }
        str++;
    }
    
    return true;
}

/**
 * Load a geometry from a WKT string
 */
static int load_from_wkt(const char* wkt, sfcgal_geometry_t** geometry) {
    *geometry = sfcgal_io_read_wkt(wkt, strlen(wkt));
    return (*geometry == NULL) ? -1 : 0;
}

/**
 * Load a geometry from a WKB string (hex)
 */
static int load_from_hex_wkb(const char* hex_wkb, sfcgal_geometry_t** geometry) {
    // Skip whitespace for accurate reading
    char* cleaned_hex = strdup(hex_wkb);
    if (!cleaned_hex) {
        return -1;
    }
    
    // Remove whitespace from hex string
    size_t i = 0, j = 0;
    while (hex_wkb[j]) {
        if (!isspace((unsigned char)hex_wkb[j])) {
            cleaned_hex[i++] = hex_wkb[j];
        }
        j++;
    }
    cleaned_hex[i] = '\0';
    
    // Now we can read the WKB directly using SFCGAL's function
    *geometry = sfcgal_io_read_wkb(cleaned_hex, strlen(cleaned_hex));
    
    free(cleaned_hex);
    return (*geometry == NULL) ? -1 : 0;
}

/**
 * Load a geometry from a file containing WKT
 */
static int load_from_file_wkt(const char* filename, sfcgal_geometry_t** geometry) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        return -1;
    }
    
    // Read file content into buffer
    char* buffer = (char*) malloc(BUFFER_SIZE);
    if (!buffer) {
        fclose(file);
        return -1;
    }
    
    size_t buffer_size = BUFFER_SIZE;
    size_t total_read = 0;
    size_t bytes_read;
    
    // Read the entire file, potentially growing the buffer
    while ((bytes_read = fread(buffer + total_read, 1, buffer_size - total_read - 1, file)) > 0) {
        total_read += bytes_read;
        
        // Check if we need to grow the buffer
        if (total_read >= buffer_size - 1) {
            buffer_size *= 2;
            char* new_buffer = (char*) realloc(buffer, buffer_size);
            if (!new_buffer) {
                free(buffer);
                fclose(file);
                return -1;
            }
            buffer = new_buffer;
        }
    }
    
    buffer[total_read] = '\0';
    fclose(file);
    
    // Now process the buffer to remove \r and \n characters
    size_t write_pos = 0;
    for (size_t read_pos = 0; read_pos < total_read; read_pos++) {
        if (buffer[read_pos] != '\r' && buffer[read_pos] != '\n') {
            buffer[write_pos++] = buffer[read_pos];
        }
    }
    buffer[write_pos] = '\0';
    
    // Now try to parse the WKT
    int result = -1;
    if( is_wkt(buffer)) {
        result = load_from_wkt(buffer, geometry);
    }
    free(buffer);
    
    return result;
}

/**
 * Load a geometry from a file containing WKB (binary)
 */
static int load_from_file_wkb(const char* filename, sfcgal_geometry_t** geometry) {
    FILE* file = fopen(filename, "rb");
    if (!file) {
        return -1;
    }
    
    // Get file size
    if (fseek(file, 0, SEEK_END) != 0) {
        fclose(file);
        return -1;
    }
    
    long file_size_long = ftell(file);
    if (file_size_long < 0) {
        fclose(file);
        return -1;
    }
    
    if (fseek(file, 0, SEEK_SET) != 0) {
        fclose(file);
        return -1;
    }
    
    // Safe conversion from long to size_t
    size_t file_size = (size_t)file_size_long;
    
    // Allocate buffer with explicit cast
    unsigned char* buffer = (unsigned char*)malloc(file_size);
    if (!buffer) {
        fclose(file);
        return -1;
    }
    
    size_t bytes_read = fread(buffer, 1, file_size, file);
    fclose(file);
    
    // Use SFCGAL's function to read WKB directly from binary
    *geometry = sfcgal_io_read_wkb((const char*)buffer, bytes_read);
    
    free(buffer);
    return (*geometry == NULL) ? -1 : 0;
}

/**
 * Load a geometry from stdin (WKT)
 */
static int load_from_stdin_wkt(sfcgal_geometry_t** geometry) {
    // Allocate initial buffer
    size_t buffer_size = BUFFER_SIZE;
    char* buffer = (char*) malloc(buffer_size);
    if (!buffer) {
        return -1;
    }
    
    size_t bytes_read = 0;
    int c;
    
    // Read input until EOF
    while ((c = getchar()) != EOF) {
        // Resize buffer if necessary
        if (bytes_read >= buffer_size - 1) {
            buffer_size *= 2;
            char* new_buffer = (char*) realloc(buffer, buffer_size);
            if (!new_buffer) {
                free(buffer);
                return -1;
            }
            buffer = new_buffer;
        }
        
        buffer[bytes_read++] = (char)c;
    }
    
    // Null-terminate the buffer
    buffer[bytes_read] = '\0';
    
    int result = load_from_wkt(buffer, geometry);
    free(buffer);
    
    return result;
}

/**
 * Load a geometry from stdin (WKB binary)
 */
static int load_from_stdin_wkb(sfcgal_geometry_t** geometry) {
    // Allocate initial buffer
    size_t buffer_size = BUFFER_SIZE;
    unsigned char* buffer = (unsigned char*) malloc(buffer_size);
    if (!buffer) {
        return -1;
    }
    
    size_t bytes_read = 0;
    int c;
    
    // Read input until EOF
    while ((c = getchar()) != EOF) {
        // Resize buffer if necessary
        if (bytes_read >= buffer_size) {
            buffer_size *= 2;
            unsigned char* new_buffer = (unsigned char*) realloc(buffer, buffer_size);
            if (!new_buffer) {
                free(buffer);
                return -1;
            }
            buffer = new_buffer;
        }
        
        buffer[bytes_read++] = (unsigned char)c;
    }
    
    // Use SFCGAL's function to read WKB directly
    *geometry = sfcgal_io_read_wkb((const char*)buffer, bytes_read);
    
    free(buffer);
    return (*geometry == NULL) ? -1 : 0;
}

/**
 * Load a geometry from a source
 */
int load_geometry(const char* source, sfcgal_geometry_t** geometry) {
    if (!source || !geometry) {
        return -1;
    }
    
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
    
    // Try to load from file
    // First try as WKT
    if (load_from_file_wkt(source, geometry) == 0) {
        return 0;
    }
    
    // Then try as WKB
    return load_from_file_wkb(source, geometry);
}

/**
 * Output a geometry as WKT
 */
static void output_as_wkt(const sfcgal_geometry_t* geometry, int precision) {
    char* wkt;
    size_t len;
    
    sfcgal_geometry_as_text_decim(geometry, precision, &wkt, &len);
    
    if (wkt) {
        printf("%s\n", wkt);
        sfcgal_free_buffer(wkt);
    }
}

/**
 * Output a geometry as WKB (hex)
 */
static void output_as_wkb(const sfcgal_geometry_t* geometry) {
    char* wkb;
    size_t len;
    
    sfcgal_geometry_as_hexwkb(geometry, &wkb, &len);
    
    if (wkb) {
        printf("%s\n", wkb);
        sfcgal_free_buffer(wkb);
    }
}

/**
 * Output a geometry as plain text
 */
static void output_as_txt(const sfcgal_geometry_t* geometry, int precision) {
    char* wkt;
    size_t len;
    
    sfcgal_geometry_as_text_decim(geometry, precision, &wkt, &len);
    
    if (wkt) {
        printf("%s\n", wkt);
        sfcgal_free_buffer(wkt);
    }
}

/**
 * Output a geometry as GeoJSON
 */
static void output_as_geojson(const sfcgal_geometry_t* geometry, int precision) {
    char* wkt;
    size_t len;
    
    // Get WKT representation first
    sfcgal_geometry_as_text_decim(geometry, precision, &wkt, &len);
    
    if (!wkt) {
        return;
    }
    
    // This is a very basic GeoJSON conversion and doesn't handle all cases properly
    // For a production implementation, a proper GeoJSON library should be used
    
    printf("{\n");
    printf("  \"type\": \"Feature\",\n");
    printf("  \"geometry\": {\n");
    
    if (strncmp(wkt, "POINT", 5) == 0) {
        printf("    \"type\": \"Point\",\n");
    } else if (strncmp(wkt, "LINESTRING", 10) == 0) {
        printf("    \"type\": \"LineString\",\n");
    } else if (strncmp(wkt, "POLYGON", 7) == 0) {
        printf("    \"type\": \"Polygon\",\n");
    } else if (strncmp(wkt, "MULTIPOINT", 10) == 0) {
        printf("    \"type\": \"MultiPoint\",\n");
    } else if (strncmp(wkt, "MULTILINESTRING", 15) == 0) {
        printf("    \"type\": \"MultiLineString\",\n");
    } else if (strncmp(wkt, "MULTIPOLYGON", 12) == 0) {
        printf("    \"type\": \"MultiPolygon\",\n");
    } else if (strncmp(wkt, "GEOMETRYCOLLECTION", 18) == 0) {
        printf("    \"type\": \"GeometryCollection\",\n");
    } else {
        printf("    \"type\": \"Unknown\",\n");
    }
    
    // For a proper implementation, we would need to parse the WKT and convert it to GeoJSON coordinates
    printf("    \"coordinates\": \"WKT conversion not fully implemented\"\n");
    printf("  },\n");
    printf("  \"properties\": {}\n");
    printf("}\n");
    
    sfcgal_free_buffer(wkt);
}

/**
 * Output a geometry in the specified format
 */
void output_geometry(const sfcgal_geometry_t* geometry, OutputFormat format, int precision) {
    if (!geometry) {
        return;
    }
    
    switch (format) {
        case FORMAT_WKT:
            output_as_wkt(geometry, precision);
            break;
        case FORMAT_WKB:
            output_as_wkb(geometry);
            break;
        case FORMAT_TXT:
            output_as_txt(geometry, precision);
            break;
        case FORMAT_GEOJSON:
            output_as_geojson(geometry, precision);
            break;
        default:
            fprintf(stderr, "Unknown output format\n");
            break;
    }
}
