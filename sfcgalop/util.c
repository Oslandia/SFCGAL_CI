/**
 * util.c - Implementation of utility functions
 */

#include "util.h"
#include "sfcgal_c.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

/**
 * Custom warning handler for SFCGAL
 */
static int warning_handler(const char* format, ...) {
    va_list args;
    
    fprintf(stderr, "SFCGAL WARNING: ");
    
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    
    fprintf(stderr, "\n");
    
    return 0;
}

/**
 * Custom error handler for SFCGAL
 */
static int error_handler(const char* format, ...) {
    va_list args;
    
    fprintf(stderr, "SFCGAL ERROR: ");
    
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    
    fprintf(stderr, "\n");
    
    return 0;
}

/**
 * Set up error handlers for SFCGAL
 */
void setup_error_handlers(void) {
    sfcgal_set_error_handlers(warning_handler, error_handler);
}
