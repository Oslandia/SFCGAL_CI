/**
 * util.c - Implementation of utility functions for SFCGAL operations
 */

#include "util.h"
#include "sfcgal_c.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

/**
 * Custom warning handler for SFCGAL
 * 
 * @param format Printf-style format string
 * @param ... Variable arguments for format string
 * @return Always returns 0
 */
static int warning_handler(const char* format, ...) {
    if (!format) {
        return 0;
    }
    
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
 * 
 * @param format Printf-style format string
 * @param ... Variable arguments for format string
 * @return Always returns 0
 */
static int error_handler(const char* format, ...) {
    if (!format) {
        return 0;
    }
    
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

/**
 * Check if a string is null or empty
 */
bool is_null_or_empty(const char* str) {
    return str == NULL || str[0] == '\0';
}

/**
 * Safe string duplication
 */
char* safe_strdup(const char* str) {
    if (!str) {
        return NULL;
    }
    
    return strdup(str);
}

/**
 * Safe allocation with zero initialization
 */
void* safe_calloc(size_t size) {
    if (size == 0) {
        return NULL;
    }
    
    return calloc(1, size);
}
