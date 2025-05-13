/**
 * util.h - Utility function declarations for SFCGAL operations
 */

#ifndef UTIL_H
#define UTIL_H

#include <stdbool.h>
#include <stddef.h>

/**
 * Set up error handlers for SFCGAL
 */
void setup_error_handlers(void);

/**
 * Check if a string is null or empty
 * 
 * @param str String to check
 * @return true if string is NULL or empty, false otherwise
 */
bool is_null_or_empty(const char* str);

/**
 * Safe string duplication
 * 
 * @param str String to duplicate
 * @return Newly allocated copy of the string or NULL on allocation failure
 */
char* safe_strdup(const char* str);

/**
 * Safe allocation with zero initialization
 * 
 * @param size Size in bytes to allocate
 * @return Pointer to allocated memory or NULL on allocation failure
 */
void* safe_calloc(size_t size);

#endif /* UTIL_H */
