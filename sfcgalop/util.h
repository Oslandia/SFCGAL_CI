/**
 * util.h - Utility function declarations for SFCGAL operations
 */

#ifndef UTIL_H
#define UTIL_H

#include <stdbool.h>
#include <stddef.h>

#if defined(_WIN32) && defined(_MSC_VER)
  #pragma warning(push)
  #pragma warning(disable : 4996) /* Disable deprecation warnings */
  #define strcasecmp _stricmp
  #define strncasecmp _strnicmp
#elif defined(__unix__) || defined(__APPLE__)
  #include <strings.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Security limits */
#define MAX_STRING_LENGTH (1024 * 1024)     /* 1MB */
#define MAX_BUFFER_SIZE (100 * 1024 * 1024) /* 100MB */

/**
 * Set up error handlers for SFCGAL
 */
void
setup_error_handlers(void);

/**
 * Check if a string is null or empty
 *
 * @param str String to check
 * @return true if string is NULL or empty, false otherwise
 */
bool
is_null_or_empty(const char *str);

/**
 * Safe string duplication (legacy interface)
 *
 * @param str String to duplicate
 * @return Newly allocated copy of the string or NULL on allocation failure
 */
char *
safe_strdup(const char *str);

/**
 * Safe allocation with zero initialization (legacy interface)
 *
 * @param size Size in bytes to allocate
 * @return Pointer to allocated memory or NULL on allocation failure
 */
void *
safe_calloc(size_t size);

/**
 * Enhanced safe string duplication with length validation
 *
 * @param str String to duplicate
 * @param max_len Maximum allowed length
 * @return Newly allocated copy of the string or NULL on error
 */
char *
safe_strdup_enhanced(const char *str, size_t max_len);

/**
 * Enhanced safe allocation with size validation
 *
 * @param size Size in bytes to allocate
 * @param max_size Maximum allowed size
 * @return Pointer to allocated memory or NULL on error
 */
void *
safe_calloc_enhanced(size_t size, size_t max_size);

/**
 * Safe string length calculation
 *
 * @param str String to measure
 * @param max_len Maximum length to check
 * @return String length or SIZE_MAX if exceeds max_len
 */
size_t
safe_string_length(const char *str, size_t max_len);

/**
 * Safe memory reallocation
 *
 * @param ptr Existing pointer (can be NULL)
 * @param new_size New size in bytes
 * @param max_size Maximum allowed size
 * @return Reallocated memory or NULL on error
 */
void *
safe_realloc_enhanced(void *ptr, size_t new_size, size_t max_size);

/**
 * Secure memory clearing
 *
 * @param ptr Pointer to memory to clear
 * @param size Number of bytes to clear
 */
void
secure_memzero(void *ptr, size_t size);

/**
 * Safe string concatenation
 *
 * @param dest Destination buffer
 * @param dest_size Size of destination buffer
 * @param src Source string
 * @return true on success, false on error
 */
bool
safe_strcat(char *dest, size_t dest_size, const char *src);

/**
 * Safe string copy
 *
 * @param dest Destination buffer
 * @param dest_size Size of destination buffer
 * @param src Source string
 * @return true on success, false on error
 */
bool
safe_strcpy(char *dest, size_t dest_size, const char *src);

/**
 * Safe integer parsing with overflow checking
 *
 * @param str String to parse
 * @param result Pointer to store result
 * @param base Numeric base (2-36)
 * @return true on success, false on error
 */
bool
safe_parse_int(const char *str, int *result, int base);

/**
 * Safe double parsing with overflow checking
 *
 * @param str String to parse
 * @param result Pointer to store result
 * @return true on success, false on error
 */
bool
safe_parse_double(const char *str, double *result);

#ifdef __cplusplus
}
#endif

#endif /* UTIL_H */
