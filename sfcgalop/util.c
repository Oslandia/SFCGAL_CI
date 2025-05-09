/**
 * util.c - Implementation of utility functions for SFCGAL operations
 */

#include "util.h"
#include "safe_string.h"
#include "sfcgal_c.h"
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(_WIN32) && defined(_MSC_VER)
  #include <windows.h>
  #pragma warning(push)
  #pragma warning(disable : 4996) /* Disable deprecation warnings */
  #define strcasecmp _stricmp
  #define strncasecmp _strnicmp
#endif

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wformat-nonliteral"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wformat-nonliteral"
#endif
/**
 * Custom warning handler for SFCGAL
 *
 * @param format Printf-style format string
 * @param ... Variable arguments for format string
 * @return Always returns 0
 */
static int
warning_handler(const char *format, ...)
{
  if (!format) {
    return 0;
  }

  /* Validate format string length to prevent DoS */
  if (safe_string_length(format, MAX_STRING_LENGTH) == SIZE_MAX) {
    fprintf(stderr, "SFCGAL WARNING: Format string too long\n");
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
static int
error_handler(const char *format, ...)
{
  if (!format) {
    return 0;
  }

  /* Validate format string length to prevent DoS */
  if (safe_string_length(format, MAX_STRING_LENGTH) == SIZE_MAX) {
    fprintf(stderr, "SFCGAL ERROR: Format string too long\n");
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
#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

/**
 * Set up error handlers for SFCGAL
 */
void
setup_error_handlers(void)
{
  sfcgal_set_error_handlers(warning_handler, error_handler);
}

/**
 * Check if a string is null or empty
 */
bool
is_null_or_empty(const char *str)
{
  return str == NULL || str[0] == '\0';
}

/**
 * Safe string duplication (legacy interface)
 */
char *
safe_strdup(const char *str)
{
  return safe_strdup_enhanced(str, MAX_STRING_LENGTH);
}

/**
 * Safe allocation with zero initialization (legacy interface)
 */
void *
safe_calloc(size_t size)
{
  return safe_calloc_enhanced(size, MAX_BUFFER_SIZE);
}

/**
 * Enhanced safe string duplication with length validation
 */
char *
safe_strdup_enhanced(const char *str, size_t max_len)
{
  if (!str) {
    return NULL;
  }

  return safe_strdup_checked(str, max_len);
}

/**
 * Enhanced safe allocation with size validation
 */
void *
safe_calloc_enhanced(size_t size, size_t max_size)
{
  if (size == 0 || size > max_size) {
    errno = EINVAL;
    return NULL;
  }

  return safe_malloc_checked(size, max_size);
}

/**
 * Safe string length calculation
 */
size_t
safe_string_length(const char *str, size_t max_len)
{
  return safe_strlen(str, max_len);
}

/**
 * Safe memory reallocation
 */
void *
safe_realloc_enhanced(void *ptr, size_t new_size, size_t max_size)
{
  return safe_realloc_checked(ptr, new_size, max_size);
}

/**
 * Secure memory clearing
 */
void
secure_memzero(void *ptr, size_t size)
{
  secure_zero(ptr, size);
}

/**
 * Safe string concatenation
 */
bool
safe_strcat(char *dest, size_t dest_size, const char *src)
{
  if (!dest || !src || dest_size == 0) {
    return false;
  }

  size_t dest_len = safe_strlen(dest, dest_size);
  size_t src_len  = safe_strlen(src, MAX_STRING_LENGTH);

  if (dest_len == SIZE_MAX || src_len == SIZE_MAX) {
    return false;
  }

  /* Check if concatenation would overflow */
  if (dest_len + src_len >= dest_size) {
    errno = ERANGE;
    return false;
  }

#if defined(_WIN32) && defined(_MSC_VER)
  errno_t err = strcat_s(dest, dest_size, src);
  return (err == 0);
#elif defined(__STDC_LIB_EXT1__)
  errno_t err = strcat_s(dest, dest_size, src);
  return (err == 0);
#else
  strcat(dest, src);
  return true;
#endif
}

/**
 * Safe string copy
 */
bool
safe_strcpy(char *dest, size_t dest_size, const char *src)
{
  if (!dest || !src || dest_size == 0) {
    return false;
  }

  size_t src_len = safe_strlen(src, MAX_STRING_LENGTH);
  if (src_len == SIZE_MAX || src_len >= dest_size) {
    errno = ERANGE;
    return false;
  }

#if defined(_WIN32) && defined(_MSC_VER)
  errno_t err = strcpy_s(dest, dest_size, src);
  return (err == 0);
#elif defined(__STDC_LIB_EXT1__)
  errno_t err = strcpy_s(dest, dest_size, src);
  return (err == 0);
#else
  strcpy(dest, src);
  return true;
#endif
}

/**
 * Safe integer parsing with overflow checking
 */
bool
safe_parse_int(const char *str, int *result, int base)
{
  if (!str || !result) {
    return false;
  }

  long long_result;
  if (!safe_strtol(str, NULL, base, &long_result)) {
    return false;
  }

  /* Check for integer overflow */
  if (long_result > INT_MAX || long_result < INT_MIN) {
    errno = ERANGE;
    return false;
  }

  *result = (int)long_result;
  return true;
}

/**
 * Safe double parsing with overflow checking
 */
bool
safe_parse_double(const char *str, double *result)
{
  if (!str || !result) {
    return false;
  }

  return safe_strtod(str, NULL, result);
}
#if defined(_WIN32) && defined(_MSC_VER)
  #pragma warning(pop)
#endif
