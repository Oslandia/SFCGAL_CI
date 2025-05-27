/**
 * safe_string.h - Secure string and memory operations
 *
 * This module provides thread-safe and secure alternatives to standard C string
 * and memory functions, with comprehensive validation and error checking.
 */

#ifndef SAFE_STRING_H
#define SAFE_STRING_H

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32) && defined(_MSC_VER)
  #pragma warning(push)
  #pragma warning(disable : 4996) /* Disable deprecation warnings */
  #define strcasecmp _stricmp
  #define strncasecmp _strnicmp
#elif defined(__unix__) || defined(__APPLE__)
  #include <strings.h>
#endif

/* Maximum safe lengths for operations */
#define SAFE_MAX_STRING_LENGTH (1024 * 1024)     /* 1MB */
#define SAFE_MAX_BUFFER_SIZE (100 * 1024 * 1024) /* 100MB */
#define SAFE_MAX_PARAM_LENGTH 4096

/* ========================================================================== */
/* STRING MANIPULATION FUNCTIONS */
/* ========================================================================== */

/**
 * Thread-safe string tokenizer
 *
 * @param str String to tokenize (NULL for subsequent calls)
 * @param delim Delimiter characters
 * @param saveptr Pointer to store tokenizer state (must persist between calls)
 * @return Pointer to next token, or NULL if no more tokens
 */
char *
safe_strtok(char *str, const char *delim, char **saveptr);

/**
 * Safe string duplication with length validation
 *
 * @param str String to duplicate
 * @param max_len Maximum allowed length
 * @return Newly allocated string or NULL on error
 */
char *
safe_strdup(const char *str, size_t max_len);

/**
 * Safe string length calculation with maximum limit
 *
 * @param str String to measure
 * @param max_len Maximum length to check
 * @return String length or SIZE_MAX if exceeds max_len
 */
size_t
safe_strlen(const char *str, size_t max_len);

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
 * Check if a string is null or empty
 *
 * @param str String to check
 * @return true if string is NULL or empty, false otherwise
 */
bool
safe_string_is_empty(const char *str);

/* ========================================================================== */
/* STRING CONVERSION FUNCTIONS */
/* ========================================================================== */

/**
 * Safe string to long conversion with error checking
 *
 * @param str String to convert
 * @param endptr Pointer to store end of conversion (can be NULL)
 * @param base Numeric base (2-36)
 * @param result Pointer to store the result
 * @return true on success, false on error
 */
bool
safe_strtol(const char *str, char **endptr, int base, long *result);

/**
 * Safe string to double conversion with error checking
 *
 * @param str String to convert
 * @param endptr Pointer to store end of conversion (can be NULL)
 * @param result Pointer to store the result
 * @return true on success, false on error
 */
bool
safe_strtod(const char *str, char **endptr, double *result);

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

/* ========================================================================== */
/* MEMORY MANAGEMENT FUNCTIONS */
/* ========================================================================== */

/**
 * Safe memory allocation with size validation
 *
 * @param size Size to allocate
 * @param max_size Maximum allowed size
 * @return Allocated memory or NULL on error
 */
void *
safe_malloc(size_t size, size_t max_size);

/**
 * Safe memory reallocation with size validation
 *
 * @param ptr Existing pointer (can be NULL)
 * @param size New size
 * @param max_size Maximum allowed size
 * @return Reallocated memory or NULL on error
 */
void *
safe_realloc(void *ptr, size_t size, size_t max_size);

/**
 * Safe allocation with zero initialization
 *
 * @param size Size in bytes to allocate
 * @param max_size Maximum allowed size
 * @return Pointer to allocated memory or NULL on error
 */
void *
safe_calloc(size_t size, size_t max_size);

/**
 * Secure memory zeroing (prevents compiler optimization)
 *
 * @param ptr Pointer to memory to zero
 * @param size Number of bytes to zero
 */
void
safe_secure_zero(void *ptr, size_t size);

/* ========================================================================== */
/* MATHEMATICAL UTILITY FUNCTIONS */
/* ========================================================================== */

/**
 * Cross-platform isfinite implementation
 * Handles the differences between MinGW, MSVC, and POSIX systems
 */
static inline int
safe_isfinite_double(double x)
{
#if defined(_WIN32) && defined(__MINGW32__)
  /* MinGW specific implementation to avoid float conversion warning */
  return !(_isnan(x) || !_finite(x));
#elif defined(_WIN32) && defined(_MSC_VER)
  /* MSVC implementation */
  return _finite(x) && !_isnan(x);
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
  /* C99 and later */
  return isfinite(x);
#else
  /* Fallback for older compilers */
  return (x == x) && (x != HUGE_VAL) && (x != -HUGE_VAL);
#endif
}

/**
 * Cross-platform isnan implementation
 */
static inline int
safe_isnan_double(double x)
{
#if defined(_WIN32) && (defined(__MINGW32__) || defined(_MSC_VER))
  return _isnan(x);
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
  return isnan(x);
#else
  /* Fallback: NaN is the only value that is not equal to itself */
  return (x != x);
#endif
}

/* ========================================================================== */
/* CONVENIENCE MACROS */
/* ========================================================================== */

/**
 * Safe string duplication with default maximum length
 */
#define safe_strdup_default(str) safe_strdup(str, SAFE_MAX_STRING_LENGTH)

/**
 * Safe memory allocation with default maximum size
 */
#define safe_malloc_default(size) safe_malloc(size, SAFE_MAX_BUFFER_SIZE)

/**
 * Safe calloc with default maximum size
 */
#define safe_calloc_default(size) safe_calloc(size, SAFE_MAX_BUFFER_SIZE)

#ifdef __cplusplus
}
#endif

#endif /* SAFE_STRING_H */
