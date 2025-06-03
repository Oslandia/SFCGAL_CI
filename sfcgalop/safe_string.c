/**
 * safe_string.c - Implementation of secure string and memory operations
 */

#include "safe_string.h"
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(_WIN32)
  #include <windows.h>
#endif

/* ========================================================================== */
/* STRING MANIPULATION FUNCTIONS */
/* ========================================================================== */

/**
 * Thread-safe string tokenizer
 * Cross-platform implementation using strtok_r (POSIX) or strtok_s (Windows)
 */
char *
safe_strtok(char *str, const char *delim, char **saveptr)
{
  if (!delim || !saveptr) {
    return NULL;
  }

#if defined(_WIN32) && defined(_MSC_VER)
  /* Use Microsoft's secure strtok_s */
  return strtok_s(str, delim, saveptr);
#elif defined(__STDC_LIB_EXT1__)
  /* Use C11 strtok_s if available */
  return strtok_s(str, delim, saveptr);
#else
  /* Use POSIX strtok_r */
  return strtok_r(str, delim, saveptr);
#endif
}

/**
 * Safe string duplication with length validation
 */
char *
safe_strdup(const char *str, size_t max_len)
{
  if (!str) {
    return NULL;
  }

  size_t len = safe_strlen(str, max_len);
  if (len == SIZE_MAX || len >= max_len) {
    errno = EINVAL;
    return NULL;
  }

  char *result = safe_malloc(len + 1, max_len + 1);
  if (!result) {
    return NULL;
  }

  memcpy(result, str, len); // flawfinder: ignore - length validated above
  result[len] = '\0';
  return result;
}

/**
 * Safe string length calculation with maximum limit
 */
size_t
safe_strlen(const char *str, size_t max_len)
{
  if (!str) {
    return 0;
  }

  size_t len = 0;
  while (len < max_len && str[len] != '\0') {
    len++;
  }

  return (len == max_len && str[len] != '\0') ? SIZE_MAX : len;
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
  size_t src_len  = safe_strlen(src, SAFE_MAX_STRING_LENGTH);

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
  strcat(dest, src); // flawfinder: ignore - bounds checked above
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

  size_t src_len = safe_strlen(src, SAFE_MAX_STRING_LENGTH);
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
  strcpy(dest, src); // flawfinder: ignore - bounds checked above
  return true;
#endif
}

/**
 * Check if a string is null or empty
 */
bool
safe_string_is_empty(const char *str)
{
  return str == NULL || str[0] == '\0';
}

/* ========================================================================== */
/* STRING CONVERSION FUNCTIONS */
/* ========================================================================== */

/**
 * Safe string to long conversion with error checking
 */
bool
safe_strtol(const char *str, char **endptr, int base, long *result)
{
  if (!str || !result) {
    errno = EINVAL;
    return false;
  }

  if (base < 2 || base > 36) {
    errno = EINVAL;
    return false;
  }

  /* Check string length to prevent DoS */
  if (safe_strlen(str, SAFE_MAX_PARAM_LENGTH) == SIZE_MAX) {
    errno = EINVAL;
    return false;
  }

  errno              = 0;
  char *local_endptr = NULL;
  long  value        = strtol(str, &local_endptr, base);

  /* Check for conversion errors */
  if (errno == ERANGE || (errno != 0 && value == 0)) {
    return false;
  }

  /* Check if no conversion was performed */
  if (local_endptr == str) {
    errno = EINVAL;
    return false;
  }

  if (endptr) {
    *endptr = local_endptr;
  }

  *result = value;
  return true;
}

/**
 * Safe string to double conversion with error checking
 */
bool
safe_strtod(const char *str, char **endptr, double *result)
{
  if (!str || !result) {
    errno = EINVAL;
    return false;
  }

  /* Check string length to prevent DoS */
  if (safe_strlen(str, SAFE_MAX_PARAM_LENGTH) == SIZE_MAX) {
    errno = EINVAL;
    return false;
  }

  errno               = 0;
  char  *local_endptr = NULL;
  double value        = strtod(str, &local_endptr);

  /* Check for conversion errors */
  if (errno == ERANGE) {
    return false;
  }

  /* Check for infinite or NaN results */
  if (!safe_isfinite_double(value)) {
    errno = ERANGE;
    return false;
  }

  /* Check if no conversion was performed */
  if (local_endptr == str) {
    errno = EINVAL;
    return false;
  }

  if (endptr) {
    *endptr = local_endptr;
  }

  *result = value;
  return true;
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

/* ========================================================================== */
/* MEMORY MANAGEMENT FUNCTIONS */
/* ========================================================================== */

/**
 * Safe memory allocation with size validation
 */
void *
safe_malloc(size_t size, size_t max_size)
{
  if (size == 0 || size > max_size) {
    errno = EINVAL;
    return NULL;
  }

  /* Check for integer overflow */
  if (size > SIZE_MAX / 2) {
    errno = ENOMEM;
    return NULL;
  }

  void *ptr = malloc(size);
  if (!ptr) {
    errno = ENOMEM;
    return NULL;
  }

  /* Zero initialize for security */
  memset(ptr, 0, size);
  return ptr;
}

/**
 * Safe memory reallocation with size validation
 */
void *
safe_realloc(void *ptr, size_t size, size_t max_size)
{
  if (size > max_size) {
    errno = EINVAL;
    return NULL;
  }

  if (size == 0) {
    free(ptr);
    return NULL;
  }

  /* Check for integer overflow */
  if (size > SIZE_MAX / 2) {
    errno = ENOMEM;
    return NULL;
  }

  void *new_ptr = realloc(ptr, size);
  if (!new_ptr) {
    errno = ENOMEM;
    return NULL;
  }

  return new_ptr;
}

/**
 * Safe allocation with zero initialization
 */
void *
safe_calloc(size_t size, size_t max_size)
{
  if (size == 0 || size > max_size) {
    errno = EINVAL;
    return NULL;
  }

  return safe_malloc(size, max_size);
}

/**
 * Secure memory zeroing (prevents compiler optimization)
 */
void
safe_secure_zero(void *ptr, size_t size)
{
  if (!ptr || size == 0) {
    return;
  }

#if defined(_WIN32) && defined(_MSC_VER)
  SecureZeroMemory(ptr, size);
#elif defined(__STDC_LIB_EXT1__)
  memset_s(ptr, size, 0, size);
#else
  /* Prevent compiler optimization with volatile */
  volatile unsigned char *p = (volatile unsigned char *)ptr;
  for (size_t i = 0; i < size; i++) {
    p[i] = 0;
  }
#endif
}
