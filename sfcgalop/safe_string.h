/**
 * safe_string.h - Thread-safe string manipulation functions
 *
 * This module provides thread-safe alternatives to standard C string functions
 * that are not thread-safe (like strtok) and includes additional safety checks.
 */

#ifndef SAFE_STRING_H
#define SAFE_STRING_H

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Maximum safe lengths for string operations */
#define MAX_SAFE_STRING_LENGTH (1024 * 1024)     /* 1MB */
#define MAX_SAFE_BUFFER_SIZE (100 * 1024 * 1024) /* 100MB */
#define MAX_PARAM_LENGTH 4096

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
safe_strdup_checked(const char *str, size_t max_len);

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
 * Secure memory zeroing (prevents compiler optimization)
 *
 * @param ptr Pointer to memory to zero
 * @param size Number of bytes to zero
 */
void
secure_zero(void *ptr, size_t size);

/**
 * Safe memory allocation with size validation
 *
 * @param size Size to allocate
 * @param max_size Maximum allowed size
 * @return Allocated memory or NULL on error
 */
void *
safe_malloc_checked(size_t size, size_t max_size);

/**
 * Safe memory reallocation with size validation
 *
 * @param ptr Existing pointer (can be NULL)
 * @param size New size
 * @param max_size Maximum allowed size
 * @return Reallocated memory or NULL on error
 */
void *
safe_realloc_checked(void *ptr, size_t size, size_t max_size);

#ifdef __cplusplus
}
#endif

#endif /* SAFE_STRING_H */
