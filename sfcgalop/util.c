/**
 * util.c - Implementation of SFCGAL-specific utility functions
 */

#include "util.h"
#include "safe_string.h"
#include "sfcgal_c.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

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

/* ========================================================================== */
/* SFCGAL ERROR HANDLERS */
/* ========================================================================== */

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
  if (safe_strlen(format, MAX_STRING_LENGTH) == SIZE_MAX) {
    fprintf(stderr, "SFCGAL WARNING: Format string too long\n");
    return 0;
  }

  va_list args;

  fprintf(stderr, "SFCGAL WARNING: ");

  va_start(args, format);
  vfprintf(stderr, format, args); // flawfinder: ignore - format validated above
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
  if (safe_strlen(format, MAX_STRING_LENGTH) == SIZE_MAX) {
    fprintf(stderr, "SFCGAL ERROR: Format string too long\n");
    return 0;
  }

  va_list args;

  fprintf(stderr, "SFCGAL ERROR: ");

  va_start(args, format);
  vfprintf(stderr, format, args); // flawfinder: ignore - format validated above
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

#if defined(_WIN32) && defined(_MSC_VER)
  #pragma warning(pop)
#endif
