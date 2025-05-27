/**
 * util.h - SFCGAL-specific utility functions
 *
 * This module contains utility functions and configurations specifically
 * for the SFCGAL operations project.
 */

#ifndef UTIL_H
#define UTIL_H

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ========================================================================== */
/* PROJECT CONSTANTS */
/* ========================================================================== */

/* Application-level limits (different from safe_string internal limits) */
#define MAX_STRING_LENGTH (1024 * 1024)     /* 1MB */
#define MAX_BUFFER_SIZE (100 * 1024 * 1024) /* 100MB */

/* ========================================================================== */
/* SFCGAL-SPECIFIC FUNCTIONS */
/* ========================================================================== */

/**
 * Set up error handlers for SFCGAL library
 * Configures the SFCGAL library to use custom warning and error handlers
 */
void
setup_error_handlers(void);

#ifdef __cplusplus
}
#endif

#if defined(_WIN32) && defined(_MSC_VER)
  #pragma warning(pop)
#endif

#endif /* UTIL_H */
