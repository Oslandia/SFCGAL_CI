/**
 * operations_validity.c - Implementation of geometry validity operations
 */

#include "operations_validity.h"
#include "../safe_string.h"
#include "operations_common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Check if a geometry is valid
 */
DEFINE_BOOLEAN_OP_SINGLE(is_valid, sfcgal_geometry_is_valid)

/**
 * Check if a polygon is planar
 */
DEFINE_BOOLEAN_OP_SINGLE(is_planar, sfcgal_geometry_is_planar)

/**
 * Get detailed information about geometry validity
 */
OperationResult
op_is_validity_detail(const char *op_arg, const sfcgal_geometry_t *geom_a,
                      const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_TEXT, .error = false};

  // Ignore op_arg and geom_b
  (void)op_arg;
  (void)geom_b;

  char              *reason   = NULL;
  sfcgal_geometry_t *location = NULL;

  int valid = sfcgal_geometry_is_valid_detail(geom_a, &reason, &location);

  if (valid < 0) {
    result.error         = true;
    result.error_message = "Error checking validity detail";
    if (reason)
      free(reason);
    if (location)
      sfcgal_geometry_delete(location);
    return result;
  }

  if (valid == 1) {
    result.text_result =
        safe_strdup("Geometry is valid", SAFE_MAX_STRING_LENGTH);
  } else {
    if (reason != NULL) {

      size_t reason_len = safe_strlen(reason, SAFE_MAX_STRING_LENGTH);
      if (reason_len != SIZE_MAX) {
        size_t needed =
            reason_len + 50; // flawfinder: ignore - safe_strlen used
        char *full_message = malloc(needed);
        if (full_message) {
          snprintf(full_message, needed, "Geometry is invalid. Reason: %s",
                   reason);
          result.text_result =
              safe_strdup(full_message, SAFE_MAX_STRING_LENGTH);
          free(full_message);
        }
      }
      free(reason);
    } else {
      result.text_result = safe_strdup(
          "Geometry is invalid. No details available.", SAFE_MAX_STRING_LENGTH);
    }
  }

  if (location)
    sfcgal_geometry_delete(location);

  return result;
}
