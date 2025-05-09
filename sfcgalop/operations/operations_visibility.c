/**
 * operations_visibility.c - Implementation of visibility operations
 */

#include "operations_visibility.h"
#include "../util.h"
#include "operations_common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Visibility from point operation
 */
DEFINE_GEOMETRY_OP_BINARY(visibility_point, sfcgal_geometry_visibility_point)

/**
 * Visibility segment operation with two point geometries
 */
OperationResult
op_visibility_segment(const char *op_arg, const sfcgal_geometry_t *geom_a,
                      const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  // We need a third geometry (point B)
  sfcgal_geometry_t *point_b = NULL;

  // Parse the third geometry from op_arg
  if (op_arg && *op_arg) {
    const char *wkt_value = NULL;
    bool        found     = false;

    wkt_value = find_parameter(op_arg, "point_b", &found);

    if (found && wkt_value) {
      point_b = sfcgal_io_read_wkt(wkt_value, strlen(wkt_value));
      free((void *)(uintptr_t)wkt_value);
    } else {
      // Try to use the full op_arg as a WKT string
      point_b = sfcgal_io_read_wkt(op_arg, strlen(op_arg));
    }

    if (!point_b) {
      result.error         = true;
      result.error_message = "Failed to parse third geometry (point B)";
      return result;
    }
  } else {
    result.error = true;
    result.error_message =
        "Visibility segment requires a third geometry (point B)";
    return result;
  }

  result.geometry_result =
      sfcgal_geometry_visibility_segment(geom_a, geom_b, point_b);

  // Free the point we created
  sfcgal_geometry_delete(point_b);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error calculating visibility segment";
  }

  return result;
}
