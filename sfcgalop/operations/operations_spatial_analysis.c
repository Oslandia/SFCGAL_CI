/**
 * operations_spatial_analysis.c - Implementation of spatial analysis operations
 */

#include "operations_spatial_analysis.h"
#include "../util.h"
#include "operations_common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Basic boolean operations defined using macros
 */
DEFINE_BOOLEAN_OP_BINARY(intersects, sfcgal_geometry_intersects)
DEFINE_BOOLEAN_OP_BINARY(intersects_3d, sfcgal_geometry_intersects_3d)
DEFINE_BOOLEAN_OP_BINARY(covers, sfcgal_geometry_covers)
DEFINE_BOOLEAN_OP_BINARY(covers_3d, sfcgal_geometry_covers_3d)
DEFINE_BOOLEAN_OP_BINARY(equals, sfcgal_geometry_is_equals)

/**
 * Almost equals with tolerance parameter
 */
OperationResult
op_almost_equals(const char *op_arg, const sfcgal_geometry_t *geom_a,
                 const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_BOOLEAN, .error = false};

  // Parse tolerance parameter
  double tolerance = 0.0;
  if (!parse_double_parameter(op_arg, "tolerance", 0, 0.0, &tolerance)) {
    result.error         = true;
    result.error_message = "Invalid tolerance value";
    return result;
  }

  int ret = sfcgal_geometry_is_almost_equals(geom_a, geom_b, tolerance);

  if (ret < 0) {
    result.error         = true;
    result.error_message = "Error checking almost equals";
    return result;
  }

  result.boolean_result = (ret == 1);
  return result;
}
