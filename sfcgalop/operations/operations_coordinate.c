/**
 * operations_coordinate.c - Implementation of coordinate manipulation
 * operations
 */

#include "operations_coordinate.h"
#include "../util.h"
#include "operations_common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * In-place operations defined using macros
 */
DEFINE_INPLACE_OP(drop_z, sfcgal_geometry_drop_z)
DEFINE_INPLACE_OP(drop_m, sfcgal_geometry_drop_m)

/**
 * Operations that return a new geometry
 */
DEFINE_GEOMETRY_OP_SINGLE(force_lhr, sfcgal_geometry_force_lhr)
DEFINE_GEOMETRY_OP_SINGLE(force_rhr, sfcgal_geometry_force_rhr)

/**
 * Add Z coordinate to geometry with given default value
 */
OperationResult
op_force_z(const char *op_arg, const sfcgal_geometry_t *geom_a,
           const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  // Ignore geom_b
  (void)geom_b;

  // Parse defaultZ parameter
  double defaultZ = 0.0;
  if (!parse_double_parameter(op_arg, "defaultZ", 0, 0.0, &defaultZ)) {
    result.error         = true;
    result.error_message = "Invalid default Z value";
    return result;
  }

  // Clone geometry first to avoid modifying input
  result.geometry_result = sfcgal_geometry_clone(geom_a);
  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Failed to clone geometry";
    return result;
  }

  int changed = sfcgal_geometry_force_z(result.geometry_result, defaultZ);

  if (changed < 0) {
    result.error         = true;
    result.error_message = "Error forcing Z coordinate";
    sfcgal_geometry_delete(result.geometry_result);
    result.geometry_result = NULL;
  }

  return result;
}

/**
 * Add M coordinate to geometry with given default value
 */
OperationResult
op_force_m(const char *op_arg, const sfcgal_geometry_t *geom_a,
           const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  // Ignore geom_b
  (void)geom_b;

  // Parse defaultM parameter
  double defaultM = 0.0;
  if (!parse_double_parameter(op_arg, "defaultM", 0, 0.0, &defaultM)) {
    result.error         = true;
    result.error_message = "Invalid default M value";
    return result;
  }

  // Clone geometry first to avoid modifying input
  result.geometry_result = sfcgal_geometry_clone(geom_a);
  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Failed to clone geometry";
    return result;
  }

  int changed = sfcgal_geometry_force_m(result.geometry_result, defaultM);

  if (changed < 0) {
    result.error         = true;
    result.error_message = "Error forcing M coordinate";
    sfcgal_geometry_delete(result.geometry_result);
    result.geometry_result = NULL;
  }

  return result;
}
