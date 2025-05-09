/**
 * operations_skeleton.c - Implementation of skeleton operations
 */

#include "operations_skeleton.h"
#include "../util.h"
#include "operations_common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Skeleton operations using macros for simple cases
 */
DEFINE_GEOMETRY_OP_SINGLE(straight_skeleton, sfcgal_geometry_straight_skeleton)
DEFINE_GEOMETRY_OP_SINGLE(straight_skeleton_distance_in_m,
                          sfcgal_geometry_straight_skeleton_distance_in_m)
DEFINE_GEOMETRY_OP_SINGLE(approximate_medial_axis,
                          sfcgal_geometry_approximate_medial_axis)

/**
 * Extrude straight skeleton operation with height parameter
 */
OperationResult
op_extrude_straight_skeleton(const char              *op_arg,
                             const sfcgal_geometry_t *geom_a,
                             const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse height parameter
  double height = 1.0;
  if (!parse_double_parameter(op_arg, "height", 0, 1.0, &height)) {
    result.error         = true;
    result.error_message = "Invalid height value";
    return result;
  }

  result.geometry_result =
      sfcgal_geometry_extrude_straight_skeleton(geom_a, height);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error extruding straight skeleton";
  }

  return result;
}

/**
 * Extrude polygon straight skeleton operation with building height and roof
 * height parameters
 */
OperationResult
op_extrude_polygon_straight_skeleton(const char              *op_arg,
                                     const sfcgal_geometry_t *geom_a,
                                     const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse parameters
  double building_height = 1.0, roof_height = 1.0;
  if (!parse_double_parameter(op_arg, "building_height", 0, 1.0,
                              &building_height) ||
      !parse_double_parameter(op_arg, "roof_height", 1, 1.0, &roof_height)) {
    result.error         = true;
    result.error_message = "Invalid parameter value(s)";
    return result;
  }

  result.geometry_result = sfcgal_geometry_extrude_polygon_straight_skeleton(
      geom_a, building_height, roof_height);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error extruding polygon straight skeleton";
  }

  return result;
}

/**
 * Straight skeleton partition operation with auto orientation parameter
 */
OperationResult
op_straight_skeleton_partition(const char              *op_arg,
                               const sfcgal_geometry_t *geom_a,
                               const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse auto_orientation parameter
  bool auto_orientation = true;
  if (!parse_bool_parameter(op_arg, "auto_orientation", 0, true,
                            &auto_orientation)) {
    result.error         = true;
    result.error_message = "Invalid auto_orientation value";
    return result;
  }

  result.geometry_result =
      sfcgal_geometry_straight_skeleton_partition(geom_a, auto_orientation);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error calculating straight skeleton partition";
  }

  return result;
}
