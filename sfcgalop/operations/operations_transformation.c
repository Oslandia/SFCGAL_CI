/**
 * operations_transformation.c - Implementation of geometry transformation
 * operations
 */

#include "operations_transformation.h"
#include "../util.h"
#include "operations_common.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Rotate operation with angle parameter
 */
OperationResult
op_rotate(const char *op_arg, const sfcgal_geometry_t *geom_a,
          const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse angle parameter
  double angle = 0.0;
  if (!parse_double_parameter(op_arg, "angle", 0, 0.0, &angle)) {
    result.error         = true;
    result.error_message = "Invalid angle value";
    return result;
  }

  result.geometry_result = sfcgal_geometry_rotate(geom_a, angle);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error rotating geometry";
  }

  return result;
}

/**
 * Rotate 2D operation with angle, cx, cy parameters
 */
OperationResult
op_rotate_2d(const char *op_arg, const sfcgal_geometry_t *geom_a,
             const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse parameters
  double angle = 0.0, cx = 0.0, cy = 0.0;

  if (!parse_double_parameter(op_arg, "angle", 0, 0.0, &angle) ||
      !parse_double_parameter(op_arg, "cx", 1, 0.0, &cx) ||
      !parse_double_parameter(op_arg, "cy", 2, 0.0, &cy)) {
    result.error         = true;
    result.error_message = "Invalid parameter value(s)";
    return result;
  }

  result.geometry_result = sfcgal_geometry_rotate_2d(geom_a, angle, cx, cy);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error rotating geometry in 2D";
  }

  return result;
}

/**
 * Rotate 3D operation with angle, ax, ay, az parameters
 */
OperationResult
op_rotate_3d(const char *op_arg, const sfcgal_geometry_t *geom_a,
             const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse parameters
  double angle = 0.0, ax = 0.0, ay = 0.0, az = 1.0;

  if (!parse_double_parameter(op_arg, "angle", 0, 0.0, &angle) ||
      !parse_double_parameter(op_arg, "ax", 1, 0.0, &ax) ||
      !parse_double_parameter(op_arg, "ay", 2, 0.0, &ay) ||
      !parse_double_parameter(op_arg, "az", 3, 1.0, &az)) {
    result.error         = true;
    result.error_message = "Invalid parameter value(s)";
    return result;
  }

  result.geometry_result = sfcgal_geometry_rotate_3d(geom_a, angle, ax, ay, az);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error rotating geometry in 3D";
  }

  return result;
}

/**
 * Rotate 3D around center operation with angle, ax, ay, az, cx, cy, cz
 * parameters
 */
OperationResult
op_rotate_3d_around_center(const char *op_arg, const sfcgal_geometry_t *geom_a,
                           const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse parameters
  double angle = 0.0, ax = 0.0, ay = 0.0, az = 1.0, cx = 0.0, cy = 0.0,
         cz = 0.0;

  if (!parse_double_parameter(op_arg, "angle", 0, 0.0, &angle) ||
      !parse_double_parameter(op_arg, "ax", 1, 0.0, &ax) ||
      !parse_double_parameter(op_arg, "ay", 2, 0.0, &ay) ||
      !parse_double_parameter(op_arg, "az", 3, 1.0, &az) ||
      !parse_double_parameter(op_arg, "cx", 4, 0.0, &cx) ||
      !parse_double_parameter(op_arg, "cy", 5, 0.0, &cy) ||
      !parse_double_parameter(op_arg, "cz", 6, 0.0, &cz)) {
    result.error         = true;
    result.error_message = "Invalid parameter value(s)";
    return result;
  }

  result.geometry_result = sfcgal_geometry_rotate_3d_around_center(
      geom_a, angle, ax, ay, az, cx, cy, cz);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error rotating geometry in 3D around center";
  }

  return result;
}

/**
 * Rotate X operation with angle parameter
 */
OperationResult
op_rotate_x(const char *op_arg, const sfcgal_geometry_t *geom_a,
            const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse angle parameter
  double angle = 0.0;
  if (!parse_double_parameter(op_arg, "angle", 0, 0.0, &angle)) {
    result.error         = true;
    result.error_message = "Invalid angle value";
    return result;
  }

  result.geometry_result = sfcgal_geometry_rotate_x(geom_a, angle);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error rotating geometry around X axis";
  }

  return result;
}

/**
 * Rotate Y operation with angle parameter
 */
OperationResult
op_rotate_y(const char *op_arg, const sfcgal_geometry_t *geom_a,
            const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse angle parameter
  double angle = 0.0;
  if (!parse_double_parameter(op_arg, "angle", 0, 0.0, &angle)) {
    result.error         = true;
    result.error_message = "Invalid angle value";
    return result;
  }

  result.geometry_result = sfcgal_geometry_rotate_y(geom_a, angle);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error rotating geometry around Y axis";
  }

  return result;
}

/**
 * Rotate Z operation with angle parameter
 */
OperationResult
op_rotate_z(const char *op_arg, const sfcgal_geometry_t *geom_a,
            const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse angle parameter
  double angle = 0.0;
  if (!parse_double_parameter(op_arg, "angle", 0, 0.0, &angle)) {
    result.error         = true;
    result.error_message = "Invalid angle value";
    return result;
  }

  result.geometry_result = sfcgal_geometry_rotate_z(geom_a, angle);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error rotating geometry around Z axis";
  }

  return result;
}

/**
 * Scale operation with scale parameter
 */
OperationResult
op_scale(const char *op_arg, const sfcgal_geometry_t *geom_a,
         const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse scale parameter
  double scale = 1.0;
  if (!parse_double_parameter(op_arg, "scale", 0, 1.0, &scale)) {
    result.error         = true;
    result.error_message = "Invalid scale value";
    return result;
  }

  result.geometry_result = sfcgal_geometry_scale(geom_a, scale);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error scaling geometry";
  }

  return result;
}

/**
 * Scale 3D operation with sx, sy, sz parameters
 */
OperationResult
op_scale_3d(const char *op_arg, const sfcgal_geometry_t *geom_a,
            const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse parameters
  double sx = 1.0, sy = 1.0, sz = 1.0;

  if (!parse_double_parameter(op_arg, "sx", 0, 1.0, &sx) ||
      !parse_double_parameter(op_arg, "sy", 1, 1.0, &sy) ||
      !parse_double_parameter(op_arg, "sz", 2, 1.0, &sz)) {
    result.error         = true;
    result.error_message = "Invalid parameter value(s)";
    return result;
  }

  result.geometry_result = sfcgal_geometry_scale_3d(geom_a, sx, sy, sz);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error scaling geometry in 3D";
  }

  return result;
}

/**
 * Scale 3D around center operation with sx, sy, sz, cx, cy, cz parameters
 */
OperationResult
op_scale_3d_around_center(const char *op_arg, const sfcgal_geometry_t *geom_a,
                          const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse parameters
  double sx = 1.0, sy = 1.0, sz = 1.0, cx = 0.0, cy = 0.0, cz = 0.0;

  if (!parse_double_parameter(op_arg, "sx", 0, 1.0, &sx) ||
      !parse_double_parameter(op_arg, "sy", 1, 1.0, &sy) ||
      !parse_double_parameter(op_arg, "sz", 2, 1.0, &sz) ||
      !parse_double_parameter(op_arg, "cx", 3, 0.0, &cx) ||
      !parse_double_parameter(op_arg, "cy", 4, 0.0, &cy) ||
      !parse_double_parameter(op_arg, "cz", 5, 0.0, &cz)) {
    result.error         = true;
    result.error_message = "Invalid parameter value(s)";
    return result;
  }

  result.geometry_result =
      sfcgal_geometry_scale_3d_around_center(geom_a, sx, sy, sz, cx, cy, cz);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error scaling geometry in 3D around center";
  }

  return result;
}

/**
 * Translate 2D operation with dx, dy parameters
 */
OperationResult
op_translate_2d(const char *op_arg, const sfcgal_geometry_t *geom_a,
                const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse parameters
  double dx = 0.0, dy = 0.0;

  if (!parse_double_parameter(op_arg, "dx", 0, 0.0, &dx) ||
      !parse_double_parameter(op_arg, "dy", 1, 0.0, &dy)) {
    result.error         = true;
    result.error_message = "Invalid parameter value(s)";
    return result;
  }

  result.geometry_result = sfcgal_geometry_translate_2d(geom_a, dx, dy);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error translating geometry in 2D";
  }

  return result;
}

/**
 * Translate 3D operation with dx, dy, dz parameters
 */
OperationResult
op_translate_3d(const char *op_arg, const sfcgal_geometry_t *geom_a,
                const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse parameters
  double dx = 0.0, dy = 0.0, dz = 0.0;

  if (!parse_double_parameter(op_arg, "dx", 0, 0.0, &dx) ||
      !parse_double_parameter(op_arg, "dy", 1, 0.0, &dy) ||
      !parse_double_parameter(op_arg, "dz", 2, 0.0, &dz)) {
    result.error         = true;
    result.error_message = "Invalid parameter value(s)";
    return result;
  }

  result.geometry_result = sfcgal_geometry_translate_3d(geom_a, dx, dy, dz);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error translating geometry in 3D";
  }

  return result;
}

/**
 * Extrude operation with parameters
 */
OperationResult
op_extrude(const char *op_arg, const sfcgal_geometry_t *geom_a,
           const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse parameters
  double dx = 0.0, dy = 0.0, dz = 1.0;
  if (!parse_double_parameter(op_arg, "dx", 0, 0.0, &dx) ||
      !parse_double_parameter(op_arg, "dy", 1, 0.0, &dy) ||
      !parse_double_parameter(op_arg, "dz", 2, 1.0, &dz)) {
    result.error         = true;
    result.error_message = "Invalid parameter value(s)";
    return result;
  }

  result.geometry_result = sfcgal_geometry_extrude(geom_a, dx, dy, dz);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error extruding geometry";
  }

  return result;
}

/**
 * Round operation with scale parameter
 */
OperationResult
op_round(const char *op_arg, const sfcgal_geometry_t *geom_a,
         const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse scale parameter
  int scale = 0;
  if (!parse_int_parameter(op_arg, "scale", 0, 0, &scale)) {
    result.error         = true;
    result.error_message = "Invalid scale value";
    return result;
  }

  result.geometry_result = sfcgal_geometry_round(geom_a, scale);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error rounding geometry";
  }

  return result;
}

/**
 * Buffer 3D operation with radius, segments, and buffer_type parameters
 */
OperationResult
op_buffer3d(const char *op_arg, const sfcgal_geometry_t *geom_a,
            const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse parameters
  double radius                      = 1.0;
  int    segments                    = 16;
  double buffer_type_val             = 0.0; // Default to SFCGAL_BUFFER3D_ROUND
  sfcgal_buffer3d_type_t buffer_type = SFCGAL_BUFFER3D_ROUND;

  if (!parse_double_parameter(op_arg, "radius", 0, 1.0, &radius) ||
      !parse_int_parameter(op_arg, "segments", 1, 16, &segments) ||
      !parse_double_parameter(op_arg, "buffer_type", 2, 0.0,
                              &buffer_type_val)) {
    result.error         = true;
    result.error_message = "Invalid parameter value(s)";
    return result;
  }

  // Convert buffer_type_val to enum
  int buffer_type_int = (int)buffer_type_val;
  if (buffer_type_int == 0) {
    buffer_type = SFCGAL_BUFFER3D_ROUND;
  } else if (buffer_type_int == 1) {
    buffer_type = SFCGAL_BUFFER3D_CYLSPHERE;
  } else if (buffer_type_int == 2) {
    buffer_type = SFCGAL_BUFFER3D_FLAT;
  } else {
    result.error         = true;
    result.error_message = "Invalid buffer_type value (must be 0, 1, or 2)";
    return result;
  }

  result.geometry_result =
      sfcgal_geometry_buffer3d(geom_a, radius, segments, buffer_type);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error calculating buffer3D";
  }

  return result;
}

/**
 * Offset polygon operation with radius parameter
 */
OperationResult
op_offset_polygon(const char *op_arg, const sfcgal_geometry_t *geom_a,
                  const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse radius parameter
  double radius = 0.0;
  if (!parse_double_parameter(op_arg, "radius", 0, 0.0, &radius)) {
    result.error         = true;
    result.error_message = "Invalid radius value";
    return result;
  }

  if (op_arg == NULL) {
    result.error         = true;
    result.error_message = "Offset polygon requires a radius parameter";
    return result;
  }

  result.geometry_result = sfcgal_geometry_offset_polygon(geom_a, radius);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error calculating offset polygon";
  }

  return result;
}

/**
 * Line substring operation with start and end parameters
 */
OperationResult
op_line_substring(const char *op_arg, const sfcgal_geometry_t *geom_a,
                  const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse parameters
  double start = 0.0, end = 1.0;
  if (!parse_double_parameter(op_arg, "start", 0, 0.0, &start) ||
      !parse_double_parameter(op_arg, "end", 1, 1.0, &end)) {
    result.error         = true;
    result.error_message = "Invalid parameter value(s)";
    return result;
  }

  if (op_arg == NULL) {
    result.error         = true;
    result.error_message = "Line substring requires start and end parameters";
    return result;
  }

  result.geometry_result = sfcgal_geometry_line_sub_string(geom_a, start, end);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error calculating line substring";
  }

  return result;
}

/**
 * Simplify operation with threshold and preserveTopology parameters
 */
OperationResult
op_simplify(const char *op_arg, const sfcgal_geometry_t *geom_a,
            const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type = RESULT_GEOMETRY, .error = false};

  (void)geom_b; // Ignore geom_b

  // Parse parameters
  double threshold         = 1.0;
  bool   preserve_topology = true;

  if (!parse_double_parameter(op_arg, "threshold", 0, 1.0, &threshold) ||
      !parse_bool_parameter(op_arg, "preserve_topology", 1, true,
                            &preserve_topology)) {
    result.error         = true;
    result.error_message = "Invalid parameter value(s)";
    return result;
  }

  result.geometry_result =
      sfcgal_geometry_simplify(geom_a, threshold, preserve_topology);

  if (!result.geometry_result) {
    result.error         = true;
    result.error_message = "Error simplifying geometry";
  }

  return result;
}
