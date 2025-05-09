/**
 * operations_common.h - Common utilities for all SFCGAL operations
 */

#ifndef OPERATIONS_COMMON_H
#define OPERATIONS_COMMON_H

#include "operations.h"

/* ========================================================================== */
/* MACROS FOR GENERATING OPERATION FUNCTIONS */
/* ========================================================================== */

/**
 * Macro for operations that return a boolean result (single geometry input)
 */
#define DEFINE_BOOLEAN_OP_SINGLE(name, func)                                   \
  OperationResult op_##name(const char              *op_arg,                   \
                            const sfcgal_geometry_t *geom_a,                   \
                            const sfcgal_geometry_t *geom_b)                   \
  {                                                                            \
    OperationResult result = {.type = RESULT_BOOLEAN, .error = false};         \
    (void)op_arg; /* Ignore op_arg */                                          \
    (void)geom_b; /* Ignore geom_b */                                          \
    int ret = func(geom_a);                                                    \
    if (ret < 0) {                                                             \
      result.error         = true;                                             \
      result.error_message = "Error executing " #name;                         \
      return result;                                                           \
    }                                                                          \
    result.boolean_result = (ret == 1);                                        \
    return result;                                                             \
  }

/**
 * Macro for operations that return a boolean result (two geometry inputs)
 */
#define DEFINE_BOOLEAN_OP_BINARY(name, func)                                   \
  OperationResult op_##name(const char              *op_arg,                   \
                            const sfcgal_geometry_t *geom_a,                   \
                            const sfcgal_geometry_t *geom_b)                   \
  {                                                                            \
    OperationResult result = {.type = RESULT_BOOLEAN, .error = false};         \
    (void)op_arg; /* Ignore op_arg */                                          \
    int ret = func(geom_a, geom_b);                                            \
    if (ret < 0) {                                                             \
      result.error         = true;                                             \
      result.error_message = "Error executing " #name;                         \
      return result;                                                           \
    }                                                                          \
    result.boolean_result = (ret == 1);                                        \
    return result;                                                             \
  }

/**
 * Macro for operations that return a numeric result (single geometry input)
 */
#define DEFINE_NUMERIC_OP_SINGLE(name, func)                                   \
  OperationResult op_##name(const char              *op_arg,                   \
                            const sfcgal_geometry_t *geom_a,                   \
                            const sfcgal_geometry_t *geom_b)                   \
  {                                                                            \
    OperationResult result = {.type = RESULT_NUMERIC, .error = false};         \
    (void)op_arg; /* Ignore op_arg */                                          \
    (void)geom_b; /* Ignore geom_b */                                          \
    double value = func(geom_a);                                               \
    if (value < 0 && !sfcgal_geometry_is_empty(geom_a)) {                      \
      result.error         = true;                                             \
      result.error_message = "Error calculating " #name;                       \
      return result;                                                           \
    }                                                                          \
    result.numeric_result = value;                                             \
    return result;                                                             \
  }

/**
 * Macro for operations that return a numeric result (two geometry inputs)
 */
#define DEFINE_NUMERIC_OP_BINARY(name, func)                                   \
  OperationResult op_##name(const char              *op_arg,                   \
                            const sfcgal_geometry_t *geom_a,                   \
                            const sfcgal_geometry_t *geom_b)                   \
  {                                                                            \
    OperationResult result = {.type = RESULT_NUMERIC, .error = false};         \
    (void)op_arg; /* Ignore op_arg */                                          \
    result.numeric_result = func(geom_a, geom_b);                              \
    return result;                                                             \
  }

/**
 * Macro for operations that return a geometry (single geometry input)
 */
#define DEFINE_GEOMETRY_OP_SINGLE(name, func)                                  \
  OperationResult op_##name(const char              *op_arg,                   \
                            const sfcgal_geometry_t *geom_a,                   \
                            const sfcgal_geometry_t *geom_b)                   \
  {                                                                            \
    OperationResult result = {.type = RESULT_GEOMETRY, .error = false};        \
    (void)op_arg; /* Ignore op_arg */                                          \
    (void)geom_b; /* Ignore geom_b */                                          \
    result.geometry_result = func(geom_a);                                     \
    if (!result.geometry_result) {                                             \
      result.error         = true;                                             \
      result.error_message = "Error calculating " #name;                       \
    }                                                                          \
    return result;                                                             \
  }

/**
 * Macro for operations that return a geometry (two geometry inputs)
 */
#define DEFINE_GEOMETRY_OP_BINARY(name, func)                                  \
  OperationResult op_##name(const char              *op_arg,                   \
                            const sfcgal_geometry_t *geom_a,                   \
                            const sfcgal_geometry_t *geom_b)                   \
  {                                                                            \
    OperationResult result = {.type = RESULT_GEOMETRY, .error = false};        \
    (void)op_arg; /* Ignore op_arg */                                          \
    result.geometry_result = func(geom_a, geom_b);                             \
    if (!result.geometry_result) {                                             \
      result.error         = true;                                             \
      result.error_message = "Error calculating " #name;                       \
    }                                                                          \
    return result;                                                             \
  }

/**
 * Macro for operations that modify a geometry in place and return a clone
 */
#define DEFINE_INPLACE_OP(name, func)                                          \
  OperationResult op_##name(const char              *op_arg,                   \
                            const sfcgal_geometry_t *geom_a,                   \
                            const sfcgal_geometry_t *geom_b)                   \
  {                                                                            \
    OperationResult result = {.type = RESULT_GEOMETRY, .error = false};        \
    (void)op_arg; /* Ignore op_arg */                                          \
    (void)geom_b; /* Ignore geom_b */                                          \
    /* Clone geometry first to avoid modifying input */                        \
    result.geometry_result = sfcgal_geometry_clone(geom_a);                    \
    if (!result.geometry_result) {                                             \
      result.error         = true;                                             \
      result.error_message = "Failed to clone geometry";                       \
      return result;                                                           \
    }                                                                          \
    int changed = func(result.geometry_result);                                \
    if (changed < 0) {                                                         \
      result.error         = true;                                             \
      result.error_message = "Error executing " #name;                         \
      sfcgal_geometry_delete(result.geometry_result);                          \
      result.geometry_result = NULL;                                           \
    }                                                                          \
    return result;                                                             \
  }

#endif /* OPERATIONS_COMMON_H */
