/**
 * operations.h - Function declarations for SFCGAL operations
 */

#ifndef OPERATIONS_H
#define OPERATIONS_H

#include "sfcgal_c.h"
#include <stdbool.h>

/**
 * Result type enumeration
 */
typedef enum {
  RESULT_GEOMETRY, // Result is a geometry
  RESULT_BOOLEAN,  // Result is a boolean
  RESULT_NUMERIC,  // Result is a numeric value
  RESULT_TEXT      // Result is a text string
} ResultType;

/**
 * Operation result structure
 */
typedef struct {
  ResultType type;
  union {
    sfcgal_geometry_t *geometry_result;
    bool               boolean_result;
    double             numeric_result;
    char              *text_result;
  };
  bool        error;
  const char *error_message;
} OperationResult;

/**
 * Parameter structure for operations
 */
typedef struct {
  const char *name;
  const char *description;
  double      default_value;
} OperationParam;

/**
 * Operation function typedef
 */
typedef OperationResult (*OperationFunction)(const char              *op_arg,
                                             const sfcgal_geometry_t *geom_a,
                                             const sfcgal_geometry_t *geom_b);

/**
 * Operation structure to map operation names to functions
 */
typedef struct {
  const char       *name;  // Name of the operation
  OperationFunction func;  // Function pointer for the operation
  bool requires_geom_b;    // Whether the operation requires a second geometry
  const char *description; // Description of the operation
  const char *category;    // Category of the operation
  ResultType  result_type; // Type of result the operation returns
  bool requires_arg; // Whether the operation requires an additional argument
  const OperationParam *params;      // Array of parameters
  int                   param_count; // Number of parameters
} Operation;

/**
 * Execute an SFCGAL operation
 *
 * @param op_name The name of the operation to execute
 * @param op_arg Additional argument for the operation
 * @param geom_a The first geometry
 * @param geom_b The second geometry (may be NULL)
 * @return OperationResult structure containing the result
 */
OperationResult
execute_operation(const char *op_name, const char *op_arg,
                  const sfcgal_geometry_t *geom_a,
                  const sfcgal_geometry_t *geom_b);

/**
 * Find an operation by name
 *
 * @param name Operation name
 * @return Pointer to the Operation or NULL if not found
 */
const Operation *
find_operation(const char *name);

/**
 * Print available operations to stdout
 */
void
print_available_operations(void);

/**
 * Print help for a specific operation
 *
 * @param name Operation name
 * @return true if the operation was found and help printed, false otherwise
 */
bool
print_operation_help(const char *name);

/**
 * Check if an operation requires arguments
 *
 * @param name Operation name
 * @return true if the operation requires arguments, false otherwise
 */
bool
operation_requires_arg(const char *name);

#endif /* OPERATIONS_H */
