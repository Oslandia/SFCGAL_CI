/**
 * operations.c - Implementation of SFCGAL operations
 */

#include "operations.h"
#include "operations/operations_params.h"
#include "operations/operations_registry.h"
#include "safe_string.h"
#include "util.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Find an operation by name
 */
const Operation *
find_operation(const char *name)
{
  if (!name) {
    return NULL;
  }

  for (int i = 0; operations[i].name != NULL; i++) {
    if (strcasecmp(operations[i].name, name) == 0) {
      return &operations[i];
    }
  }

  return NULL;
}

/**
 * Execute an SFCGAL operation
 */
OperationResult
execute_operation(const char *op_name, const char *op_arg,
                  const sfcgal_geometry_t *geom_a,
                  const sfcgal_geometry_t *geom_b)
{
  OperationResult result = {.type           = RESULT_BOOLEAN,
                            .boolean_result = false,
                            .error          = true,
                            .error_message  = "Unknown operation"};

  if (!op_name) {
    return result;
  }

  // Find operation
  const Operation *op = find_operation(op_name);
  if (!op) {
    return result;
  }

  // Check if operation requires geom_b
  if (op->requires_geom_b && !geom_b) {
    result.error_message = "Operation requires a second geometry";
    return result;
  }

  // Execute operation
  return op->func(op_arg, geom_a, geom_b);
}

/**
 * Check if an operation requires arguments
 */
bool
operation_requires_arg(const char *name)
{
  const Operation *op = find_operation(name);

  if (!op) {
    // Operation not found, assume it doesn't require an argument
    return false;
  }

  return op->requires_arg;
}

/**
 * Helper function to get a string representation of the input type for an
 * operation
 *
 * @param op Operation to get input type for
 * @return String describing the input type(s)
 */
static const char *
get_input_type_str(const Operation *op)
{
  if (op->requires_geom_b) {
    return op->requires_arg ? "A, B, params" : "A, B";
  } else {
    return op->requires_arg ? "A, params" : "A";
  }
}

/**
 * Helper function to get a string representation of the output type for an
 * operation
 *
 * @param op Operation to get output type for
 * @return String describing the output type
 */
static const char *
get_output_type_str(const Operation *op)
{
  switch (op->result_type) {
  case RESULT_GEOMETRY:
    return "G";
  case RESULT_BOOLEAN:
    return "B";
  case RESULT_NUMERIC:
    return "D";
  case RESULT_TEXT:
    return "T";
  default:
    return "?";
  }
}

/**
 * Get array of unique categories
 *
 * @param num_categories Pointer to store the number of categories
 * @return Array of category strings (must be freed by caller)
 */
static char **
get_unique_categories(int *num_categories)
{
  if (!num_categories) {
    return NULL;
  }

  // First, count unique categories
  int max_categories = 0;
  for (int i = 0; operations[i].name != NULL; i++) {
    max_categories++;
  }

  // We know max_categories is positive since there are operations defined
  char **categories = (char **)malloc((size_t)max_categories * sizeof(char *));
  if (!categories) {
    *num_categories = 0;
    return NULL;
  }

  *num_categories = 0;

  // Add each unique category
  for (int i = 0; operations[i].name != NULL; i++) {
    bool found = false;

    // Check if category is already in list
    for (int j = 0; j < *num_categories; j++) {
      if (strcmp(categories[j], operations[i].category) == 0) {
        found = true;
        break;
      }
    }

    // If not found, add it
    if (!found) {
      categories[*num_categories] =
          safe_strdup(operations[i].category, SAFE_MAX_STRING_LENGTH);
      (*num_categories)++;
    }
  }

  return categories;
}

/**
 * Calculate the maximum widths for each column in the operations table
 */
static void
calculate_column_widths(size_t *max_category_width, size_t *max_operation_width,
                        size_t *max_input_width, size_t *max_output_width,
                        size_t *max_description_width)
{
  if (!max_category_width || !max_operation_width || !max_input_width ||
      !max_output_width || !max_description_width) {
    return;
  }

  // Initialize with minimum widths
  *max_category_width    = 8;  // "Category"
  *max_operation_width   = 9;  // "Operation"
  *max_input_width       = 5;  // "Input"
  *max_output_width      = 6;  // "Output"
  *max_description_width = 11; // "Description"

  // Get unique categories
  int    num_categories = 0;
  char **categories     = get_unique_categories(&num_categories);

  if (!categories) {
    return;
  }

  // Check category width
  for (int i = 0; i < num_categories; i++) {
    size_t category_len = safe_strlen(categories[i], SAFE_MAX_STRING_LENGTH);
    if (category_len == SIZE_MAX) {
      continue; // Skip invalid category names
    }
    if (category_len > *max_category_width) {
      *max_category_width = category_len;
    }
  }

  // Free categories now that we're done with them
  free(categories);

  // Check operation width, input width, output width, and description width
  for (int i = 0; operations[i].name != NULL; i++) {
    // Operation name
    size_t name_len = safe_strlen(operations[i].name, SAFE_MAX_STRING_LENGTH);
    if (name_len == SIZE_MAX) {
      continue; // Skip invalid operation names
    }
    if (name_len > *max_operation_width) {
      *max_operation_width = name_len;
    }

    // Input type
    const char *input_type = get_input_type_str(&operations[i]);
    size_t      input_len  = safe_strlen(input_type, SAFE_MAX_STRING_LENGTH);
    if (input_len == SIZE_MAX) {
      continue; // Skip invalid input type strings
    }
    if (input_len > *max_input_width) {
      *max_input_width = input_len;
    }

    // Output type (single character, no need to check)

    // Description
    size_t desc_len =
        safe_strlen(operations[i].description, SAFE_MAX_STRING_LENGTH);
    if (desc_len == SIZE_MAX) {
      continue; // Skip invalid descriptions
    }
    if (desc_len > *max_description_width) {
      *max_description_width = desc_len;
    }
  }

  // Add padding
  *max_category_width += 2;
  *max_operation_width += 2;
  *max_input_width += 2;
  *max_output_width += 2;
  *max_description_width += 2;
}

/**
 * Print available operations to stdout in a well-formatted table
 */
void
print_available_operations(void)
{
  // Calculate column widths
  size_t cat_width, op_width, in_width, out_width, desc_width;
  calculate_column_widths(&cat_width, &op_width, &in_width, &out_width,
                          &desc_width);

  // Calculate total width for horizontal lines
  size_t total_width = cat_width + op_width + in_width + out_width +
                       desc_width + 6; // 6 for the vertical bars

  // Print top border
  printf("┌");
  for (size_t i = 0; i < cat_width; i++)
    printf("─");
  printf("┬");
  for (size_t i = 0; i < op_width; i++)
    printf("─");
  printf("┬");
  for (size_t i = 0; i < in_width; i++)
    printf("─");
  printf("┬");
  for (size_t i = 0; i < out_width; i++)
    printf("─");
  printf("┬");
  for (size_t i = 0; i < desc_width; i++)
    printf("─");
  printf("┐\n");

  // Print header
  printf("│ %-*s│ %-*s│ %-*s│ %-*s│ %-*s│\n", (int)cat_width - 1, "Category",
         (int)op_width - 1, "Operation", (int)in_width - 1, "Input",
         (int)out_width - 1, "Output", (int)desc_width - 1, "Description");

  // Print header separator
  printf("├");
  for (size_t i = 0; i < cat_width; i++)
    printf("─");
  printf("┼");
  for (size_t i = 0; i < op_width; i++)
    printf("─");
  printf("┼");
  for (size_t i = 0; i < in_width; i++)
    printf("─");
  printf("┼");
  for (size_t i = 0; i < out_width; i++)
    printf("─");
  printf("┼");
  for (size_t i = 0; i < desc_width; i++)
    printf("─");
  printf("┤\n");

  // Get unique categories
  int    num_categories = 0;
  char **categories     = get_unique_categories(&num_categories);

  if (!categories) {
    printf("│ Error generating categories table%*s│\n", (int)(total_width - 34),
           "");
    printf("└");
    for (size_t i = 0; i < total_width - 2; i++)
      printf("─");
    printf("┘\n");
    return;
  }

  // Print operations by category
  for (int cat_idx = 0; cat_idx < num_categories; cat_idx++) {
    const char *category = categories[cat_idx];
    int         printed  = 0;

    // Print operations in this category
    for (int i = 0; operations[i].name != NULL; i++) {
      if (strcmp(operations[i].category, category) == 0) {
        if (printed == 0) {
          // First operation in category - print category name
          printf("│ %-*s│ %-*s│ %-*s│ %-*s│ %-*s│\n", (int)cat_width - 1,
                 category, (int)op_width - 1, operations[i].name,
                 (int)in_width - 1, get_input_type_str(&operations[i]),
                 (int)out_width - 1, get_output_type_str(&operations[i]),
                 (int)desc_width - 1, operations[i].description);
        } else {
          // Subsequent operations - leave category column blank
          printf("│ %-*s│ %-*s│ %-*s│ %-*s│ %-*s│\n", (int)cat_width - 1, "",
                 (int)op_width - 1, operations[i].name, (int)in_width - 1,
                 get_input_type_str(&operations[i]), (int)out_width - 1,
                 get_output_type_str(&operations[i]), (int)desc_width - 1,
                 operations[i].description);
        }
        printed++;
      }
    }

    // Print separator between categories, unless it's the last category
    if (cat_idx < num_categories - 1) {
      printf("├");
      for (size_t i = 0; i < cat_width; i++)
        printf("─");
      printf("┼");
      for (size_t i = 0; i < op_width; i++)
        printf("─");
      printf("┼");
      for (size_t i = 0; i < in_width; i++)
        printf("─");
      printf("┼");
      for (size_t i = 0; i < out_width; i++)
        printf("─");
      printf("┼");
      for (size_t i = 0; i < desc_width; i++)
        printf("─");
      printf("┤\n");
    }
  }

  // Print table footer
  printf("└");
  for (size_t i = 0; i < cat_width; i++)
    printf("─");
  printf("┴");
  for (size_t i = 0; i < op_width; i++)
    printf("─");
  printf("┴");
  for (size_t i = 0; i < in_width; i++)
    printf("─");
  printf("┴");
  for (size_t i = 0; i < out_width; i++)
    printf("─");
  printf("┴");
  for (size_t i = 0; i < desc_width; i++)
    printf("─");
  printf("┘\n");

  // Legend
  printf("\nLegend:\n");
  printf("  A, B: Input geometries      G: Output geometry      D: Numeric "
         "result\n");
  printf("  T: Text result              params: Parameter list\n");

  // Free allocated memory
  free(categories);
}

/**
 * Print help for a specific operation
 */
bool
print_operation_help(const char *name)
{
  if (!name) {
    return false;
  }

  const Operation *op = find_operation(name);
  if (!op) {
    printf("Operation '%s' not found.\n", name);
    return false;
  }

  printf("Operation: %s\n", op->name);
  printf("Category: %s\n", op->category);
  printf("Description: %s\n", op->description);
  printf("Input: %s\n", get_input_type_str(op));
  printf("Output: %s (%s)\n", get_output_type_str(op),
         op->result_type == RESULT_GEOMETRY  ? "Geometry"
         : op->result_type == RESULT_BOOLEAN ? "Boolean"
         : op->result_type == RESULT_NUMERIC ? "Numeric"
         : op->result_type == RESULT_TEXT    ? "Text"
                                             : "Unknown");

  if (op->param_count > 0 && op->params) {
    printf("\nParameters:\n");
    for (int i = 0; i < op->param_count; i++) {
      printf("  %s: %s (default: %.6g)\n", op->params[i].name,
             op->params[i].description, op->params[i].default_value);
    }
    printf("\nParameters can be specified by position (comma-separated) or by "
           "name (name=value).\n");
    printf("Example: %s ", op->name);

    // Generate example parameter string
    for (int i = 0; i < op->param_count; i++) {
      // Use something other than default for the example
      double example_value = op->params[i].default_value + 1.0;
      if (i > 0) {
        printf(",");
      }
      printf("%s=%.1f", op->params[i].name, example_value);
    }
    printf("\n");
  }

  return true;
}
