/**
 * operations_params.h - Parameter parsing for SFCGAL operations (ENHANCED
 * VERSION)
 */

#ifndef OPERATIONS_PARAMS_H
#define OPERATIONS_PARAMS_H

#include <stdbool.h>

/**
 * Check if an argument string is in the format name=value
 *
 * @param arg Argument string to check
 * @param name Parameter name to look for
 * @param value_start Pointer to store the start of value part
 * @return true if the argument is in format name=value, false otherwise
 */
bool
is_named_parameter(const char *arg, const char *name, const char **value_start);

/**
 * Find a parameter value in a comma-separated list of arguments
 *
 * @param arg_list Full argument string with comma-separated parameters
 * @param name Parameter name to find
 * @param found Pointer to store whether parameter was found or not
 * @return Parameter value as string (must be freed by caller), or NULL if not
 * found
 */
const char *
find_parameter(const char *arg_list, const char *name, bool *found);

/**
 * Parse a double parameter from argument string
 *
 * @param arg_list Full argument string (comma-separated or single value)
 * @param name Parameter name to look for
 * @param param_index Position of parameter in positional arguments (0-based)
 * @param default_value Default value to use if parameter not found
 * @param result Pointer to store the resulting value
 * @return true if parameter was successfully parsed, false on error
 */
bool
parse_double_parameter(const char *arg_list, const char *name, int param_index,
                       double default_value, double *result);

/**
 * Parse an integer parameter from argument string
 *
 * @param arg_list Full argument string (comma-separated or single value)
 * @param name Parameter name to look for
 * @param param_index Position of parameter in positional arguments (0-based)
 * @param default_value Default value to use if parameter not found
 * @param result Pointer to store the resulting value
 * @return true if parameter was successfully parsed, false on error
 */
bool
parse_int_parameter(const char *arg_list, const char *name, int param_index,
                    int default_value, int *result);

/**
 * Parse a boolean parameter from argument string
 *
 * @param arg_list Full argument string (comma-separated or single value)
 * @param name Parameter name to look for
 * @param param_index Position of parameter in positional arguments (0-based)
 * @param default_value Default value to use if parameter not found
 * @param result Pointer to store the resulting value
 * @return true if parameter was successfully parsed, false on error
 */
bool
parse_bool_parameter(const char *arg_list, const char *name, int param_index,
                     bool default_value, bool *result);

/**
 * Parse a boolean value from a string
 *
 * @param value_str String value to parse
 * @param result Pointer to store the resulting boolean value
 * @return true if parsing was successful, false on error
 */
bool
parse_bool_value(const char *value_str, bool *result);

#endif /* OPERATIONS_PARAMS_H */
