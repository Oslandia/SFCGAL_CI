/**
 * operations_params.c - Implementation of parameter parsing functions
 */

#include "operations_params.h"
#include "../safe_string.h"
#include "../util.h"
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Check if an argument string is in the format name=value
 */
bool
is_named_parameter(const char *arg, const char *name, const char **value_start)
{
  if (!arg || !name || !value_start) {
    return false;
  }

  /* Validate input lengths */
  size_t arg_len  = safe_strlen(arg, MAX_PARAM_LENGTH);
  size_t name_len = safe_strlen(name, MAX_PARAM_LENGTH);

  if (arg_len == SIZE_MAX || name_len == SIZE_MAX || name_len == 0) {
    return false;
  }

  if (arg_len <= name_len + 1) { /* +1 for '=' */
    return false;
  }

  if (strncmp(arg, name, name_len) == 0 && arg[name_len] == '=') {
    *value_start = arg + name_len + 1;
    return true;
  }

  return false;
}

/**
 * Find a parameter value in a comma-separated list of arguments
 */
const char *
find_parameter(const char *arg_list, const char *name, bool *found)
{
  if (!arg_list || !name || !found) {
    return NULL;
  }

  *found = false;

  /* Validate input lengths */
  size_t arg_list_len = safe_strlen(arg_list, MAX_PARAM_LENGTH);
  if (arg_list_len == SIZE_MAX) {
    return NULL;
  }

  /* Make a copy of the argument list since safe_strtok modifies it */
  char *arg_copy = safe_strdup_checked(arg_list, MAX_PARAM_LENGTH);
  if (!arg_copy) {
    return NULL;
  }

  /* Parse comma-separated list using thread-safe tokenizer */
  char       *saveptr     = NULL;
  char       *token       = safe_strtok(arg_copy, ",", &saveptr);
  const char *value_start = NULL;
  char       *result      = NULL;

  while (token) {
    /* Check if this token is in name=value format for our parameter */
    if (is_named_parameter(token, name, &value_start)) {
      *found = true;
      result = safe_strdup_checked(value_start, MAX_PARAM_LENGTH);
      break;
    }
    token = safe_strtok(NULL, ",", &saveptr);
  }

  free(arg_copy);
  return result;
}

/**
 * Parse a double parameter from argument string
 */
bool
parse_double_parameter(const char *arg_list, const char *name, int param_index,
                       double default_value, double *result)
{
  if (!result) {
    return false;
  }

  /* Validate param_index */
  if (param_index < 0 || param_index > 100) { /* Reasonable limit */
    *result = default_value;
    return true;
  }

  /* Set default value initially */
  *result = default_value;

  if (!arg_list) {
    return true; /* No arguments, use default */
  }

  /* Validate input length */
  if (safe_strlen(arg_list, MAX_PARAM_LENGTH) == SIZE_MAX) {
    return false;
  }

  /* Check if we have a simple value (no commas, no equals) */
  if (strchr(arg_list, ',') == NULL && strchr(arg_list, '=') == NULL) {
    /* If this is the first parameter (index 0), use it directly */
    if (param_index == 0) {
      return safe_strtod(arg_list, NULL, result);
    }
    return true; /* Not this parameter, use default */
  }

  /* Otherwise, we need to check for named parameter or positional parameter */
  bool        found     = false;
  const char *value_str = find_parameter(arg_list, name, &found);

  if (found && value_str) {
    /* Parse named parameter */
    bool success = safe_strtod(value_str, NULL, result);
    free((void *)(uintptr_t)value_str);
    return success;
  } else if (!found) {
    /* Try to find by position */
    char *arg_copy = safe_strdup_checked(arg_list, MAX_PARAM_LENGTH);
    if (!arg_copy) {
      return false;
    }

    char *saveptr       = NULL;
    char *token         = safe_strtok(arg_copy, ",", &saveptr);
    int   current_index = 0;

    /* Skip to the desired position */
    while (token && current_index < param_index) {
      token = safe_strtok(NULL, ",", &saveptr);
      current_index++;
    }

    bool success = true;
    /* If we found the position and it's not a named parameter */
    if (token && strchr(token, '=') == NULL) {
      success = safe_strtod(token, NULL, result);
    }

    free(arg_copy);
    return success;
  }

  return true; /* Not found, use default */
}

/**
 * Helper function to parse boolean value from string (forward declaration)
 */
static bool
parse_boolean_value(const char *value_str, bool *result);

/**
 * Parse an integer parameter from argument string
 */
bool
parse_int_parameter(const char *arg_list, const char *name, int param_index,
                    int default_value, int *result)
{
  if (!result) {
    return false;
  }

  /* Parse as double first to check for overflow */
  double double_result;
  bool   success = parse_double_parameter(arg_list, name, param_index,
                                          (double)default_value, &double_result);

  if (success) {
    /* Check for integer overflow */
    if (double_result > (double)INT_MAX || double_result < (double)INT_MIN) {
      errno   = ERANGE;
      *result = default_value;
      return false;
    }

    /* Convert to int with proper bounds checking */
    *result = (int)double_result;
  } else {
    *result = default_value;
  }

  return success;
}

/**
 * Parse a boolean parameter from argument string
 */
bool
parse_bool_parameter(const char *arg_list, const char *name, int param_index,
                     bool default_value, bool *result)
{
  if (!result) {
    return false;
  }

  /* Validate param_index */
  if (param_index < 0 || param_index > 100) { /* Reasonable limit */
    *result = default_value;
    return true;
  }

  /* Set default value initially */
  *result = default_value;

  if (!arg_list) {
    return true; /* No arguments, use default */
  }

  /* Validate input length */
  if (safe_strlen(arg_list, MAX_PARAM_LENGTH) == SIZE_MAX) {
    return false;
  }

  /* Check if we have a simple value (no commas, no equals) */
  if (strchr(arg_list, ',') == NULL && strchr(arg_list, '=') == NULL) {
    /* If this is the first parameter (index 0), use it directly */
    if (param_index == 0) {
      return parse_boolean_value(arg_list, result);
    }
    return true; /* Not this parameter, use default */
  }

  /* Otherwise, we need to check for named parameter or positional parameter */
  bool        found     = false;
  const char *value_str = find_parameter(arg_list, name, &found);

  if (found && value_str) {
    /* Parse named parameter */
    bool success = parse_boolean_value(value_str, result);
    free((void *)(uintptr_t)value_str);
    return success;
  } else if (!found) {
    /* Try to find by position */
    char *arg_copy = safe_strdup_checked(arg_list, MAX_PARAM_LENGTH);
    if (!arg_copy) {
      return false;
    }

    char *saveptr       = NULL;
    char *token         = safe_strtok(arg_copy, ",", &saveptr);
    int   current_index = 0;

    /* Skip to the desired position */
    while (token && current_index < param_index) {
      token = safe_strtok(NULL, ",", &saveptr);
      current_index++;
    }

    bool success = true;
    /* If we found the position and it's not a named parameter */
    if (token && strchr(token, '=') == NULL) {
      success = parse_boolean_value(token, result);
    }

    free(arg_copy);
    return success;
  }

  return true; /* Not found, use default */
}

/**
 * Helper function to parse boolean value from string
 */
static bool
parse_boolean_value(const char *value_str, bool *result)
{
  if (!value_str || !result) {
    return false;
  }

  /* Check for boolean string values */
  if (strcmp(value_str, "1") == 0 || strcasecmp(value_str, "true") == 0 ||
      strcasecmp(value_str, "yes") == 0 || strcasecmp(value_str, "on") == 0) {
    *result = true;
    return true;
  } else if (strcmp(value_str, "0") == 0 ||
             strcasecmp(value_str, "false") == 0 ||
             strcasecmp(value_str, "no") == 0 ||
             strcasecmp(value_str, "off") == 0) {
    *result = false;
    return true;
  }

  return false; /* Unrecognized value */
}
