/**
 * operations_params.c - Implementation of parameter parsing functions
 */

#include "operations_params.h"
#include "../util.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/**
 * Check if an argument string is in the format name=value
 */
bool is_named_parameter(const char* arg, const char* name, const char** value_start) {
    if (!arg || !name || !value_start) {
        return false;
    }
    
    size_t name_len = strlen(name);
    if (strncmp(arg, name, name_len) == 0 && arg[name_len] == '=') {
        *value_start = arg + name_len + 1;
        return true;
    }
    
    return false;
}

/**
 * Find a parameter value in a comma-separated list of arguments
 */
const char* find_parameter(const char* arg_list, const char* name, bool* found) {
    if (!arg_list || !name || !found) {
        *found = false;
        return NULL;
    }
    
    *found = false;
    
    // Make a copy of the argument list since strtok modifies it
    char* arg_copy = strdup(arg_list);
    if (!arg_copy) {
        return NULL;
    }
    
    // Parse comma-separated list
    char* token = strtok(arg_copy, ",");
    const char* value_start = NULL;
    
    while (token) {
        // Check if this token is in name=value format for our parameter
        if (is_named_parameter(token, name, &value_start)) {
            *found = true;
            break;
        }
        token = strtok(NULL, ",");
    }
    
    // If found, duplicate the value part
    const char* result = NULL;
    if (*found) {
        result = strdup(value_start);
    }
    
    free(arg_copy);
    return result;
}

/**
 * Parse a double parameter from argument string
 */
bool parse_double_parameter(const char* arg_list, const char* name, 
                          int param_index, double default_value, double* result) {
    if (!result) {
        return false;
    }
    
    // Set default value initially
    *result = default_value;
    
    if (!arg_list) {
        return true; // No arguments, use default
    }
    
    // First check if we have a simple value (no commas, no equals)
    if (strchr(arg_list, ',') == NULL && strchr(arg_list, '=') == NULL) {
        // If this is the first parameter (index 0), use it directly
        if (param_index == 0) {
            char* endptr;
            *result = strtod(arg_list, &endptr);
            return (endptr != arg_list); // Success if strtod parsed something
        }
        return true; // Not this parameter, use default
    }
    
    // Otherwise, we need to check for named parameter or positional parameter
    bool found = false;
    const char* value_str = find_parameter(arg_list, name, &found);
    
    if (found && value_str) {
        // Parse named parameter
        char* endptr;
        *result = strtod(value_str, &endptr);
        free((void*)value_str);
        return (endptr != value_str); // Success if strtod parsed something
    } else if (!found) {
        // Try to find by position
        char* arg_copy = strdup(arg_list);
        if (!arg_copy) {
            return false;
        }
        
        char* token = strtok(arg_copy, ",");
        int current_index = 0;
        
        // Skip to the desired position
        while (token && current_index < param_index) {
            token = strtok(NULL, ",");
            current_index++;
        }
        
        // If we found the position and it's not a named parameter
        if (token && strchr(token, '=') == NULL) {
            char* endptr;
            *result = strtod(token, &endptr);
            free(arg_copy);
            return (endptr != token); // Success if strtod parsed something
        }
        
        free(arg_copy);
    }
    
    return true; // Not found, use default
}

/**
 * Parse an integer parameter from argument string
 */
bool parse_int_parameter(const char* arg_list, const char* name,
                       int param_index, int default_value, int* result) {
    if (!result) {
        return false;
    }
    
    // Parse as double first
    double double_result;
    bool success = parse_double_parameter(arg_list, name, param_index, 
                                         (double)default_value, &double_result);
    
    if (success) {
        *result = (int)double_result;
    } else {
        *result = default_value;
    }
    
    return success;
}

/**
 * Parse a boolean parameter from argument string
 */
bool parse_bool_parameter(const char* arg_list, const char* name,
                        int param_index, bool default_value, bool* result) {
    if (!result) {
        return false;
    }
    
    // Set default value initially
    *result = default_value;
    
    if (!arg_list) {
        return true; // No arguments, use default
    }
    
    // First check if we have a simple value (no commas, no equals)
    if (strchr(arg_list, ',') == NULL && strchr(arg_list, '=') == NULL) {
        // If this is the first parameter (index 0), use it directly
        if (param_index == 0) {
            // Check for boolean string values
            if (strcmp(arg_list, "1") == 0 ||
                strcasecmp(arg_list, "true") == 0 ||
                strcasecmp(arg_list, "yes") == 0 ||
                strcasecmp(arg_list, "on") == 0) {
                *result = true;
                return true;
            } else if (strcmp(arg_list, "0") == 0 ||
                       strcasecmp(arg_list, "false") == 0 ||
                       strcasecmp(arg_list, "no") == 0 ||
                       strcasecmp(arg_list, "off") == 0) {
                *result = false;
                return true;
            }
            return false; // Unrecognized value
        }
        return true; // Not this parameter, use default
    }
    
    // Otherwise, we need to check for named parameter or positional parameter
    bool found = false;
    const char* value_str = find_parameter(arg_list, name, &found);
    
    if (found && value_str) {
        // Check for boolean string values
        if (strcmp(value_str, "1") == 0 ||
            strcasecmp(value_str, "true") == 0 ||
            strcasecmp(value_str, "yes") == 0 ||
            strcasecmp(value_str, "on") == 0) {
            *result = true;
            free((void*)value_str);
            return true;
        } else if (strcmp(value_str, "0") == 0 ||
                  strcasecmp(value_str, "false") == 0 ||
                  strcasecmp(value_str, "no") == 0 ||
                  strcasecmp(value_str, "off") == 0) {
            *result = false;
            free((void*)value_str);
            return true;
        }
        free((void*)value_str);
        return false; // Unrecognized value
    } else if (!found) {
        // Try to find by position
        char* arg_copy = strdup(arg_list);
        if (!arg_copy) {
            return false;
        }
        
        char* token = strtok(arg_copy, ",");
        int current_index = 0;
        
        // Skip to the desired position
        while (token && current_index < param_index) {
            token = strtok(NULL, ",");
            current_index++;
        }
        
        // If we found the position and it's not a named parameter
        if (token && strchr(token, '=') == NULL) {
            // Check for boolean string values
            if (strcmp(token, "1") == 0 ||
                strcasecmp(token, "true") == 0 ||
                strcasecmp(token, "yes") == 0 ||
                strcasecmp(token, "on") == 0) {
                *result = true;
                free(arg_copy);
                return true;
            } else if (strcmp(token, "0") == 0 ||
                      strcasecmp(token, "false") == 0 ||
                      strcasecmp(token, "no") == 0 ||
                      strcasecmp(token, "off") == 0) {
                *result = false;
                free(arg_copy);
                return true;
            }
            free(arg_copy);
            return false; // Unrecognized value
        }
        
        free(arg_copy);
    }
    
    return true; // Not found, use default
}
