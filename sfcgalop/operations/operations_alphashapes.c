/**
 * operations_alphashapes.c - Implementation of alpha shapes operations
 */

#include "operations_alphashapes.h"
#include "operations_common.h"
#include "../util.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

#if !defined(_MSC_VER)
/**
 * Alpha shapes operation with alpha and allow_holes parameters
 */
OperationResult op_alpha_shapes(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Parse parameters
    double alpha = 0.0;
    bool allow_holes = false;
    
    if (!parse_double_parameter(op_arg, "alpha", 0, 0.0, &alpha) ||
        !parse_bool_parameter(op_arg, "allow_holes", 1, false, &allow_holes)) {
        result.error = true;
        result.error_message = "Invalid parameter value(s)";
        return result;
    }
    
    result.geometry_result = sfcgal_geometry_alpha_shapes(geom_a, alpha, allow_holes);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error calculating alpha shapes";
    }
    
    return result;
}

/**
 * Optimal alpha shapes operation with allow_holes and nb_components parameters
 */
OperationResult op_optimal_alpha_shapes(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Parse parameters
    bool allow_holes = false;
    int nb_components_int = 1;
    size_t nb_components = 1;
    
    if (!parse_bool_parameter(op_arg, "allow_holes", 0, false, &allow_holes) ||
        !parse_int_parameter(op_arg, "nb_components", 1, 1, &nb_components_int)) {
        result.error = true;
        result.error_message = "Invalid parameter value(s)";
        return result;
    }
    
    if (nb_components_int < 1) {
        result.error = true;
        result.error_message = "nb_components must be positive";
        return result;
    }
    nb_components = (size_t)nb_components_int;
    
    result.geometry_result = sfcgal_geometry_optimal_alpha_shapes(geom_a, allow_holes, nb_components);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error calculating optimal alpha shapes";
    }
    
    return result;
}
#endif

/**
 * Alpha wrapping 3D operation with relativeAlpha and relativeOffset parameters
 */
OperationResult op_alpha_wrapping_3d(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Parse parameters
    int relative_alpha_int = 100;
    int relative_offset_int = 100;
    size_t relative_alpha = 100;
    size_t relative_offset = 100;
    
    if (!parse_int_parameter(op_arg, "relative_alpha", 0, 100, &relative_alpha_int) ||
        !parse_int_parameter(op_arg, "relative_offset", 1, 100, &relative_offset_int)) {
        result.error = true;
        result.error_message = "Invalid parameter value(s)";
        return result;
    }
    
    if (relative_alpha_int < 1 || relative_offset_int < 1) {
        result.error = true;
        result.error_message = "relative_alpha and relative_offset must be positive";
        return result;
    }
    relative_alpha = (size_t)relative_alpha_int;
    relative_offset = (size_t)relative_offset_int;
    
    result.geometry_result = sfcgal_geometry_alpha_wrapping_3d(geom_a, relative_alpha, relative_offset);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error calculating alpha wrapping 3D";
    }
    
    return result;
}
