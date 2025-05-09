/**
 * operations.c - Implementation of SFCGAL operations
 */

#include "operations.h"
#include "util.h"
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* ========================================================================== */
/* MACROS FOR GENERATING OPERATION FUNCTIONS */
/* ========================================================================== */

/**
 * Macro for operations that return a boolean result (single geometry input)
 */
#define DEFINE_BOOLEAN_OP_SINGLE(name, func) \
static OperationResult op_##name(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) { \
    OperationResult result = { \
        .type = RESULT_BOOLEAN, \
        .error = false \
    }; \
    (void)op_arg; /* Ignore op_arg */ \
    (void)geom_b; /* Ignore geom_b */ \
    int ret = func(geom_a); \
    if (ret < 0) { \
        result.error = true; \
        result.error_message = "Error executing " #name; \
        return result; \
    } \
    result.boolean_result = (ret == 1); \
    return result; \
}

/**
 * Macro for operations that return a boolean result (two geometry inputs)
 */
#define DEFINE_BOOLEAN_OP_BINARY(name, func) \
static OperationResult op_##name(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) { \
    OperationResult result = { \
        .type = RESULT_BOOLEAN, \
        .error = false \
    }; \
    (void)op_arg; /* Ignore op_arg */ \
    int ret = func(geom_a, geom_b); \
    if (ret < 0) { \
        result.error = true; \
        result.error_message = "Error executing " #name; \
        return result; \
    } \
    result.boolean_result = (ret == 1); \
    return result; \
}

/**
 * Macro for operations that return a numeric result (single geometry input)
 */
#define DEFINE_NUMERIC_OP_SINGLE(name, func) \
static OperationResult op_##name(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) { \
    OperationResult result = { \
        .type = RESULT_NUMERIC, \
        .error = false \
    }; \
    (void)op_arg; /* Ignore op_arg */ \
    (void)geom_b; /* Ignore geom_b */ \
    double value = func(geom_a); \
    if (value < 0 && !sfcgal_geometry_is_empty(geom_a)) { \
        result.error = true; \
        result.error_message = "Error calculating " #name; \
        return result; \
    } \
    result.numeric_result = value; \
    return result; \
}

/**
 * Macro for operations that return a numeric result (two geometry inputs)
 */
#define DEFINE_NUMERIC_OP_BINARY(name, func) \
static OperationResult op_##name(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) { \
    OperationResult result = { \
        .type = RESULT_NUMERIC, \
        .error = false \
    }; \
    (void)op_arg; /* Ignore op_arg */ \
    result.numeric_result = func(geom_a, geom_b); \
    return result; \
}

/**
 * Macro for operations that return a geometry (single geometry input)
 */
#define DEFINE_GEOMETRY_OP_SINGLE(name, func) \
static OperationResult op_##name(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) { \
    OperationResult result = { \
        .type = RESULT_GEOMETRY, \
        .error = false \
    }; \
    (void)op_arg; /* Ignore op_arg */ \
    (void)geom_b; /* Ignore geom_b */ \
    result.geometry_result = func(geom_a); \
    if (!result.geometry_result) { \
        result.error = true; \
        result.error_message = "Error calculating " #name; \
    } \
    return result; \
}

/**
 * Macro for operations that return a geometry (two geometry inputs)
 */
#define DEFINE_GEOMETRY_OP_BINARY(name, func) \
static OperationResult op_##name(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) { \
    OperationResult result = { \
        .type = RESULT_GEOMETRY, \
        .error = false \
    }; \
    (void)op_arg; /* Ignore op_arg */ \
    result.geometry_result = func(geom_a, geom_b); \
    if (!result.geometry_result) { \
        result.error = true; \
        result.error_message = "Error calculating " #name; \
    } \
    return result; \
}

/**
 * Macro for operations that modify a geometry in place and return a clone
 */
#define DEFINE_INPLACE_OP(name, func) \
static OperationResult op_##name(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) { \
    OperationResult result = { \
        .type = RESULT_GEOMETRY, \
        .error = false \
    }; \
    (void)op_arg; /* Ignore op_arg */ \
    (void)geom_b; /* Ignore geom_b */ \
    /* Clone geometry first to avoid modifying input */ \
    result.geometry_result = sfcgal_geometry_clone(geom_a); \
    if (!result.geometry_result) { \
        result.error = true; \
        result.error_message = "Failed to clone geometry"; \
        return result; \
    } \
    int changed = func(result.geometry_result); \
    if (changed < 0) { \
        result.error = true; \
        result.error_message = "Error executing " #name; \
        sfcgal_geometry_delete(result.geometry_result); \
        result.geometry_result = NULL; \
    } \
    return result; \
}

/* ========================================================================== */
/* DEFINE OPERATIONS USING MACROS */
/* ========================================================================== */

// Boolean operations - single geometry input
DEFINE_BOOLEAN_OP_SINGLE(is_valid, sfcgal_geometry_is_valid)
DEFINE_BOOLEAN_OP_SINGLE(is_planar, sfcgal_geometry_is_planar)

// Boolean operations - two geometry inputs
DEFINE_BOOLEAN_OP_BINARY(intersects, sfcgal_geometry_intersects)
DEFINE_BOOLEAN_OP_BINARY(intersects_3d, sfcgal_geometry_intersects_3d)
DEFINE_BOOLEAN_OP_BINARY(covers, sfcgal_geometry_covers)
DEFINE_BOOLEAN_OP_BINARY(covers_3d, sfcgal_geometry_covers_3d)
DEFINE_BOOLEAN_OP_BINARY(equals, sfcgal_geometry_is_equals)

// Numeric operations - single geometry input
DEFINE_NUMERIC_OP_SINGLE(volume, sfcgal_geometry_volume)
DEFINE_NUMERIC_OP_SINGLE(area, sfcgal_geometry_area)
DEFINE_NUMERIC_OP_SINGLE(area_3d, sfcgal_geometry_area_3d)
DEFINE_NUMERIC_OP_SINGLE(length, sfcgal_geometry_length)
DEFINE_NUMERIC_OP_SINGLE(length_3d, sfcgal_geometry_length_3d)
DEFINE_NUMERIC_OP_SINGLE(orientation, sfcgal_geometry_orientation)

// Numeric operations - two geometry inputs
DEFINE_NUMERIC_OP_BINARY(distance, sfcgal_geometry_distance)
DEFINE_NUMERIC_OP_BINARY(distance_3d, sfcgal_geometry_distance_3d)

// Geometry operations - single geometry input
DEFINE_GEOMETRY_OP_SINGLE(convexhull, sfcgal_geometry_convexhull)
DEFINE_GEOMETRY_OP_SINGLE(convexhull_3d, sfcgal_geometry_convexhull_3d)
DEFINE_GEOMETRY_OP_SINGLE(tesselate, sfcgal_geometry_tesselate)
DEFINE_GEOMETRY_OP_SINGLE(triangulate_2dz, sfcgal_geometry_triangulate_2dz)
DEFINE_GEOMETRY_OP_SINGLE(make_solid, sfcgal_geometry_make_solid)
DEFINE_GEOMETRY_OP_SINGLE(force_lhr, sfcgal_geometry_force_lhr)
DEFINE_GEOMETRY_OP_SINGLE(force_rhr, sfcgal_geometry_force_rhr)
DEFINE_GEOMETRY_OP_SINGLE(straight_skeleton, sfcgal_geometry_straight_skeleton)
DEFINE_GEOMETRY_OP_SINGLE(straight_skeleton_distance_in_m, sfcgal_geometry_straight_skeleton_distance_in_m)
DEFINE_GEOMETRY_OP_SINGLE(approximate_medial_axis, sfcgal_geometry_approximate_medial_axis)
DEFINE_GEOMETRY_OP_SINGLE(envelope, sfcgal_geometry_envelope)
DEFINE_GEOMETRY_OP_SINGLE(envelope_3d, sfcgal_geometry_envelope_3d)
DEFINE_GEOMETRY_OP_SINGLE(boundary, sfcgal_geometry_boundary)
DEFINE_GEOMETRY_OP_SINGLE(centroid, sfcgal_geometry_centroid)
DEFINE_GEOMETRY_OP_SINGLE(centroid_3d, sfcgal_geometry_centroid_3d)
DEFINE_GEOMETRY_OP_SINGLE(y_monotone_partition, sfcgal_y_monotone_partition_2)
DEFINE_GEOMETRY_OP_SINGLE(approx_convex_partition, sfcgal_approx_convex_partition_2)
DEFINE_GEOMETRY_OP_SINGLE(greene_approx_convex_partition, sfcgal_greene_approx_convex_partition_2)
DEFINE_GEOMETRY_OP_SINGLE(optimal_convex_partition, sfcgal_optimal_convex_partition_2)

// Geometry operations - two geometry inputs
DEFINE_GEOMETRY_OP_BINARY(intersection, sfcgal_geometry_intersection)
DEFINE_GEOMETRY_OP_BINARY(intersection_3d, sfcgal_geometry_intersection_3d)
DEFINE_GEOMETRY_OP_BINARY(difference, sfcgal_geometry_difference)
DEFINE_GEOMETRY_OP_BINARY(difference_3d, sfcgal_geometry_difference_3d)
DEFINE_GEOMETRY_OP_BINARY(union, sfcgal_geometry_union)
DEFINE_GEOMETRY_OP_BINARY(union_3d, sfcgal_geometry_union_3d)
DEFINE_GEOMETRY_OP_BINARY(minkowski_sum, sfcgal_geometry_minkowski_sum)
DEFINE_GEOMETRY_OP_BINARY(visibility_point, sfcgal_geometry_visibility_point)

// In-place operations
DEFINE_INPLACE_OP(drop_z, sfcgal_geometry_drop_z)
DEFINE_INPLACE_OP(drop_m, sfcgal_geometry_drop_m)

/* ========================================================================== */
/* CUSTOM OPERATION IMPLEMENTATIONS (for operations that don't fit the macros) */
/* ========================================================================== */

/**
 * Get detailed information about geometry validity
 */
static OperationResult op_is_complexity_detail(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_TEXT,
        .error = false
    };
    
    // Ignore op_arg and geom_b
    (void)op_arg;
    (void)geom_b;
    
    char* reason = NULL;
    sfcgal_geometry_t* location = NULL;
    
    int valid = sfcgal_geometry_is_valid_detail(geom_a, &reason, &location);
    
    if (valid < 0) {
        result.error = true;
        result.error_message = "Error checking validity detail";
        if (reason) free(reason);
        if (location) sfcgal_geometry_delete(location);
        return result;
    }
    
    if (valid == 1) {
        result.text_result = strdup("Geometry is valid");
    } else {
        if (reason != NULL) {
            char* full_message = malloc(strlen(reason) + 50);
            if (full_message) {
                sprintf(full_message, "Geometry is invalid. Reason: %s", reason);
                result.text_result = full_message;
            } else {
                result.error = true;
                result.error_message = "Memory allocation failed";
            }
            free(reason);
        } else {
            result.text_result = strdup("Geometry is invalid. No details available.");
        }
    }
    
    if (location) sfcgal_geometry_delete(location);
    
    return result;
}

/**
 * Get number of geometries in a collection
 */
static OperationResult op_num_geometries(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_NUMERIC,
        .error = false
    };
    
    // Ignore op_arg and geom_b
    (void)op_arg;
    (void)geom_b;
    
    size_t num = sfcgal_geometry_num_geometries(geom_a);
    result.numeric_result = (double)num;
    
    return result;
}

/**
 * Add Z coordinate to geometry with given default value
 */
static OperationResult op_force_z(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    // Ignore geom_b
    (void)geom_b;
    
    // Get default Z value from op_arg
    double defaultZ = 0.0;
    if (op_arg) {
        char* endptr;
        defaultZ = strtod(op_arg, &endptr);
        if (endptr == op_arg) {
            result.error = true;
            result.error_message = "Invalid default Z value";
            return result;
        }
    }
    
    // Clone geometry first to avoid modifying input
    result.geometry_result = sfcgal_geometry_clone(geom_a);
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Failed to clone geometry";
        return result;
    }
    
    int changed = sfcgal_geometry_force_z(result.geometry_result, defaultZ);
    
    if (changed < 0) {
        result.error = true;
        result.error_message = "Error forcing Z coordinate";
        sfcgal_geometry_delete(result.geometry_result);
        result.geometry_result = NULL;
    }
    
    return result;
}

/**
 * Add M coordinate to geometry with given default value
 */
static OperationResult op_force_m(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    // Ignore geom_b
    (void)geom_b;
    
    // Get default M value from op_arg
    double defaultM = 0.0;
    if (op_arg) {
        char* endptr;
        defaultM = strtod(op_arg, &endptr);
        if (endptr == op_arg) {
            result.error = true;
            result.error_message = "Invalid default M value";
            return result;
        }
    }
    
    // Clone geometry first to avoid modifying input
    result.geometry_result = sfcgal_geometry_clone(geom_a);
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Failed to clone geometry";
        return result;
    }
    
    int changed = sfcgal_geometry_force_m(result.geometry_result, defaultM);
    
    if (changed < 0) {
        result.error = true;
        result.error_message = "Error forcing M coordinate";
        sfcgal_geometry_delete(result.geometry_result);
        result.geometry_result = NULL;
    }
    
    return result;
}

/**
 * Almost equals with tolerance parameter
 */
static OperationResult op_almost_equals(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_BOOLEAN,
        .error = false
    };
    
    // Parse tolerance from op_arg
    double tolerance = 0.0;
    if (op_arg) {
        char* endptr;
        tolerance = strtod(op_arg, &endptr);
        if (endptr == op_arg) {
            result.error = true;
            result.error_message = "Invalid tolerance value";
            return result;
        }
    }
    
    int ret = sfcgal_geometry_is_almost_equals(geom_a, geom_b, tolerance);
    
    if (ret < 0) {
        result.error = true;
        result.error_message = "Error checking almost equals";
        return result;
    }
    
    result.boolean_result = (ret == 1);
    return result;
}

/**
 * Round operation with scale parameter
 */
static OperationResult op_round(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Parse scale from op_arg
    int scale = 0;
    if (op_arg) {
        char* endptr;
        long scale_long = strtol(op_arg, &endptr, 10);
        if (endptr == op_arg || scale_long < INT_MIN || scale_long > INT_MAX) {
            result.error = true;
            result.error_message = "Invalid scale value";
            return result;
        }
        scale = (int)scale_long;
    }
    
    result.geometry_result = sfcgal_geometry_round(geom_a, scale);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error rounding geometry";
    }
    
    return result;
}

/**
 * Extrude operation with parameters
 */
static OperationResult op_extrude(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Default values
    double dx = 0.0, dy = 0.0, dz = 1.0;
    
    // Parse parameters from op_arg
    if (op_arg) {
        // Format should be: "dx,dy,dz"
        char* args_copy = strdup(op_arg);
        if (!args_copy) {
            result.error = true;
            result.error_message = "Memory allocation failed";
            return result;
        }
        
        char* token = strtok(args_copy, ",");
        if (token) {
            dx = strtod(token, NULL);
            token = strtok(NULL, ",");
            if (token) {
                dy = strtod(token, NULL);
                token = strtok(NULL, ",");
                if (token) {
                    dz = strtod(token, NULL);
                }
            }
        }
        
        free(args_copy);
    }
    
    result.geometry_result = sfcgal_geometry_extrude(geom_a, dx, dy, dz);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error extruding geometry";
    }
    
    return result;
}

/**
 * Offset polygon operation with radius parameter
 */
static OperationResult op_offset_polygon(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Parse radius from op_arg
    double radius = 0.0;
    if (op_arg) {
        char* endptr;
        radius = strtod(op_arg, &endptr);
        if (endptr == op_arg) {
            result.error = true;
            result.error_message = "Invalid radius value";
            return result;
        }
    } else {
        result.error = true;
        result.error_message = "Offset polygon requires a radius parameter";
        return result;
    }
    
    result.geometry_result = sfcgal_geometry_offset_polygon(geom_a, radius);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error calculating offset polygon";
    }
    
    return result;
}

/**
 * Line substring operation with start and end parameters
 */
static OperationResult op_line_substring(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Default values
    double start = 0.0, end = 1.0;
    
    // Parse parameters from op_arg
    if (op_arg) {
        // Format should be: "start,end"
        char* args_copy = strdup(op_arg);
        if (!args_copy) {
            result.error = true;
            result.error_message = "Memory allocation failed";
            return result;
        }
        
        char* token = strtok(args_copy, ",");
        if (token) {
            start = strtod(token, NULL);
            token = strtok(NULL, ",");
            if (token) {
                end = strtod(token, NULL);
            }
        }
        
        free(args_copy);
    } else {
        result.error = true;
        result.error_message = "Line substring requires start and end parameters";
        return result;
    }
    
    result.geometry_result = sfcgal_geometry_line_sub_string(geom_a, start, end);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error calculating line substring";
    }
    
    return result;
}

/**
 * Extrude straight skeleton operation with height parameter
 */
static OperationResult op_extrude_straight_skeleton(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Parse height from op_arg
    double height = 1.0;
    if (op_arg) {
        char* endptr;
        height = strtod(op_arg, &endptr);
        if (endptr == op_arg) {
            result.error = true;
            result.error_message = "Invalid height value";
            return result;
        }
    }
    
    result.geometry_result = sfcgal_geometry_extrude_straight_skeleton(geom_a, height);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error extruding straight skeleton";
    }
    
    return result;
}

/**
 * Extrude polygon straight skeleton operation with building height and roof height parameters
 */
static OperationResult op_extrude_polygon_straight_skeleton(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Default values
    double building_height = 1.0, roof_height = 1.0;
    
    // Parse parameters from op_arg
    if (op_arg) {
        // Format should be: "building_height,roof_height"
        char* args_copy = strdup(op_arg);
        if (!args_copy) {
            result.error = true;
            result.error_message = "Memory allocation failed";
            return result;
        }
        
        char* token = strtok(args_copy, ",");
        if (token) {
            building_height = strtod(token, NULL);
            token = strtok(NULL, ",");
            if (token) {
                roof_height = strtod(token, NULL);
            }
        }
        
        free(args_copy);
    } else {
        result.error = true;
        result.error_message = "Extrude polygon straight skeleton requires building height and roof height parameters";
        return result;
    }
    
    result.geometry_result = sfcgal_geometry_extrude_polygon_straight_skeleton(geom_a, building_height, roof_height);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error extruding polygon straight skeleton";
    }
    
    return result;
}

/**
 * Straight skeleton partition operation with auto orientation parameter
 */
static OperationResult op_straight_skeleton_partition(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Parse auto orientation from op_arg
    bool auto_orientation = true;
    if (op_arg) {
        if (strcmp(op_arg, "0") == 0 || 
            strcmp(op_arg, "false") == 0 || 
            strcmp(op_arg, "False") == 0 || 
            strcmp(op_arg, "FALSE") == 0) {
            auto_orientation = false;
        }
    }
    
    result.geometry_result = sfcgal_geometry_straight_skeleton_partition(geom_a, auto_orientation);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error calculating straight skeleton partition";
    }
    
    return result;
}

#if !defined(_MSC_VER)
/**
 * Alpha shapes operation with alpha and allow_holes parameters
 */
static OperationResult op_alpha_shapes(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Default values
    double alpha = 0.0;
    bool allow_holes = false;
    
    // Parse parameters from op_arg
    if (op_arg) {
        // Format should be: "alpha,allow_holes"
        char* args_copy = strdup(op_arg);
        if (!args_copy) {
            result.error = true;
            result.error_message = "Memory allocation failed";
            return result;
        }
        
        char* token = strtok(args_copy, ",");
        if (token) {
            alpha = strtod(token, NULL);
            token = strtok(NULL, ",");
            if (token) {
                allow_holes = (strcmp(token, "1") == 0 || 
                              strcmp(token, "true") == 0 || 
                              strcmp(token, "True") == 0 || 
                              strcmp(token, "TRUE") == 0);
            }
        }
        
        free(args_copy);
    } else {
        result.error = true;
        result.error_message = "Alpha shapes requires alpha and allow_holes parameters";
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
static OperationResult op_optimal_alpha_shapes(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Default values
    bool allow_holes = false;
    size_t nb_components = 1;
    
    // Parse parameters from op_arg
    if (op_arg) {
        // Format should be: "allow_holes,nb_components"
        char* args_copy = strdup(op_arg);
        if (!args_copy) {
            result.error = true;
            result.error_message = "Memory allocation failed";
            return result;
        }
        
        char* token = strtok(args_copy, ",");
        if (token) {
            allow_holes = (strcmp(token, "1") == 0 || 
                          strcmp(token, "true") == 0 || 
                          strcmp(token, "True") == 0 || 
                          strcmp(token, "TRUE") == 0);
            token = strtok(NULL, ",");
            if (token) {
                long components = strtol(token, NULL, 10);
                if (components > 0) {
                    nb_components = (size_t)components;
                }
            }
        }
        
        free(args_copy);
    } else {
        result.error = true;
        result.error_message = "Optimal alpha shapes requires allow_holes and nb_components parameters";
        return result;
    }
    
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
static OperationResult op_alpha_wrapping_3d(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Default values
    size_t relative_alpha = 100;
    size_t relative_offset = 100;
    
    // Parse parameters from op_arg
    if (op_arg) {
        // Format should be: "relative_alpha,relative_offset"
        char* args_copy = strdup(op_arg);
        if (!args_copy) {
            result.error = true;
            result.error_message = "Memory allocation failed";
            return result;
        }
        
        char* token = strtok(args_copy, ",");
        if (token) {
            long alpha = strtol(token, NULL, 10);
            if (alpha > 0) {
                relative_alpha = (size_t)alpha;
            }
            token = strtok(NULL, ",");
            if (token) {
                long offset = strtol(token, NULL, 10);
                if (offset > 0) {
                    relative_offset = (size_t)offset;
                }
            }
        }
        
        free(args_copy);
    } else {
        result.error = true;
        result.error_message = "Alpha wrapping 3D requires relative_alpha and relative_offset parameters";
        return result;
    }
    
    result.geometry_result = sfcgal_geometry_alpha_wrapping_3d(geom_a, relative_alpha, relative_offset);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error calculating alpha wrapping 3D";
    }
    
    return result;
}

/**
 * Visibility segment operation with two point geometries 
 */
static OperationResult op_visibility_segment(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)op_arg; // Ignore op_arg
    
    // We need a third geometry (point B)
    sfcgal_geometry_t* point_b = NULL;
    
    // Parse the third geometry from op_arg
    if (op_arg && *op_arg) {
        point_b = sfcgal_io_read_wkt(op_arg, strlen(op_arg));
        if (!point_b) {
            result.error = true;
            result.error_message = "Failed to parse third geometry (point B)";
            return result;
        }
    } else {
        result.error = true;
        result.error_message = "Visibility segment requires a third geometry (point B)";
        return result;
    }
    
    result.geometry_result = sfcgal_geometry_visibility_segment(geom_a, geom_b, point_b);
    
    // Free the point we created
    sfcgal_geometry_delete(point_b);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error calculating visibility segment";
    }
    
    return result;
}

/**
 * Buffer 3D operation with radius, segments, and buffer_type parameters
 */
static OperationResult op_buffer3d(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Default values
    double radius = 1.0;
    int segments = 16;
    sfcgal_buffer3d_type_t buffer_type = SFCGAL_BUFFER3D_ROUND;
    
    // Parse parameters from op_arg
    if (op_arg) {
        // Format should be: "radius,segments,buffer_type"
        char* args_copy = strdup(op_arg);
        if (!args_copy) {
            result.error = true;
            result.error_message = "Memory allocation failed";
            return result;
        }
        
        char* token = strtok(args_copy, ",");
        if (token) {
            radius = strtod(token, NULL);
            token = strtok(NULL, ",");
            if (token) {
                long segs = strtol(token, NULL, 10);
                if (segs > 0) {
                    segments = (int)segs;
                }
                token = strtok(NULL, ",");
                if (token) {
                    if (strcmp(token, "round") == 0 || strcmp(token, "ROUND") == 0 || strcmp(token, "0") == 0) {
                        buffer_type = SFCGAL_BUFFER3D_ROUND;
                    } else if (strcmp(token, "cylsphere") == 0 || strcmp(token, "CYLSPHERE") == 0 || strcmp(token, "1") == 0) {
                        buffer_type = SFCGAL_BUFFER3D_CYLSPHERE;
                    } else if (strcmp(token, "flat") == 0 || strcmp(token, "FLAT") == 0 || strcmp(token, "2") == 0) {
                        buffer_type = SFCGAL_BUFFER3D_FLAT;
                    }
                }
            }
        }
        
        free(args_copy);
    } else {
        result.error = true;
        result.error_message = "Buffer3D requires radius, segments, and buffer_type parameters";
        return result;
    }
    
    result.geometry_result = sfcgal_geometry_buffer3d(geom_a, radius, segments, buffer_type);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error calculating buffer3D";
    }
    
    return result;
}

/**
 * Simplify operation with threshold and preserveTopology parameters
 */
static OperationResult op_simplify(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Default values
    double threshold = 1.0;
    bool preserve_topology = true;
    
    // Parse parameters from op_arg
    if (op_arg) {
        // Format should be: "threshold,preserve_topology"
        char* args_copy = strdup(op_arg);
        if (!args_copy) {
            result.error = true;
            result.error_message = "Memory allocation failed";
            return result;
        }
        
        char* token = strtok(args_copy, ",");
        if (token) {
            threshold = strtod(token, NULL);
            token = strtok(NULL, ",");
            if (token) {
                preserve_topology = (strcmp(token, "1") == 0 || 
                                   strcmp(token, "true") == 0 || 
                                   strcmp(token, "True") == 0 || 
                                   strcmp(token, "TRUE") == 0);
            }
        }
        
        free(args_copy);
    } else {
        result.error = true;
        result.error_message = "Simplify requires threshold and preserve_topology parameters";
        return result;
    }
    
    result.geometry_result = sfcgal_geometry_simplify(geom_a, threshold, preserve_topology);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error simplifying geometry";
    }
    
    return result;
}

/**
 * Rotate operation with angle parameter
 */
static OperationResult op_rotate(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Parse angle from op_arg
    double angle = 0.0;
    if (op_arg) {
        char* endptr;
        angle = strtod(op_arg, &endptr);
        if (endptr == op_arg) {
            result.error = true;
            result.error_message = "Invalid angle value";
            return result;
        }
    } else {
        result.error = true;
        result.error_message = "Rotate requires an angle parameter";
        return result;
    }
    
    result.geometry_result = sfcgal_geometry_rotate(geom_a, angle);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error rotating geometry";
    }
    
    return result;
}

/**
 * Rotate 2D operation with angle, cx, cy parameters
 */
static OperationResult op_rotate_2d(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Default values
    double angle = 0.0, cx = 0.0, cy = 0.0;
    
    // Parse parameters from op_arg
    if (op_arg) {
        // Format should be: "angle,cx,cy"
        char* args_copy = strdup(op_arg);
        if (!args_copy) {
            result.error = true;
            result.error_message = "Memory allocation failed";
            return result;
        }
        
        char* token = strtok(args_copy, ",");
        if (token) {
            angle = strtod(token, NULL);
            token = strtok(NULL, ",");
            if (token) {
                cx = strtod(token, NULL);
                token = strtok(NULL, ",");
                if (token) {
                    cy = strtod(token, NULL);
                }
            }
        }
        
        free(args_copy);
    } else {
        result.error = true;
        result.error_message = "Rotate 2D requires angle, cx, cy parameters";
        return result;
    }
    
    result.geometry_result = sfcgal_geometry_rotate_2d(geom_a, angle, cx, cy);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error rotating geometry in 2D";
    }
    
    return result;
}

/**
 * Rotate 3D operation with angle, ax, ay, az parameters
 */
static OperationResult op_rotate_3d(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Default values
    double angle = 0.0, ax = 0.0, ay = 0.0, az = 1.0;
    
    // Parse parameters from op_arg
    if (op_arg) {
        // Format should be: "angle,ax,ay,az"
        char* args_copy = strdup(op_arg);
        if (!args_copy) {
            result.error = true;
            result.error_message = "Memory allocation failed";
            return result;
        }
        
        char* token = strtok(args_copy, ",");
        if (token) {
            angle = strtod(token, NULL);
            token = strtok(NULL, ",");
            if (token) {
                ax = strtod(token, NULL);
                token = strtok(NULL, ",");
                if (token) {
                    ay = strtod(token, NULL);
                    token = strtok(NULL, ",");
                    if (token) {
                        az = strtod(token, NULL);
                    }
                }
            }
        }
        
        free(args_copy);
    } else {
        result.error = true;
        result.error_message = "Rotate 3D requires angle, ax, ay, az parameters";
        return result;
    }
    
    result.geometry_result = sfcgal_geometry_rotate_3d(geom_a, angle, ax, ay, az);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error rotating geometry in 3D";
    }
    
    return result;
}

/**
 * Rotate 3D around center operation with angle, ax, ay, az, cx, cy, cz parameters
 */
static OperationResult op_rotate_3d_around_center(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Default values
    double angle = 0.0, ax = 0.0, ay = 0.0, az = 1.0, cx = 0.0, cy = 0.0, cz = 0.0;
    
    // Parse parameters from op_arg
    if (op_arg) {
        // Format should be: "angle,ax,ay,az,cx,cy,cz"
        char* args_copy = strdup(op_arg);
        if (!args_copy) {
            result.error = true;
            result.error_message = "Memory allocation failed";
            return result;
        }
        
        char* token = strtok(args_copy, ",");
        if (token) {
            angle = strtod(token, NULL);
            token = strtok(NULL, ",");
            if (token) {
                ax = strtod(token, NULL);
                token = strtok(NULL, ",");
                if (token) {
                    ay = strtod(token, NULL);
                    token = strtok(NULL, ",");
                    if (token) {
                        az = strtod(token, NULL);
                        token = strtok(NULL, ",");
                        if (token) {
                            cx = strtod(token, NULL);
                            token = strtok(NULL, ",");
                            if (token) {
                                cy = strtod(token, NULL);
                                token = strtok(NULL, ",");
                                if (token) {
                                    cz = strtod(token, NULL);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        free(args_copy);
    } else {
        result.error = true;
        result.error_message = "Rotate 3D around center requires angle, ax, ay, az, cx, cy, cz parameters";
        return result;
    }
    
    result.geometry_result = sfcgal_geometry_rotate_3d_around_center(geom_a, angle, ax, ay, az, cx, cy, cz);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error rotating geometry in 3D around center";
    }
    
    return result;
}

/**
 * Rotate X operation with angle parameter
 */
static OperationResult op_rotate_x(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Parse angle from op_arg
    double angle = 0.0;
    if (op_arg) {
        char* endptr;
        angle = strtod(op_arg, &endptr);
        if (endptr == op_arg) {
            result.error = true;
            result.error_message = "Invalid angle value";
            return result;
        }
    } else {
        result.error = true;
        result.error_message = "Rotate X requires an angle parameter";
        return result;
    }
    
    result.geometry_result = sfcgal_geometry_rotate_x(geom_a, angle);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error rotating geometry around X axis";
    }
    
    return result;
}

/**
 * Rotate Y operation with angle parameter
 */
static OperationResult op_rotate_y(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Parse angle from op_arg
    double angle = 0.0;
    if (op_arg) {
        char* endptr;
        angle = strtod(op_arg, &endptr);
        if (endptr == op_arg) {
            result.error = true;
            result.error_message = "Invalid angle value";
            return result;
        }
    } else {
        result.error = true;
        result.error_message = "Rotate Y requires an angle parameter";
        return result;
    }
    
    result.geometry_result = sfcgal_geometry_rotate_y(geom_a, angle);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error rotating geometry around Y axis";
    }
    
    return result;
}

/**
 * Rotate Z operation with angle parameter
 */
static OperationResult op_rotate_z(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Parse angle from op_arg
    double angle = 0.0;
    if (op_arg) {
        char* endptr;
        angle = strtod(op_arg, &endptr);
        if (endptr == op_arg) {
            result.error = true;
            result.error_message = "Invalid angle value";
            return result;
        }
    } else {
        result.error = true;
        result.error_message = "Rotate Z requires an angle parameter";
        return result;
    }
    
    result.geometry_result = sfcgal_geometry_rotate_z(geom_a, angle);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error rotating geometry around Z axis";
    }
    
    return result;
}

/**
 * Scale operation with scale parameter
 */
static OperationResult op_scale(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Parse scale from op_arg
    double scale = 1.0;
    if (op_arg) {
        char* endptr;
        scale = strtod(op_arg, &endptr);
        if (endptr == op_arg) {
            result.error = true;
            result.error_message = "Invalid scale value";
            return result;
        }
    } else {
        result.error = true;
        result.error_message = "Scale requires a scale parameter";
        return result;
    }
    
    result.geometry_result = sfcgal_geometry_scale(geom_a, scale);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error scaling geometry";
    }
    
    return result;
}

/**
 * Scale 3D operation with sx, sy, sz parameters
 */
static OperationResult op_scale_3d(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Default values
    double sx = 1.0, sy = 1.0, sz = 1.0;
    
    // Parse parameters from op_arg
    if (op_arg) {
        // Format should be: "sx,sy,sz"
        char* args_copy = strdup(op_arg);
        if (!args_copy) {
            result.error = true;
            result.error_message = "Memory allocation failed";
            return result;
        }
        
        char* token = strtok(args_copy, ",");
        if (token) {
            sx = strtod(token, NULL);
            token = strtok(NULL, ",");
            if (token) {
                sy = strtod(token, NULL);
                token = strtok(NULL, ",");
                if (token) {
                    sz = strtod(token, NULL);
                }
            }
        }
        
        free(args_copy);
    } else {
        result.error = true;
        result.error_message = "Scale 3D requires sx, sy, sz parameters";
        return result;
    }
    
    result.geometry_result = sfcgal_geometry_scale_3d(geom_a, sx, sy, sz);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error scaling geometry in 3D";
    }
    
    return result;
}

/**
 * Scale 3D around center operation with sx, sy, sz, cx, cy, cz parameters
 */
static OperationResult op_scale_3d_around_center(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Default values
    double sx = 1.0, sy = 1.0, sz = 1.0, cx = 0.0, cy = 0.0, cz = 0.0;
    
    // Parse parameters from op_arg
    if (op_arg) {
        // Format should be: "sx,sy,sz,cx,cy,cz"
        char* args_copy = strdup(op_arg);
        if (!args_copy) {
            result.error = true;
            result.error_message = "Memory allocation failed";
            return result;
        }
        
        char* token = strtok(args_copy, ",");
        if (token) {
            sx = strtod(token, NULL);
            token = strtok(NULL, ",");
            if (token) {
                sy = strtod(token, NULL);
                token = strtok(NULL, ",");
                if (token) {
                    sz = strtod(token, NULL);
                    token = strtok(NULL, ",");
                    if (token) {
                        cx = strtod(token, NULL);
                        token = strtok(NULL, ",");
                        if (token) {
                            cy = strtod(token, NULL);
                            token = strtok(NULL, ",");
                            if (token) {
                                cz = strtod(token, NULL);
                            }
                        }
                    }
                }
            }
        }
        
        free(args_copy);
    } else {
        result.error = true;
        result.error_message = "Scale 3D around center requires sx, sy, sz, cx, cy, cz parameters";
        return result;
    }
    
    result.geometry_result = sfcgal_geometry_scale_3d_around_center(geom_a, sx, sy, sz, cx, cy, cz);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error scaling geometry in 3D around center";
    }
    
    return result;
}

/**
 * Translate 2D operation with dx, dy parameters
 */
static OperationResult op_translate_2d(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Default values
    double dx = 0.0, dy = 0.0;
    
    // Parse parameters from op_arg
    if (op_arg) {
        // Format should be: "dx,dy"
        char* args_copy = strdup(op_arg);
        if (!args_copy) {
            result.error = true;
            result.error_message = "Memory allocation failed";
            return result;
        }
        
        char* token = strtok(args_copy, ",");
        if (token) {
            dx = strtod(token, NULL);
            token = strtok(NULL, ",");
            if (token) {
                dy = strtod(token, NULL);
            }
        }
        
        free(args_copy);
    } else {
        result.error = true;
        result.error_message = "Translate 2D requires dx, dy parameters";
        return result;
    }
    
    result.geometry_result = sfcgal_geometry_translate_2d(geom_a, dx, dy);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error translating geometry in 2D";
    }
    
    return result;
}

/**
 * Translate 3D operation with dx, dy, dz parameters
 */
static OperationResult op_translate_3d(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_GEOMETRY,
        .error = false
    };
    
    (void)geom_b; // Ignore geom_b
    
    // Default values
    double dx = 0.0, dy = 0.0, dz = 0.0;
    
    // Parse parameters from op_arg
    if (op_arg) {
        // Format should be: "dx,dy,dz"
        char* args_copy = strdup(op_arg);
        if (!args_copy) {
            result.error = true;
            result.error_message = "Memory allocation failed";
            return result;
        }
        
        char* token = strtok(args_copy, ",");
        if (token) {
            dx = strtod(token, NULL);
            token = strtok(NULL, ",");
            if (token) {
                dy = strtod(token, NULL);
                token = strtok(NULL, ",");
                if (token) {
                    dz = strtod(token, NULL);
                }
            }
        }
        
        free(args_copy);
    } else {
        result.error = true;
        result.error_message = "Translate 3D requires dx, dy, dz parameters";
        return result;
    }
    
    result.geometry_result = sfcgal_geometry_translate_3d(geom_a, dx, dy, dz);
    
    if (!result.geometry_result) {
        result.error = true;
        result.error_message = "Error translating geometry in 3D";
    }
    
    return result;
}

/* ========================================================================== */
/* OPERATIONS REGISTRY */
/* ========================================================================== */

// Array of supported operations
static const Operation operations[] = {
    /* Validity operations */
    {"is_valid", op_is_valid, false, "Test if geometry is valid", "Validity", RESULT_BOOLEAN, false},
    {"is_complexity_detail", op_is_complexity_detail, false, "Get detailed information about geometry validity", "Validity", RESULT_TEXT, false},
    {"is_planar", op_is_planar, false, "Test if polygon is planar", "Validity", RESULT_BOOLEAN, false},
    
    /* Coordinate handling operations */
    {"drop_z", op_drop_z, false, "Remove Z coordinate from geometry", "Coordinate Handling", RESULT_GEOMETRY, false},
    {"drop_m", op_drop_m, false, "Remove M coordinate from geometry", "Coordinate Handling", RESULT_GEOMETRY, false},
    {"force_z", op_force_z, false, "Add Z coordinate with given default value", "Coordinate Handling", RESULT_GEOMETRY, true},
    {"force_m", op_force_m, false, "Add M coordinate with given default value", "Coordinate Handling", RESULT_GEOMETRY, true},
    
    /* Information operations */
    {"num_geometries", op_num_geometries, false, "Get number of geometries in a collection", "Information", RESULT_NUMERIC, false},
    {"orientation", op_orientation, false, "Get orientation of polygon (-1=CCW, 1=CW, 0=invalid)", "Information", RESULT_NUMERIC, false},
    
    /* Spatial analysis operations */
    {"intersects", op_intersects, true, "Check if geometries intersect", "Spatial Analysis", RESULT_BOOLEAN, false},
    {"intersects_3d", op_intersects_3d, true, "Check if geometries intersect in 3D", "Spatial Analysis", RESULT_BOOLEAN, false},
    {"covers", op_covers, true, "Check if geometry A covers geometry B", "Spatial Analysis", RESULT_BOOLEAN, false},
    {"covers_3d", op_covers_3d, true, "Check if geometry A covers geometry B in 3D", "Spatial Analysis", RESULT_BOOLEAN, false},
    {"equals", op_equals, true, "Check if geometries are equal", "Spatial Analysis", RESULT_BOOLEAN, false},
    {"almost_equals", op_almost_equals, true, "Check if geometries are almost equal (with tolerance)", "Spatial Analysis", RESULT_BOOLEAN, true},
    
    /* Construction operations */
    {"intersection", op_intersection, true, "Get intersection of geometries", "Construction", RESULT_GEOMETRY, false},
    {"intersection_3d", op_intersection_3d, true, "Get 3D intersection of geometries", "Construction", RESULT_GEOMETRY, false},
    {"difference", op_difference, true, "Get difference of geometries", "Construction", RESULT_GEOMETRY, false},
    {"difference_3d", op_difference_3d, true, "Get 3D difference of geometries", "Construction", RESULT_GEOMETRY, false},
    {"union", op_union, true, "Get union of geometries", "Construction", RESULT_GEOMETRY, false},
    {"union_3d", op_union_3d, true, "Get 3D union of geometries", "Construction", RESULT_GEOMETRY, false},
    {"convexhull", op_convexhull, false, "Get convex hull of geometry", "Construction", RESULT_GEOMETRY, false},
    {"convexhull_3d", op_convexhull_3d, false, "Get 3D convex hull of geometry", "Construction", RESULT_GEOMETRY, false},
    {"make_solid", op_make_solid, false, "Convert polyhedral surface to solid", "Construction", RESULT_GEOMETRY, false},
    {"boundary", op_boundary, false, "Get boundary of geometry", "Construction", RESULT_GEOMETRY, false},
    {"centroid", op_centroid, false, "Get centroid of geometry", "Construction", RESULT_GEOMETRY, false},
    {"centroid_3d", op_centroid_3d, false, "Get 3D centroid of geometry", "Construction", RESULT_GEOMETRY, false},
    {"envelope", op_envelope, false, "Get 2D envelope (bounding box) of geometry", "Construction", RESULT_GEOMETRY, false},
    {"envelope_3d", op_envelope_3d, false, "Get 3D envelope (bounding box) of geometry", "Construction", RESULT_GEOMETRY, false},
    {"minkowski_sum", op_minkowski_sum, true, "Get Minkowski sum of geometries", "Construction", RESULT_GEOMETRY, false},
    
    /* Tesselation operations */
    {"tesselate", op_tesselate, false, "Tesselate geometry into triangles", "Tesselation", RESULT_GEOMETRY, false},
    {"triangulate_2dz", op_triangulate_2dz, false, "Triangulate 2D geometry preserving Z", "Tesselation", RESULT_GEOMETRY, false},
    {"y_monotone_partition", op_y_monotone_partition, false, "Partition into y-monotone polygons", "Tesselation", RESULT_GEOMETRY, false},
    {"approx_convex_partition", op_approx_convex_partition, false, "Partition into approx. convex polygons", "Tesselation", RESULT_GEOMETRY, false},
    {"greene_approx_convex_partition", op_greene_approx_convex_partition, false, "Greene's approx. convex partition", "Tesselation", RESULT_GEOMETRY, false},
    {"optimal_convex_partition", op_optimal_convex_partition, false, "Optimal convex partition", "Tesselation", RESULT_GEOMETRY, false},
    
    /* Orientation operations */
    {"force_lhr", op_force_lhr, false, "Force left-hand rule orientation", "Orientation", RESULT_GEOMETRY, false},
    {"force_rhr", op_force_rhr, false, "Force right-hand rule orientation", "Orientation", RESULT_GEOMETRY, false},
    
    /* Metrics operations */
    {"distance", op_distance, true, "Calculate distance between geometries", "Metrics", RESULT_NUMERIC, false},
    {"distance_3d", op_distance_3d, true, "Calculate 3D distance between geometries", "Metrics", RESULT_NUMERIC, false},
    {"volume", op_volume, false, "Calculate volume of 3D geometry", "Metrics", RESULT_NUMERIC, false},
    {"area", op_area, false, "Calculate area of 2D geometry", "Metrics", RESULT_NUMERIC, false},
    {"area_3d", op_area_3d, false, "Calculate area of 3D geometry", "Metrics", RESULT_NUMERIC, false},
    {"length", op_length, false, "Calculate length of geometry", "Metrics", RESULT_NUMERIC, false},
    {"length_3d", op_length_3d, false, "Calculate 3D length of geometry", "Metrics", RESULT_NUMERIC, false},
    
    /* Advanced operations */
    {"extrude", op_extrude, false, "Extrude 2D geometry to 3D", "Advanced", RESULT_GEOMETRY, true},
    {"buffer3d", op_buffer3d, false, "Create 3D buffer around geometry", "Advanced", RESULT_GEOMETRY, true},
    {"round", op_round, false, "Round coordinates to given precision", "Advanced", RESULT_GEOMETRY, true},
    {"offset_polygon", op_offset_polygon, false, "Create offset polygon", "Advanced", RESULT_GEOMETRY, true},
    {"line_substring", op_line_substring, false, "Get substring of linestring", "Advanced", RESULT_GEOMETRY, true},
    {"simplify", op_simplify, false, "Simplify geometry", "Advanced", RESULT_GEOMETRY, true},
    
    /* Skeleton operations */
    {"straight_skeleton", op_straight_skeleton, false, "Compute straight skeleton", "Skeleton", RESULT_GEOMETRY, false},
    {"straight_skeleton_distance_in_m", op_straight_skeleton_distance_in_m, false, "Compute straight skeleton with distances", "Skeleton", RESULT_GEOMETRY, false},
    {"extrude_straight_skeleton", op_extrude_straight_skeleton, false, "Extrude using straight skeleton", "Skeleton", RESULT_GEOMETRY, true},
    {"extrude_polygon_straight_skeleton", op_extrude_polygon_straight_skeleton, false, "Extrude polygon with roof", "Skeleton", RESULT_GEOMETRY, true},
    {"approximate_medial_axis", op_approximate_medial_axis, false, "Compute approximate medial axis", "Skeleton", RESULT_GEOMETRY, false},
    {"straight_skeleton_partition", op_straight_skeleton_partition, false, "Partition using straight skeleton", "Skeleton", RESULT_GEOMETRY, true},
    
    /* Alpha shape operations */
#if !defined(_MSC_VER)
    {"alpha_shapes", op_alpha_shapes, false, "Compute alpha shapes", "Alpha Shapes", RESULT_GEOMETRY, true},
    {"optimal_alpha_shapes", op_optimal_alpha_shapes, false, "Compute optimal alpha shapes", "Alpha Shapes", RESULT_GEOMETRY, true},
#endif
    {"alpha_wrapping_3d", op_alpha_wrapping_3d, false, "Compute 3D alpha wrapping", "Alpha Shapes", RESULT_GEOMETRY, true},
    
    /* Visibility operations */
    {"visibility_point", op_visibility_point, true, "Compute visibility polygon from point", "Visibility", RESULT_GEOMETRY, false},
    {"visibility_segment", op_visibility_segment, true, "Compute visibility polygon from segment", "Visibility", RESULT_GEOMETRY, true},
    
    /* Transformation operations */
    {"rotate", op_rotate, false, "Rotate geometry", "Transformation", RESULT_GEOMETRY, true},
    {"rotate_2d", op_rotate_2d, false, "Rotate geometry in 2D around point", "Transformation", RESULT_GEOMETRY, true},
    {"rotate_3d", op_rotate_3d, false, "Rotate geometry in 3D around axis", "Transformation", RESULT_GEOMETRY, true},
    {"rotate_3d_around_center", op_rotate_3d_around_center, false, "Rotate geometry in 3D around axis and center", "Transformation", RESULT_GEOMETRY, true},
    {"rotate_x", op_rotate_x, false, "Rotate geometry around X axis", "Transformation", RESULT_GEOMETRY, true},
    {"rotate_y", op_rotate_y, false, "Rotate geometry around Y axis", "Transformation", RESULT_GEOMETRY, true},
    {"rotate_z", op_rotate_z, false, "Rotate geometry around Z axis", "Transformation", RESULT_GEOMETRY, true},
    {"scale", op_scale, false, "Scale geometry uniformly", "Transformation", RESULT_GEOMETRY, true},
    {"scale_3d", op_scale_3d, false, "Scale geometry in 3D", "Transformation", RESULT_GEOMETRY, true},
    {"scale_3d_around_center", op_scale_3d_around_center, false, "Scale geometry in 3D around center", "Transformation", RESULT_GEOMETRY, true},
    {"translate_2d", op_translate_2d, false, "Translate geometry in 2D", "Transformation", RESULT_GEOMETRY, true},
    {"translate_3d", op_translate_3d, false, "Translate geometry in 3D", "Transformation", RESULT_GEOMETRY, true},
    
    {NULL, NULL, false, NULL, NULL, RESULT_BOOLEAN, false} // sentinel value
};

/**
 * Find an operation by name
 */
const Operation* find_operation(const char* name) {
    for (int i = 0; operations[i].name != NULL; i++) {
        if (strcmp(operations[i].name, name) == 0) {
            return &operations[i];
        }
    }
    return NULL;
}

/**
 * Helper function to get a string representation of the input type for an operation
 */
static const char* get_input_type_str(const Operation* op) {
    static char buffer[10];
    
    if (op->requires_geom_b) {
        if (op->requires_arg) {
            strcpy(buffer, "A, B, N");
        } else {
            strcpy(buffer, "A, B");
        }
    } else {
        if (op->requires_arg) {
            strcpy(buffer, "A, N");
        } else {
            strcpy(buffer, "A");
        }
    }
    
    return buffer;
}

/**
 * Helper function to get a string representation of the output type for an operation
 */
static const char* get_output_type_str(const Operation* op) {
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
 * Execute an SFCGAL operation
 */
OperationResult execute_operation(const char* op_name, const char* op_arg, 
                                 sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
    OperationResult result = {
        .type = RESULT_BOOLEAN,
        .boolean_result = false,
        .error = true,
        .error_message = "Unknown operation"
    };
    
    // Find operation
    const Operation* op = find_operation(op_name);
    if (op == NULL) {
        return result;
    }
    
    // Check if operation requires geom_b
    if (op->requires_geom_b && geom_b == NULL) {
        result.error_message = "Operation requires a second geometry";
        return result;
    }
    
    // Execute operation
    return op->func(op_arg, geom_a, geom_b);
}

/**
 * Get array of unique categories
 */
static char** get_unique_categories(int* num_categories) {
    // First, count unique categories
    int max_categories = 0;
    for (int i = 0; operations[i].name != NULL; i++) {
        max_categories++;
    }
    
    // We know max_categories is positive since there are operations defined
    char** categories = (char**)malloc((size_t)max_categories * sizeof(char*));
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
            categories[*num_categories] = (char*)operations[i].category;
            (*num_categories)++;
        }
    }
    
    return categories;
}

/**
 * Calculate the maximum widths for each column in the operations table
 * 
 * @param max_category_width Pointer to store the maximum category width
 * @param max_operation_width Pointer to store the maximum operation name width
 * @param max_input_width Pointer to store the maximum input type width
 * @param max_output_width Pointer to store the maximum output type width
 * @param max_description_width Pointer to store the maximum description width
 */
static void calculate_column_widths(size_t* max_category_width, 
                                   size_t* max_operation_width,
                                   size_t* max_input_width,
                                   size_t* max_output_width,
                                   size_t* max_description_width) {
    // Initialize with minimum widths
    *max_category_width = 8;      // "Category"
    *max_operation_width = 9;     // "Operation"
    *max_input_width = 5;         // "Input"
    *max_output_width = 6;        // "Output"
    *max_description_width = 11;  // "Description"
    
    // Get unique categories
    int num_categories = 0;
    char** categories = get_unique_categories(&num_categories);
    
    if (!categories) {
        return;
    }
    
    // Check category width
    for (int i = 0; i < num_categories; i++) {
        size_t category_len = strlen(categories[i]);
        if (category_len > *max_category_width) {
            *max_category_width = category_len;
        }
    }
    
    // Free categories now that we're done with them
    free(categories);
    
    // Check operation width, input width, output width, and description width
    for (int i = 0; operations[i].name != NULL; i++) {
        // Operation name
        size_t name_len = strlen(operations[i].name);
        if (name_len > *max_operation_width) {
            *max_operation_width = name_len;
        }
        
        // Input type
        const char* input_type = get_input_type_str(&operations[i]);
        size_t input_len = strlen(input_type);
        if (input_len > *max_input_width) {
            *max_input_width = input_len;
        }
        
        // Output type (single character, no need to check)
        
        // Description
        size_t desc_len = strlen(operations[i].description);
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
 * Print available operations to stdout in a well-formatted table, dynamically generated
 * from the operations array
 */
void print_available_operations(void) {
    // Calculate column widths
    size_t cat_width, op_width, in_width, out_width, desc_width;
    calculate_column_widths(&cat_width, &op_width, &in_width, &out_width, &desc_width);
    
    // Calculate total width for horizontal lines
    size_t total_width = cat_width + op_width + in_width + out_width + desc_width + 6; // 6 for the vertical bars
    
    // Print top border
    printf("┌");
    for (size_t i = 0; i < cat_width; i++) printf("─");
    printf("┬");
    for (size_t i = 0; i < op_width; i++) printf("─");
    printf("┬");
    for (size_t i = 0; i < in_width; i++) printf("─");
    printf("┬");
    for (size_t i = 0; i < out_width; i++) printf("─");
    printf("┬");
    for (size_t i = 0; i < desc_width; i++) printf("─");
    printf("┐\n");
    
    // Print header
    printf("│ %-*s│ %-*s│ %-*s│ %-*s│ %-*s│\n",
           (int)cat_width-1, "Category",
           (int)op_width-1, "Operation",
           (int)in_width-1, "Input",
           (int)out_width-1, "Output",
           (int)desc_width-1, "Description");
    
    // Print header separator
    printf("├");
    for (size_t i = 0; i < cat_width; i++) printf("─");
    printf("┼");
    for (size_t i = 0; i < op_width; i++) printf("─");
    printf("┼");
    for (size_t i = 0; i < in_width; i++) printf("─");
    printf("┼");
    for (size_t i = 0; i < out_width; i++) printf("─");
    printf("┼");
    for (size_t i = 0; i < desc_width; i++) printf("─");
    printf("┤\n");
    
    // Get unique categories
    int num_categories = 0;
    char** categories = get_unique_categories(&num_categories);
    
    if (!categories) {
        printf("│ Error generating categories table%*s│\n", (int)(total_width - 34), "");
        printf("└");
        for (size_t i = 0; i < total_width - 2; i++) printf("─");
        printf("┘\n");
        return;
    }
    
    // Print operations by category
    for (int cat_idx = 0; cat_idx < num_categories; cat_idx++) {
        const char* category = categories[cat_idx];
        int printed = 0;
        
        // Print operations in this category
        for (int i = 0; operations[i].name != NULL; i++) {
            if (strcmp(operations[i].category, category) == 0) {
                if (printed == 0) {
                    // First operation in category - print category name
                    printf("│ %-*s│ %-*s│ %-*s│ %-*s│ %-*s│\n",
                           (int)cat_width-1, category,
                           (int)op_width-1, operations[i].name,
                           (int)in_width-1, get_input_type_str(&operations[i]),
                           (int)out_width-1, get_output_type_str(&operations[i]),
                           (int)desc_width-1, operations[i].description);
                } else {
                    // Subsequent operations - leave category column blank
                    printf("│ %-*s│ %-*s│ %-*s│ %-*s│ %-*s│\n",
                           (int)cat_width-1, "",
                           (int)op_width-1, operations[i].name,
                           (int)in_width-1, get_input_type_str(&operations[i]),
                           (int)out_width-1, get_output_type_str(&operations[i]),
                           (int)desc_width-1, operations[i].description);
                }
                printed++;
            }
        }
        
        // Print separator between categories, unless it's the last category
        if (cat_idx < num_categories - 1) {
            printf("├");
            for (size_t i = 0; i < cat_width; i++) printf("─");
            printf("┼");
            for (size_t i = 0; i < op_width; i++) printf("─");
            printf("┼");
            for (size_t i = 0; i < in_width; i++) printf("─");
            printf("┼");
            for (size_t i = 0; i < out_width; i++) printf("─");
            printf("┼");
            for (size_t i = 0; i < desc_width; i++) printf("─");
            printf("┤\n");
        }
    }
    
    // Print table footer
    printf("└");
    for (size_t i = 0; i < cat_width; i++) printf("─");
    printf("┴");
    for (size_t i = 0; i < op_width; i++) printf("─");
    printf("┴");
    for (size_t i = 0; i < in_width; i++) printf("─");
    printf("┴");
    for (size_t i = 0; i < out_width; i++) printf("─");
    printf("┴");
    for (size_t i = 0; i < desc_width; i++) printf("─");
    printf("┘\n");
    
    // Legend
    printf("\nLegend:\n");
    printf("  A, B: Input geometries      G: Output geometry      D: Numeric result\n");
    printf("  T: Text result              N: Numeric argument\n");
    
    // Free allocated memory
    free(categories);
}
