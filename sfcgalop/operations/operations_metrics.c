/**
 * operations_metrics.c - Implementation of geometry metrics operations
 */

#include "operations_metrics.h"
#include "operations_common.h"
#include "../util.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/**
 * Numeric operations - single geometry input
 */
DEFINE_NUMERIC_OP_SINGLE(volume, sfcgal_geometry_volume)
DEFINE_NUMERIC_OP_SINGLE(area, sfcgal_geometry_area)
DEFINE_NUMERIC_OP_SINGLE(area_3d, sfcgal_geometry_area_3d)
DEFINE_NUMERIC_OP_SINGLE(length, sfcgal_geometry_length)
DEFINE_NUMERIC_OP_SINGLE(length_3d, sfcgal_geometry_length_3d)
DEFINE_NUMERIC_OP_SINGLE(orientation, sfcgal_geometry_orientation)

/**
 * Numeric operations - two geometry inputs
 */
DEFINE_NUMERIC_OP_BINARY(distance, sfcgal_geometry_distance)
DEFINE_NUMERIC_OP_BINARY(distance_3d, sfcgal_geometry_distance_3d)

/**
 * Get number of geometries in a collection
 */
OperationResult op_num_geometries(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b) {
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
