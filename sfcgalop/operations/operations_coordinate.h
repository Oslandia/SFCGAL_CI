/**
 * operations_coordinate.h - Coordinate manipulation operations
 */

#ifndef OPERATIONS_COORDINATE_H
#define OPERATIONS_COORDINATE_H

#include "operations.h"
#include "operations_params.h"

/**
 * Remove Z coordinate from geometry
 */
OperationResult op_drop_z(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Remove M coordinate from geometry
 */
OperationResult op_drop_m(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Add Z coordinate with given default value
 */
OperationResult op_force_z(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Add M coordinate with given default value
 */
OperationResult op_force_m(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Force left-hand rule orientation
 */
OperationResult op_force_lhr(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Force right-hand rule orientation
 */
OperationResult op_force_rhr(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

#endif /* OPERATIONS_COORDINATE_H */
