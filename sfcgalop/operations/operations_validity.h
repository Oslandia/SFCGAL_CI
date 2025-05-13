/**
 * operations_validity.h - Geometry validity operations
 */

#ifndef OPERATIONS_VALIDITY_H
#define OPERATIONS_VALIDITY_H

#include "operations.h"
#include "operations_params.h"

/**
 * Check if a geometry is valid
 */
OperationResult op_is_valid(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Get detailed information about geometry validity
 */
OperationResult op_is_validity_detail(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Check if a polygon is planar
 */
OperationResult op_is_planar(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

#endif /* OPERATIONS_VALIDITY_H */
