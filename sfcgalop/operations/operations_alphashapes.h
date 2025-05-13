/**
 * operations_alphashapes.h - Alpha shapes operations for SFCGAL
 */

#ifndef OPERATIONS_ALPHASHAPES_H
#define OPERATIONS_ALPHASHAPES_H

#include "operations.h"
#include "operations_params.h"

#if !defined(_MSC_VER)
/**
 * Compute alpha shapes
 */
OperationResult op_alpha_shapes(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Compute optimal alpha shapes
 */
OperationResult op_optimal_alpha_shapes(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);
#endif

/**
 * Compute 3D alpha wrapping
 */
OperationResult op_alpha_wrapping_3d(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

#endif /* OPERATIONS_ALPHASHAPES_H */
