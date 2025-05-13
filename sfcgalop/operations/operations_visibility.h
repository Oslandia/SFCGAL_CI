/**
 * operations_visibility.h - Visibility operations for SFCGAL
 */

#ifndef OPERATIONS_VISIBILITY_H
#define OPERATIONS_VISIBILITY_H

#include "operations.h"
#include "operations_params.h"

/**
 * Compute visibility polygon from point
 */
OperationResult op_visibility_point(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Compute visibility polygon from segment
 */
OperationResult op_visibility_segment(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

#endif /* OPERATIONS_VISIBILITY_H */
