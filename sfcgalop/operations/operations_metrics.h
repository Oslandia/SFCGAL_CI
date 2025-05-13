/**
 * operations_metrics.h - Geometry metrics operations
 */

#ifndef OPERATIONS_METRICS_H
#define OPERATIONS_METRICS_H

#include "operations.h"
#include "operations_params.h"

/**
 * Calculate distance between geometries
 */
OperationResult op_distance(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Calculate 3D distance between geometries
 */
OperationResult op_distance_3d(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Calculate volume of 3D geometry
 */
OperationResult op_volume(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Calculate area of 2D geometry
 */
OperationResult op_area(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Calculate area of 3D geometry
 */
OperationResult op_area_3d(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Calculate length of geometry
 */
OperationResult op_length(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Calculate 3D length of geometry
 */
OperationResult op_length_3d(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Get orientation of polygon (-1=CCW, 1=CW, 0=invalid)
 */
OperationResult op_orientation(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Get number of geometries in a collection
 */
OperationResult op_num_geometries(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

#endif /* OPERATIONS_METRICS_H */
