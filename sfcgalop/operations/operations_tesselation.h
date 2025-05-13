/**
 * operations_tesselation.h - Geometry tesselation operations
 */

#ifndef OPERATIONS_TESSELATION_H
#define OPERATIONS_TESSELATION_H

#include "operations.h"
#include "operations_params.h"

/**
 * Tesselate geometry into triangles
 */
OperationResult op_tesselate(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Triangulate 2D geometry preserving Z
 */
OperationResult op_triangulate_2dz(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Partition into y-monotone polygons
 */
OperationResult op_y_monotone_partition(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Partition into approx. convex polygons
 */
OperationResult op_approx_convex_partition(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Greene's approx. convex partition
 */
OperationResult op_greene_approx_convex_partition(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Optimal convex partition
 */
OperationResult op_optimal_convex_partition(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

#endif /* OPERATIONS_TESSELATION_H */
