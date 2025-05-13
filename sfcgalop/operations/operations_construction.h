/**
 * operations_construction.h - Geometry construction operations
 */

#ifndef OPERATIONS_CONSTRUCTION_H
#define OPERATIONS_CONSTRUCTION_H

#include "operations.h"
#include "operations_params.h"

/**
 * Get intersection of geometries
 */
OperationResult op_intersection(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Get 3D intersection of geometries
 */
OperationResult op_intersection_3d(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Get difference of geometries
 */
OperationResult op_difference(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Get 3D difference of geometries
 */
OperationResult op_difference_3d(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Get union of geometries
 */
OperationResult op_union(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Get 3D union of geometries
 */
OperationResult op_union_3d(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Get convex hull of geometry
 */
OperationResult op_convexhull(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Get 3D convex hull of geometry
 */
OperationResult op_convexhull_3d(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Convert polyhedral surface to solid
 */
OperationResult op_make_solid(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Get boundary of geometry
 */
OperationResult op_boundary(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Get centroid of geometry
 */
OperationResult op_centroid(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Get 3D centroid of geometry
 */
OperationResult op_centroid_3d(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Get 2D envelope (bounding box) of geometry
 */
OperationResult op_envelope(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Get 3D envelope (bounding box) of geometry
 */
OperationResult op_envelope_3d(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

/**
 * Get Minkowski sum of geometries
 */
OperationResult op_minkowski_sum(const char* op_arg, sfcgal_geometry_t* geom_a, sfcgal_geometry_t* geom_b);

#endif /* OPERATIONS_CONSTRUCTION_H */
