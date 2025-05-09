/**
 * operations_spatial_analysis.h - Spatial analysis operations
 */

#ifndef OPERATIONS_SPATIAL_ANALYSIS_H
#define OPERATIONS_SPATIAL_ANALYSIS_H

#include "operations.h"
#include "operations_params.h"

/**
 * Check if geometries intersect
 */
OperationResult
op_intersects(const char *op_arg, const sfcgal_geometry_t *geom_a,
              const sfcgal_geometry_t *geom_b);

/**
 * Check if geometries intersect in 3D
 */
OperationResult
op_intersects_3d(const char *op_arg, const sfcgal_geometry_t *geom_a,
                 const sfcgal_geometry_t *geom_b);

/**
 * Check if geometry A covers geometry B
 */
OperationResult
op_covers(const char *op_arg, const sfcgal_geometry_t *geom_a,
          const sfcgal_geometry_t *geom_b);

/**
 * Check if geometry A covers geometry B in 3D
 */
OperationResult
op_covers_3d(const char *op_arg, const sfcgal_geometry_t *geom_a,
             const sfcgal_geometry_t *geom_b);

/**
 * Check if geometries are equal
 */
OperationResult
op_equals(const char *op_arg, const sfcgal_geometry_t *geom_a,
          const sfcgal_geometry_t *geom_b);

/**
 * Check if geometries are almost equal with tolerance
 */
OperationResult
op_almost_equals(const char *op_arg, const sfcgal_geometry_t *geom_a,
                 const sfcgal_geometry_t *geom_b);

#endif /* OPERATIONS_SPATIAL_ANALYSIS_H */
