/**
 * operations_skeleton.h - Skeleton operations for SFCGAL
 */

#ifndef OPERATIONS_SKELETON_H
#define OPERATIONS_SKELETON_H

#include "operations.h"
#include "operations_params.h"

/**
 * Compute straight skeleton
 */
OperationResult
op_straight_skeleton(const char *op_arg, const sfcgal_geometry_t *geom_a,
                     const sfcgal_geometry_t *geom_b);

/**
 * Compute straight skeleton with distances
 */
OperationResult
op_straight_skeleton_distance_in_m(const char              *op_arg,
                                   const sfcgal_geometry_t *geom_a,
                                   const sfcgal_geometry_t *geom_b);

/**
 * Compute approximate medial axis
 */
OperationResult
op_approximate_medial_axis(const char *op_arg, const sfcgal_geometry_t *geom_a,
                           const sfcgal_geometry_t *geom_b);

/**
 * Extrude using straight skeleton
 */
OperationResult
op_extrude_straight_skeleton(const char              *op_arg,
                             const sfcgal_geometry_t *geom_a,
                             const sfcgal_geometry_t *geom_b);

/**
 * Extrude polygon with roof
 */
OperationResult
op_extrude_polygon_straight_skeleton(const char              *op_arg,
                                     const sfcgal_geometry_t *geom_a,
                                     const sfcgal_geometry_t *geom_b);

/**
 * Partition using straight skeleton
 */
OperationResult
op_straight_skeleton_partition(const char              *op_arg,
                               const sfcgal_geometry_t *geom_a,
                               const sfcgal_geometry_t *geom_b);

#endif /* OPERATIONS_SKELETON_H */
