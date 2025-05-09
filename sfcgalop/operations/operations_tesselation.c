/**
 * operations_tesselation.c - Implementation of geometry tesselation operations
 */

#include "operations_tesselation.h"
#include "../util.h"
#include "operations_common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Tesselation operations - single geometry input
 */
DEFINE_GEOMETRY_OP_SINGLE(tesselate, sfcgal_geometry_tesselate)
DEFINE_GEOMETRY_OP_SINGLE(triangulate_2dz, sfcgal_geometry_triangulate_2dz)
DEFINE_GEOMETRY_OP_SINGLE(y_monotone_partition, sfcgal_y_monotone_partition_2)
DEFINE_GEOMETRY_OP_SINGLE(approx_convex_partition,
                          sfcgal_approx_convex_partition_2)
DEFINE_GEOMETRY_OP_SINGLE(greene_approx_convex_partition,
                          sfcgal_greene_approx_convex_partition_2)
DEFINE_GEOMETRY_OP_SINGLE(optimal_convex_partition,
                          sfcgal_optimal_convex_partition_2)
