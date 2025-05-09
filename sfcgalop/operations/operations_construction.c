/**
 * operations_construction.c - Implementation of geometry construction
 * operations
 */

#include "operations_construction.h"
#include "../util.h"
#include "operations_common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Geometry operations - single geometry input
 */
DEFINE_GEOMETRY_OP_SINGLE(convexhull, sfcgal_geometry_convexhull)
DEFINE_GEOMETRY_OP_SINGLE(convexhull_3d, sfcgal_geometry_convexhull_3d)
DEFINE_GEOMETRY_OP_SINGLE(make_solid, sfcgal_geometry_make_solid)
DEFINE_GEOMETRY_OP_SINGLE(boundary, sfcgal_geometry_boundary)
DEFINE_GEOMETRY_OP_SINGLE(centroid, sfcgal_geometry_centroid)
DEFINE_GEOMETRY_OP_SINGLE(centroid_3d, sfcgal_geometry_centroid_3d)
DEFINE_GEOMETRY_OP_SINGLE(envelope, sfcgal_geometry_envelope)
DEFINE_GEOMETRY_OP_SINGLE(envelope_3d, sfcgal_geometry_envelope_3d)

/**
 * Geometry operations - two geometry inputs
 */
DEFINE_GEOMETRY_OP_BINARY(intersection, sfcgal_geometry_intersection)
DEFINE_GEOMETRY_OP_BINARY(intersection_3d, sfcgal_geometry_intersection_3d)
DEFINE_GEOMETRY_OP_BINARY(difference, sfcgal_geometry_difference)
DEFINE_GEOMETRY_OP_BINARY(difference_3d, sfcgal_geometry_difference_3d)
DEFINE_GEOMETRY_OP_BINARY(union, sfcgal_geometry_union)
DEFINE_GEOMETRY_OP_BINARY(union_3d, sfcgal_geometry_union_3d)
DEFINE_GEOMETRY_OP_BINARY(minkowski_sum, sfcgal_geometry_minkowski_sum)
