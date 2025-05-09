/**
 * operations_registry.c - Implementation of the operations registry (SECURITY
 * HARDENED)
 */

#include "operations_registry.h"
#include "operations_alphashapes.h"
#include "operations_construction.h"
#include "operations_coordinate.h"
#include "operations_metrics.h"
#include "operations_skeleton.h"
#include "operations_spatial_analysis.h"
#include "operations_tesselation.h"
#include "operations_transformation.h"
#include "operations_validity.h"
#include "operations_visibility.h"
#include <limits.h>
#include <string.h>

/* ========================================================================== */
/* PARAMETER DECLARATIONS */
/* ========================================================================== */

// Parameter declarations for force_z
const OperationParam force_z_params[] = {{"defaultZ", "Default Z value", 0.0}};

// Parameter declarations for force_m
const OperationParam force_m_params[] = {{"defaultM", "Default M value", 0.0}};

// Parameter declarations for almost_equals
const OperationParam almost_equals_params[] = {
    {"tolerance", "Tolerance for equality comparison", 0.0}};

// Parameter declarations for extrude
const OperationParam extrude_params[] = {
    {"dx", "X component of extrusion vector", 0.0},
    {"dy", "Y component of extrusion vector", 0.0},
    {"dz", "Z component of extrusion vector", 1.0}};

// Parameter declarations for buffer3d
const OperationParam buffer3d_params[] = {
    {"radius", "Buffer radius", 1.0},
    {"segments", "Number of segments for curves", 16},
    {"buffer_type", "Buffer type (0=round, 1=cylsphere, 2=flat)", 0.0}};

// Parameter declarations for round
const OperationParam round_params[] = {
    {"scale", "Scale factor for rounding", 0}};

// Parameter declarations for offset_polygon
const OperationParam offset_polygon_params[] = {
    {"radius", "Offset radius", 1.0}};

// Parameter declarations for line_substring
const OperationParam line_substring_params[] = {
    {"start", "Start parameter (0.0 to 1.0)", 0.0},
    {"end", "End parameter (0.0 to 1.0)", 1.0}};

// Parameter declarations for simplify
const OperationParam simplify_params[] = {
    {"threshold", "Simplification threshold", 1.0},
    {"preserve_topology", "Preserve topology during simplification", 1.0}};

// Parameter declarations for extrude_straight_skeleton
const OperationParam extrude_straight_skeleton_params[] = {
    {"height", "Extrusion height", 1.0}};

// Parameter declarations for extrude_polygon_straight_skeleton
const OperationParam extrude_polygon_straight_skeleton_params[] = {
    {"building_height", "Building height", 1.0},
    {"roof_height", "Roof height", 1.0}};

// Parameter declarations for straight_skeleton_partition
const OperationParam straight_skeleton_partition_params[] = {
    {"auto_orientation", "Auto-orientation flag (0=false, 1=true)", 1.0}};

// Parameter declarations for alpha_shapes
const OperationParam alpha_shapes_params[] = {
    {"alpha", "Alpha parameter", 0.0},
    {"allow_holes", "Allow holes flag (0=false, 1=true)", 0.0}};

// Parameter declarations for optimal_alpha_shapes
const OperationParam optimal_alpha_shapes_params[] = {
    {"allow_holes", "Allow holes flag (0=false, 1=true)", 0.0},
    {"nb_components", "Number of components", 1.0}};

// Parameter declarations for alpha_wrapping_3d
const OperationParam alpha_wrapping_3d_params[] = {
    {"relative_alpha", "Relative alpha parameter", 100.0},
    {"relative_offset", "Relative offset parameter", 100.0}};

// Parameter declarations for visibility_segment
const OperationParam visibility_segment_params[] = {
    {"point_b", "Second point for segment (WKT)", 0.0}};

// Parameter declarations for rotate
const OperationParam rotate_params[] = {
    {"angle", "Rotation angle in radians", 0.0}};

// Parameter declarations for rotate_2d
const OperationParam rotate_2d_params[] = {
    {"angle", "Rotation angle in radians", 0.0},
    {"cx", "X coordinate of rotation center", 0.0},
    {"cy", "Y coordinate of rotation center", 0.0}};

// Parameter declarations for rotate_3d
const OperationParam rotate_3d_params[] = {
    {"angle", "Rotation angle in radians", 0.0},
    {"ax", "X component of rotation axis", 0.0},
    {"ay", "Y component of rotation axis", 0.0},
    {"az", "Z component of rotation axis", 1.0}};

// Parameter declarations for rotate_3d_around_center
const OperationParam rotate_3d_around_center_params[] = {
    {"angle", "Rotation angle in radians", 0.0},
    {"ax", "X component of rotation axis", 0.0},
    {"ay", "Y component of rotation axis", 0.0},
    {"az", "Z component of rotation axis", 1.0},
    {"cx", "X coordinate of rotation center", 0.0},
    {"cy", "Y coordinate of rotation center", 0.0},
    {"cz", "Z coordinate of rotation center", 0.0}};

// Parameter declarations for rotate_x
const OperationParam rotate_x_params[] = {
    {"angle", "Rotation angle in radians", 0.0}};

// Parameter declarations for rotate_y
const OperationParam rotate_y_params[] = {
    {"angle", "Rotation angle in radians", 0.0}};

// Parameter declarations for rotate_z
const OperationParam rotate_z_params[] = {
    {"angle", "Rotation angle in radians", 0.0}};

// Parameter declarations for scale
const OperationParam scale_params[] = {{"scale", "Scale factor", 1.0}};

// Parameter declarations for scale_3d
const OperationParam scale_3d_params[] = {{"sx", "X scale factor", 1.0},
                                          {"sy", "Y scale factor", 1.0},
                                          {"sz", "Z scale factor", 1.0}};

// Parameter declarations for scale_3d_around_center
const OperationParam scale_3d_around_center_params[] = {
    {"sx", "X scale factor", 1.0},
    {"sy", "Y scale factor", 1.0},
    {"sz", "Z scale factor", 1.0},
    {"cx", "X coordinate of scaling center", 0.0},
    {"cy", "Y coordinate of scaling center", 0.0},
    {"cz", "Z coordinate of scaling center", 0.0}};

// Parameter declarations for translate_2d
const OperationParam translate_2d_params[] = {{"dx", "X translation", 0.0},
                                              {"dy", "Y translation", 0.0}};

// Parameter declarations for translate_3d
const OperationParam translate_3d_params[] = {{"dx", "X translation", 0.0},
                                              {"dy", "Y translation", 0.0},
                                              {"dz", "Z translation", 0.0}};

/* ========================================================================== */
/* OPERATIONS REGISTRY */
/* ========================================================================== */

// Array of supported operations
const Operation operations[] = {
    /* Validity operations */
    {"is_valid", op_is_valid, false, "Test if geometry is valid", "Validity",
     RESULT_BOOLEAN, false, NULL, 0},
    {"is_validity_detail", op_is_validity_detail, false,
     "Get detailed information about geometry validity", "Validity",
     RESULT_TEXT, false, NULL, 0},
    {"is_planar", op_is_planar, false, "Test if polygon is planar", "Validity",
     RESULT_BOOLEAN, false, NULL, 0},

    /* Coordinate handling operations */
    {"drop_z", op_drop_z, false, "Remove Z coordinate from geometry",
     "Coordinate Handling", RESULT_GEOMETRY, false, NULL, 0},
    {"drop_m", op_drop_m, false, "Remove M coordinate from geometry",
     "Coordinate Handling", RESULT_GEOMETRY, false, NULL, 0},
    {"force_z", op_force_z, false, "Add Z coordinate with given default value",
     "Coordinate Handling", RESULT_GEOMETRY, true, force_z_params, 1},
    {"force_m", op_force_m, false, "Add M coordinate with given default value",
     "Coordinate Handling", RESULT_GEOMETRY, true, force_m_params, 1},

    /* Information operations */
    {"num_geometries", op_num_geometries, false,
     "Get number of geometries in a collection", "Information", RESULT_NUMERIC,
     false, NULL, 0},
    {"orientation", op_orientation, false,
     "Get orientation of polygon (-1=CCW, 1=CW, 0=invalid)", "Information",
     RESULT_NUMERIC, false, NULL, 0},

    /* Spatial analysis operations */
    {"intersects", op_intersects, true, "Check if geometries intersect",
     "Spatial Analysis", RESULT_BOOLEAN, false, NULL, 0},
    {"intersects_3d", op_intersects_3d, true,
     "Check if geometries intersect in 3D", "Spatial Analysis", RESULT_BOOLEAN,
     false, NULL, 0},
    {"covers", op_covers, true, "Check if geometry A covers geometry B",
     "Spatial Analysis", RESULT_BOOLEAN, false, NULL, 0},
    {"covers_3d", op_covers_3d, true,
     "Check if geometry A covers geometry B in 3D", "Spatial Analysis",
     RESULT_BOOLEAN, false, NULL, 0},
    {"equals", op_equals, true, "Check if geometries are equal",
     "Spatial Analysis", RESULT_BOOLEAN, false, NULL, 0},
    {"almost_equals", op_almost_equals, true,
     "Check if geometries are almost equal (with tolerance)",
     "Spatial Analysis", RESULT_BOOLEAN, true, almost_equals_params, 1},

    /* Construction operations */
    {"intersection", op_intersection, true, "Get intersection of geometries",
     "Construction", RESULT_GEOMETRY, false, NULL, 0},
    {"intersection_3d", op_intersection_3d, true,
     "Get 3D intersection of geometries", "Construction", RESULT_GEOMETRY,
     false, NULL, 0},
    {"difference", op_difference, true, "Get difference of geometries",
     "Construction", RESULT_GEOMETRY, false, NULL, 0},
    {"difference_3d", op_difference_3d, true, "Get 3D difference of geometries",
     "Construction", RESULT_GEOMETRY, false, NULL, 0},
    {"union", op_union, true, "Get union of geometries", "Construction",
     RESULT_GEOMETRY, false, NULL, 0},
    {"union_3d", op_union_3d, true, "Get 3D union of geometries",
     "Construction", RESULT_GEOMETRY, false, NULL, 0},
    {"convexhull", op_convexhull, false, "Get convex hull of geometry",
     "Construction", RESULT_GEOMETRY, false, NULL, 0},
    {"convexhull_3d", op_convexhull_3d, false, "Get 3D convex hull of geometry",
     "Construction", RESULT_GEOMETRY, false, NULL, 0},
    {"make_solid", op_make_solid, false, "Convert polyhedral surface to solid",
     "Construction", RESULT_GEOMETRY, false, NULL, 0},
    {"boundary", op_boundary, false, "Get boundary of geometry", "Construction",
     RESULT_GEOMETRY, false, NULL, 0},
    {"centroid", op_centroid, false, "Get centroid of geometry", "Construction",
     RESULT_GEOMETRY, false, NULL, 0},
    {"centroid_3d", op_centroid_3d, false, "Get 3D centroid of geometry",
     "Construction", RESULT_GEOMETRY, false, NULL, 0},
    {"envelope", op_envelope, false,
     "Get 2D envelope (bounding box) of geometry", "Construction",
     RESULT_GEOMETRY, false, NULL, 0},
    {"envelope_3d", op_envelope_3d, false,
     "Get 3D envelope (bounding box) of geometry", "Construction",
     RESULT_GEOMETRY, false, NULL, 0},
    {"minkowski_sum", op_minkowski_sum, true, "Get Minkowski sum of geometries",
     "Construction", RESULT_GEOMETRY, false, NULL, 0},

    /* Tesselation operations */
    {"tesselate", op_tesselate, false, "Tesselate geometry into triangles",
     "Tesselation", RESULT_GEOMETRY, false, NULL, 0},
    {"triangulate_2dz", op_triangulate_2dz, false,
     "Triangulate 2D geometry preserving Z", "Tesselation", RESULT_GEOMETRY,
     false, NULL, 0},
    {"y_monotone_partition", op_y_monotone_partition, false,
     "Partition into y-monotone polygons", "Tesselation", RESULT_GEOMETRY,
     false, NULL, 0},
    {"approx_convex_partition", op_approx_convex_partition, false,
     "Partition into approx. convex polygons", "Tesselation", RESULT_GEOMETRY,
     false, NULL, 0},
    {"greene_approx_convex_partition", op_greene_approx_convex_partition, false,
     "Greene's approx. convex partition", "Tesselation", RESULT_GEOMETRY, false,
     NULL, 0},
    {"optimal_convex_partition", op_optimal_convex_partition, false,
     "Optimal convex partition", "Tesselation", RESULT_GEOMETRY, false, NULL,
     0},

    /* Orientation operations */
    {"force_lhr", op_force_lhr, false, "Force left-hand rule orientation",
     "Orientation", RESULT_GEOMETRY, false, NULL, 0},
    {"force_rhr", op_force_rhr, false, "Force right-hand rule orientation",
     "Orientation", RESULT_GEOMETRY, false, NULL, 0},

    /* Metrics operations */
    {"distance", op_distance, true, "Calculate distance between geometries",
     "Metrics", RESULT_NUMERIC, false, NULL, 0},
    {"distance_3d", op_distance_3d, true,
     "Calculate 3D distance between geometries", "Metrics", RESULT_NUMERIC,
     false, NULL, 0},
    {"volume", op_volume, false, "Calculate volume of 3D geometry", "Metrics",
     RESULT_NUMERIC, false, NULL, 0},
    {"area", op_area, false, "Calculate area of 2D geometry", "Metrics",
     RESULT_NUMERIC, false, NULL, 0},
    {"area_3d", op_area_3d, false, "Calculate area of 3D geometry", "Metrics",
     RESULT_NUMERIC, false, NULL, 0},
    {"length", op_length, false, "Calculate length of geometry", "Metrics",
     RESULT_NUMERIC, false, NULL, 0},
    {"length_3d", op_length_3d, false, "Calculate 3D length of geometry",
     "Metrics", RESULT_NUMERIC, false, NULL, 0},

    /* Advanced operations */
    {"extrude", op_extrude, false, "Extrude 2D geometry to 3D", "Advanced",
     RESULT_GEOMETRY, true, extrude_params, 3},
    {"buffer3d", op_buffer3d, false, "Create 3D buffer around geometry",
     "Advanced", RESULT_GEOMETRY, true, buffer3d_params, 3},
    {"round", op_round, false, "Round coordinates to given precision",
     "Advanced", RESULT_GEOMETRY, true, round_params, 1},
    {"offset_polygon", op_offset_polygon, false, "Create offset polygon",
     "Advanced", RESULT_GEOMETRY, true, offset_polygon_params, 1},
    {"line_substring", op_line_substring, false, "Get substring of linestring",
     "Advanced", RESULT_GEOMETRY, true, line_substring_params, 2},
    {"simplify", op_simplify, false, "Simplify geometry", "Advanced",
     RESULT_GEOMETRY, true, simplify_params, 2},

    /* Skeleton operations */
    {"straight_skeleton", op_straight_skeleton, false,
     "Compute straight skeleton", "Skeleton", RESULT_GEOMETRY, false, NULL, 0},
    {"straight_skeleton_distance_in_m", op_straight_skeleton_distance_in_m,
     false, "Compute straight skeleton with distances", "Skeleton",
     RESULT_GEOMETRY, false, NULL, 0},
    {"extrude_straight_skeleton", op_extrude_straight_skeleton, false,
     "Extrude using straight skeleton", "Skeleton", RESULT_GEOMETRY, true,
     extrude_straight_skeleton_params, 1},
    {"extrude_polygon_straight_skeleton", op_extrude_polygon_straight_skeleton,
     false, "Extrude polygon with roof", "Skeleton", RESULT_GEOMETRY, true,
     extrude_polygon_straight_skeleton_params, 2},
    {"approximate_medial_axis", op_approximate_medial_axis, false,
     "Compute approximate medial axis", "Skeleton", RESULT_GEOMETRY, false,
     NULL, 0},
    {"straight_skeleton_partition", op_straight_skeleton_partition, false,
     "Partition using straight skeleton", "Skeleton", RESULT_GEOMETRY, true,
     straight_skeleton_partition_params, 1},

/* Alpha shape operations - Fix for _MSC_VER undefined warning */
#if !defined(_MSC_VER)
    {"alpha_shapes", op_alpha_shapes, false, "Compute alpha shapes",
     "Alpha Shapes", RESULT_GEOMETRY, true, alpha_shapes_params, 2},
    {"optimal_alpha_shapes", op_optimal_alpha_shapes, false,
     "Compute optimal alpha shapes", "Alpha Shapes", RESULT_GEOMETRY, true,
     optimal_alpha_shapes_params, 2},
#endif /* !defined(_MSC_VER) */
    {"alpha_wrapping_3d", op_alpha_wrapping_3d, false,
     "Compute 3D alpha wrapping", "Alpha Shapes", RESULT_GEOMETRY, true,
     alpha_wrapping_3d_params, 2},

    /* Visibility operations */
    {"visibility_point", op_visibility_point, true,
     "Compute visibility polygon from point", "Visibility", RESULT_GEOMETRY,
     false, NULL, 0},
    {"visibility_segment", op_visibility_segment, true,
     "Compute visibility polygon from segment", "Visibility", RESULT_GEOMETRY,
     true, visibility_segment_params, 1},

    /* Transformation operations */
    {"rotate", op_rotate, false, "Rotate geometry", "Transformation",
     RESULT_GEOMETRY, true, rotate_params, 1},
    {"rotate_2d", op_rotate_2d, false, "Rotate geometry in 2D around point",
     "Transformation", RESULT_GEOMETRY, true, rotate_2d_params, 3},
    {"rotate_3d", op_rotate_3d, false, "Rotate geometry in 3D around axis",
     "Transformation", RESULT_GEOMETRY, true, rotate_3d_params, 4},
    {"rotate_3d_around_center", op_rotate_3d_around_center, false,
     "Rotate geometry in 3D around axis and center", "Transformation",
     RESULT_GEOMETRY, true, rotate_3d_around_center_params, 7},
    {"rotate_x", op_rotate_x, false, "Rotate geometry around X axis",
     "Transformation", RESULT_GEOMETRY, true, rotate_x_params, 1},
    {"rotate_y", op_rotate_y, false, "Rotate geometry around Y axis",
     "Transformation", RESULT_GEOMETRY, true, rotate_y_params, 1},
    {"rotate_z", op_rotate_z, false, "Rotate geometry around Z axis",
     "Transformation", RESULT_GEOMETRY, true, rotate_z_params, 1},
    {"scale", op_scale, false, "Scale geometry uniformly", "Transformation",
     RESULT_GEOMETRY, true, scale_params, 1},
    {"scale_3d", op_scale_3d, false, "Scale geometry in 3D", "Transformation",
     RESULT_GEOMETRY, true, scale_3d_params, 3},
    {"scale_3d_around_center", op_scale_3d_around_center, false,
     "Scale geometry in 3D around center", "Transformation", RESULT_GEOMETRY,
     true, scale_3d_around_center_params, 6},
    {"translate_2d", op_translate_2d, false, "Translate geometry in 2D",
     "Transformation", RESULT_GEOMETRY, true, translate_2d_params, 2},
    {"translate_3d", op_translate_3d, false, "Translate geometry in 3D",
     "Transformation", RESULT_GEOMETRY, true, translate_3d_params, 3},

    {NULL, NULL, false, NULL, NULL, RESULT_BOOLEAN, false, NULL,
     0} /* sentinel value */
};
