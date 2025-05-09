/**
 * operations_transformation.h - Geometry transformation operations
 */

#ifndef OPERATIONS_TRANSFORMATION_H
#define OPERATIONS_TRANSFORMATION_H

#include "operations.h"
#include "operations_params.h"

/**
 * Rotation operations
 */
OperationResult
op_rotate(const char *op_arg, const sfcgal_geometry_t *geom_a,
          const sfcgal_geometry_t *geom_b);
OperationResult
op_rotate_2d(const char *op_arg, const sfcgal_geometry_t *geom_a,
             const sfcgal_geometry_t *geom_b);
OperationResult
op_rotate_3d(const char *op_arg, const sfcgal_geometry_t *geom_a,
             const sfcgal_geometry_t *geom_b);
OperationResult
op_rotate_3d_around_center(const char *op_arg, const sfcgal_geometry_t *geom_a,
                           const sfcgal_geometry_t *geom_b);
OperationResult
op_rotate_x(const char *op_arg, const sfcgal_geometry_t *geom_a,
            const sfcgal_geometry_t *geom_b);
OperationResult
op_rotate_y(const char *op_arg, const sfcgal_geometry_t *geom_a,
            const sfcgal_geometry_t *geom_b);
OperationResult
op_rotate_z(const char *op_arg, const sfcgal_geometry_t *geom_a,
            const sfcgal_geometry_t *geom_b);

/**
 * Scale operations
 */
OperationResult
op_scale(const char *op_arg, const sfcgal_geometry_t *geom_a,
         const sfcgal_geometry_t *geom_b);
OperationResult
op_scale_3d(const char *op_arg, const sfcgal_geometry_t *geom_a,
            const sfcgal_geometry_t *geom_b);
OperationResult
op_scale_3d_around_center(const char *op_arg, const sfcgal_geometry_t *geom_a,
                          const sfcgal_geometry_t *geom_b);

/**
 * Translation operations
 */
OperationResult
op_translate_2d(const char *op_arg, const sfcgal_geometry_t *geom_a,
                const sfcgal_geometry_t *geom_b);
OperationResult
op_translate_3d(const char *op_arg, const sfcgal_geometry_t *geom_a,
                const sfcgal_geometry_t *geom_b);

/**
 * Other transformation operations
 */
OperationResult
op_extrude(const char *op_arg, const sfcgal_geometry_t *geom_a,
           const sfcgal_geometry_t *geom_b);
OperationResult
op_buffer3d(const char *op_arg, const sfcgal_geometry_t *geom_a,
            const sfcgal_geometry_t *geom_b);
OperationResult
op_round(const char *op_arg, const sfcgal_geometry_t *geom_a,
         const sfcgal_geometry_t *geom_b);
OperationResult
op_offset_polygon(const char *op_arg, const sfcgal_geometry_t *geom_a,
                  const sfcgal_geometry_t *geom_b);
OperationResult
op_line_substring(const char *op_arg, const sfcgal_geometry_t *geom_a,
                  const sfcgal_geometry_t *geom_b);
OperationResult
op_simplify(const char *op_arg, const sfcgal_geometry_t *geom_a,
            const sfcgal_geometry_t *geom_b);

#endif /* OPERATIONS_TRANSFORMATION_H */
