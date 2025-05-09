/**
 * operations_registry.h - Registry of operations and parameters
 */

#ifndef OPERATIONS_REGISTRY_H
#define OPERATIONS_REGISTRY_H

#include "operations.h"

/**
 * Parameter declarations for various operations
 */
extern const OperationParam force_z_params[];
extern const OperationParam force_m_params[];
extern const OperationParam almost_equals_params[];
extern const OperationParam extrude_params[];
extern const OperationParam buffer3d_params[];
extern const OperationParam round_params[];
extern const OperationParam offset_polygon_params[];
extern const OperationParam line_substring_params[];
extern const OperationParam simplify_params[];
extern const OperationParam extrude_straight_skeleton_params[];
extern const OperationParam extrude_polygon_straight_skeleton_params[];
extern const OperationParam straight_skeleton_partition_params[];
extern const OperationParam alpha_shapes_params[];
extern const OperationParam optimal_alpha_shapes_params[];
extern const OperationParam alpha_wrapping_3d_params[];
extern const OperationParam visibility_segment_params[];
extern const OperationParam rotate_params[];
extern const OperationParam rotate_2d_params[];
extern const OperationParam rotate_3d_params[];
extern const OperationParam rotate_3d_around_center_params[];
extern const OperationParam rotate_x_params[];
extern const OperationParam rotate_y_params[];
extern const OperationParam rotate_z_params[];
extern const OperationParam scale_params[];
extern const OperationParam scale_3d_params[];
extern const OperationParam scale_3d_around_center_params[];
extern const OperationParam translate_2d_params[];
extern const OperationParam translate_3d_params[];

/**
 * The global operations registry
 */
extern const Operation operations[];

/**
 * Initialize the operations registry
 */
void
init_operations_registry(void);

#endif /* OPERATIONS_REGISTRY_H */
