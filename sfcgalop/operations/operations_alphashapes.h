/**
 * operations_alphashapes.h - Alpha shapes operations for SFCGAL (SECURITY
 * HARDENED)
 */

#ifndef OPERATIONS_ALPHASHAPES_H
#define OPERATIONS_ALPHASHAPES_H

#include "operations.h"
#include "operations_params.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Fix for _MSC_VER undefined warning */
#if !defined(_MSC_VER)
/**
 * Compute alpha shapes
 *
 * @param op_arg Operation arguments (alpha=value,allow_holes=bool)
 * @param geom_a Input geometry
 * @param geom_b Unused (should be NULL)
 * @return Operation result containing computed alpha shapes
 */
OperationResult
op_alpha_shapes(const char *op_arg, const sfcgal_geometry_t *geom_a,
                const sfcgal_geometry_t *geom_b);

/**
 * Compute optimal alpha shapes
 *
 * @param op_arg Operation arguments (allow_holes=bool,nb_components=int)
 * @param geom_a Input geometry
 * @param geom_b Unused (should be NULL)
 * @return Operation result containing computed optimal alpha shapes
 */
OperationResult
op_optimal_alpha_shapes(const char *op_arg, const sfcgal_geometry_t *geom_a,
                        const sfcgal_geometry_t *geom_b);
#endif /* !defined(_MSC_VER) */

/**
 * Compute 3D alpha wrapping
 *
 * @param op_arg Operation arguments (relative_alpha=int,relative_offset=int)
 * @param geom_a Input geometry
 * @param geom_b Unused (should be NULL)
 * @return Operation result containing computed 3D alpha wrapping
 */
OperationResult
op_alpha_wrapping_3d(const char *op_arg, const sfcgal_geometry_t *geom_a,
                     const sfcgal_geometry_t *geom_b);

#ifdef __cplusplus
}
#endif

#endif /* OPERATIONS_ALPHASHAPES_H */
