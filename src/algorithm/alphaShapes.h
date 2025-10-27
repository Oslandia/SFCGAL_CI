// Copyright (c) 2022-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_ALPHASHAPES_H_
#define SFCGAL_ALGORITHM_ALPHASHAPES_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Geometry.h"

namespace SFCGAL::algorithm {

/**
 * Compute the 2D alpha shapes for a geometry
 * https://doc.cgal.org/latest/Alpha_shapes_2/index.html#Chapter_2D_Alpha_Shapes
 * @param geometry the input geometry
 * @param alpha the alpha value (default: 1)
 * @param allow_holes whether to allow holes in the result (default: false)
 * @return the computed alpha shapes as a unique_ptr<Geometry>
 * @since 1.4.1
 */
SFCGAL_API auto
alphaShapes(const Geometry &geometry, double alpha = 1,
            bool allow_holes = false) -> std::unique_ptr<Geometry>;

/**
 * Compute the optimal 2D alpha shapes for a geometry
 * https://doc.cgal.org/latest/Alpha_shapes_2/index.html#Chapter_2D_Alpha_Shapes
 * @param geometry the input geometry
 * @param allow_holes whether to allow holes in the result (default: false)
 * @param nb_components the number of components (default: 1)
 * @return the computed optimal alpha shapes as a unique_ptr<Geometry>
 * @since 1.4.1
 */
SFCGAL_API auto
optimal_alpha_shapes(const Geometry &geometry, bool allow_holes = false,
                     size_t nb_components = 1) -> std::unique_ptr<Geometry>;

} // namespace SFCGAL::algorithm

#endif
