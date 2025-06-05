// Copyright (c) 2022-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_ALPHASHAPES_H_
#define SFCGAL_ALGORITHM_ALPHASHAPES_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Geometry.h"

namespace SFCGAL {
namespace algorithm {

/**
 * Compute the 2D alpha shapes for a geometry
 * https://doc.cgal.org/latest/Alpha_shapes_2/index.html#Chapter_2D_Alpha_Shapes
 * @since 1.4.1
 */
SFCGAL_API auto
alphaShapes(const Geometry &g, double alpha = 1, bool allow_holes = false)
    -> std::unique_ptr<Geometry>;

/**
 * Compute the optimal 2D alpha shapes for a geometry
 * https://doc.cgal.org/latest/Alpha_shapes_2/index.html#Chapter_2D_Alpha_Shapes
 * @since 1.4.1
 */
SFCGAL_API auto
optimal_alpha_shapes(const Geometry &g, bool allow_holes = false,
                     size_t nb_components = 1) -> std::unique_ptr<Geometry>;

} // namespace algorithm
} // namespace SFCGAL

#endif
