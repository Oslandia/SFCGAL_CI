// Copyright (c) 2025-2025, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_SIMPLIFICATION_H
#define SFCGAL_ALGORITHM_SIMPLIFICATION_H

#include <SFCGAL/Exception.h>
#include <SFCGAL/Geometry.h>
#include <SFCGAL/config.h>

namespace SFCGAL {
namespace algorithm {

/**
 * @brief Simplifies a geometry using the CGAL algorithm
 *  https://doc.cgal.org/latest/Polyline_simplification_2/index.html
 *
 * @param geometry The geometry to simplify
 * @param threshold The distance threshold for simplification
 * @param preserveTopology Whether to preserve topology during simplification
 * @return A simplified copy of the input geometry
 */
SFCGAL_API auto
simplify(const Geometry &geometry, double threshold, bool preserveTopology)
    -> std::unique_ptr<Geometry>;

} // namespace algorithm
} // namespace SFCGAL

#endif
