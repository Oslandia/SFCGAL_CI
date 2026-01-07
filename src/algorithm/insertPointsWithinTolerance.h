// Copyright (c) 2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_INSERTPOINTSWITHINTOLERANCE_H_
#define SFCGAL_ALGORITHM_INSERTPOINTSWITHINTOLERANCE_H_

#include "SFCGAL/config.h"
#include "SFCGAL/numeric.h"

#include <memory>

namespace SFCGAL {
class Geometry;

namespace algorithm {

/**
 * @brief Inserts points from sourceGeometry into baseGeometry within a given
 * tolerance
 *
 * This function takes a base geometry and adds points from the source geometry
 * to the base geometry where they are within the specified tolerance distance
 * from the base geometry's segments. This is useful for topological editing
 * operations where you want to densify a geometry with points from another.
 *
 * @param baseGeometry The geometry to be modified by adding points
 * @param sourceGeometry The geometry providing points to insert
 * @param tolerance The maximum distance for a point to be considered for
 * insertion
 * @return A new geometry with the points from sourceGeometry inserted into
 * baseGeometry
 */
SFCGAL_API auto
insertPointsWithinTolerance(const Geometry &baseGeometry,
                            const Geometry &sourceGeometry,
                            double          tolerance = EPSILON)
    -> std::unique_ptr<Geometry>;

} // namespace algorithm
} // namespace SFCGAL

#endif
