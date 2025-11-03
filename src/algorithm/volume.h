// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_VOLUME_H_
#define SFCGAL_ALGORITHM_VOLUME_H_

#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/config.h"

namespace SFCGAL::algorithm {

/**
 * @brief Calculate the volume of a closed PolyhedralSurface
 * @param surface The surface whose volume to calculate
 * @param noValidityCheck Tag to disable validity checking
 * @return The volume of the surface (exact arithmetic)
 * @throws GeometryInvalidityException if the surface is not closed
 * @warning The surface must be closed and properly oriented
 */
SFCGAL_API auto
volume(const PolyhedralSurface &surface,
       NoValidityCheck noValidityCheck = NoValidityCheck()) -> Kernel::FT;

/**
 * @brief Calculate the volume of a Solid
 * @param solid The solid whose volume to calculate
 * @param noValidityCheck Tag to disable validity checking
 * @return The volume of the solid (exact arithmetic)
 * @note Outer shell contributes positively, inner shells (holes) contribute
 * negatively
 */
SFCGAL_API auto
volume(const Solid &solid, NoValidityCheck noValidityCheck = NoValidityCheck())
    -> Kernel::FT;

/**
 * @brief Calculate the volume of any geometry
 * @param g The geometry whose volume to calculate
 * @return The volume of the geometry (0 for 2D geometries)
 * @note Returns 0 for geometries of dimension < 3
 * @note For MultiSolid and GeometryCollection, returns the sum of volumes
 */
SFCGAL_API auto
volume(const Geometry &g) -> Kernel::FT;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_VOLUME_H_
