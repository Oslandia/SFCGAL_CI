// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TRANSLATE_H_
#define SFCGAL_TRANSLATE_H_

#include "SFCGAL/export.h"

#include "SFCGAL/Kernel.h"

namespace SFCGAL {
class Geometry;
}

namespace SFCGAL::algorithm {

/**
 * @brief translate a geometry from a given vector
 * @param geometry The geometry to translate
 * @param vector The 3D translation vector
 */
SFCGAL_API void
translate(Geometry &geometry, const Kernel::Vector_3 &vector);
/**
 * @brief translate a geometry from a given vector
 * @param geometry The geometry to translate
 * @param vector The 2D translation vector
 */
SFCGAL_API void
translate(Geometry &geometry, const Kernel::Vector_2 &vector);
/**
 * @brief translate a geometry from a given vector
 * @param geometry The geometry to translate
 * @param dx Translation distance in X direction
 * @param dy Translation distance in Y direction
 * @param dz Translation distance in Z direction
 */
SFCGAL_API void
translate(Geometry &geometry, const Kernel::FT &dx, const Kernel::FT &dy,
          const Kernel::FT &dz);

/**
 * @brief translate a geometry from double-coordinates
 * @param geometry The geometry to translate
 * @param dx Translation distance in X direction
 * @param dy Translation distance in Y direction
 * @param dz Translation distance in Z direction
 */
SFCGAL_API void
translate(Geometry &geometry, const double &dx, const double &dy,
          const double &dz);

} // namespace SFCGAL::algorithm

#endif
