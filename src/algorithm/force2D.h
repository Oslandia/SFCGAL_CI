// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _SFCGAL_ALGORITHM_FORCE2D_H_
#define _SFCGAL_ALGORITHM_FORCE2D_H_

#include <SFCGAL/config.h>

#include <SFCGAL/Kernel.h>

namespace SFCGAL {
class Geometry;
}

namespace SFCGAL {
namespace algorithm {

/**
 * @brief force a geometry to be 2D (project on O,x,y)
 * @warning ignore empty geometries
 */
SFCGAL_API void
force2D(Geometry &g);

} // namespace algorithm
} // namespace SFCGAL

#endif
