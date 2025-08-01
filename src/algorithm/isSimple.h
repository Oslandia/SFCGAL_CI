// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_ISSIMPLE_H_
#define SFCGAL_ALGORITHM_ISSIMPLE_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Simplicity.h"

namespace SFCGAL {

/**
 * Functions used to assert for geometry simplicity
 */
void SFCGAL_API
SFCGAL_ASSERT_GEOMETRY_SIMPLICITY(const Geometry &g);

namespace algorithm {

/**
 * Check simplicity of a geometry
 */
SFCGAL_API const Simplicity
isSimple(const Geometry &g, const double &toleranceAbs = 1e-9);

} // namespace algorithm
} // namespace SFCGAL

#endif
