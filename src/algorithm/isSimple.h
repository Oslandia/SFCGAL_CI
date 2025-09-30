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
 * @brief Assert that a geometry is simple
 * @param g The geometry to check for simplicity
 * @return void (throws exception if geometry is not simple)
 */
void SFCGAL_API
SFCGAL_ASSERT_GEOMETRY_SIMPLICITY(const Geometry &g);

namespace algorithm {

/**
 * Check simplicity of a geometry
 * @note When applied to NURBSCurve geometries, simplicity
 *   is internally checked on a LineString obtained via
 *   toLineString() with its default parameters.
 */
SFCGAL_API const Simplicity
isSimple(const Geometry &g, const double &toleranceAbs = 1e-9);

} // namespace algorithm
} // namespace SFCGAL

#endif
