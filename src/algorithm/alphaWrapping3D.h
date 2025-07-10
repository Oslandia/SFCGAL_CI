// Copyright (c) 2024-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef _SFCGAL_ALGORITHM_ALPHAWRAPPING3D_H_
#define _SFCGAL_ALGORITHM_ALPHAWRAPPING3D_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/PolyhedralSurface.h"

namespace SFCGAL::algorithm {

/**
 * Computes the 3D alpha wrapping of a geometry
 * https://doc.cgal.org/latest/Alpha_wrap_3/index.html
 * @since 2.1
 * @param geom input geometry
 * @param relativeAlpha This parameter is used to determine which features will
 *  appear in the output A small relativeAlpha will produce an output less
 *  complex but less faithful to the input
 * @param relativeOffset  This parameter controls the tightness of the result
 *  A large relativeOffset parameter will tend to better preserve sharp features
 *  as projection If this parameter is equal to 0, it is computed from the alpha
 *  parameter
 * @return A PolyhedralSurface representing the 3D alpha wrapping of the
 * geometry
 */
SFCGAL_API std::unique_ptr<PolyhedralSurface>
           alphaWrapping3D(const Geometry &geom, size_t relativeAlpha,
                           size_t relativeOffset = 0);

} // namespace SFCGAL::algorithm

#endif
