/**
 *   SFCGAL
 *
 *   Copyright (C) 2024 SFCGAL Contributors and Oslandia <infos@oslandia.com>
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Library General Public
 *   License as published by the Free Software Foundation; either
 *   version 2 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Library General Public License for more details.

 *   You should have received a copy of the GNU Library General Public
 *   License along with this library; if not, see
 <http://www.gnu.org/licenses/>.
 */

#ifndef _SFCGAL_ALGORITHM_ALPHASHAPES3D_H_
#define _SFCGAL_ALGORITHM_ALPHASHAPES3D_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/PolyhedralSurface.h"

namespace SFCGAL::algorithm {

/**
 * Computes the 3D alpha wrapping of a geometry
 * https://doc.cgal.org/latest/Alpha_wrap_3/index.html
 * @ingroup public_api
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
