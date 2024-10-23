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
#include "SFCGAL/config.h"

namespace SFCGAL::algorithm {

/**
 * Alpha shape calculation mode
 */
SFCGAL_API enum class AlphaShape3DMode {
  GENERAL,    ///< All faces from alpha complex (default)
  REGULARIZED ///< Only faces from regular alpha complex
};

/**
 * Compute the 3D alpha shape for a geometry with optimal alpha value
 * https://doc.cgal.org/latest/Alpha_shapes_3/index.html
 * @ingroup public_api
 * @since 2.1
 * @param geom input geometry
 * @param mode GENERAL or REGULARIZED mode
 * @return A PolyhedralSurface representing the optimal alpha shape
 */
SFCGAL_API std::unique_ptr<PolyhedralSurface>
           alphaShapes3D(const Geometry  &geom,
                         AlphaShape3DMode mode = AlphaShape3DMode::GENERAL);

} // namespace SFCGAL::algorithm

#endif
