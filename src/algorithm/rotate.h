// Copyright (c) 2024-2024, SFCGAL Contributors and Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_ROTATE_H_
#define SFCGAL_ALGORITHM_ROTATE_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/export.h"

namespace SFCGAL::algorithm {

/**
 * @brief Rotate a geometry in 2D around the origin (0,0)
 * @param geometry The geometry to rotate
 * @param angle Rotation angle in radians
 */
SFCGAL_API void
rotate(Geometry &geometry, const Kernel::FT &angle);

/**
 * @brief Rotate a geometry in 2D around a specified point
 * @param geometry The geometry to rotate
 * @param angle Rotation angle in radians
 * @param origin Point of origin for the rotation
 */
SFCGAL_API void
rotate(Geometry &geometry, const Kernel::FT &angle, const Point &origin);

/**
 * @brief Rotate a geometry in 3D around a specified axis and origin
 * @param geometry The geometry to rotate
 * @param angle Rotation angle in radians
 * @param axis The axis of rotation
 * @param origin Point of origin for the rotation
 */
SFCGAL_API void
rotate(Geometry &geometry, const Kernel::FT &angle,
       const Kernel::Vector_3 &axis, const Point &origin = Point(0, 0, 0));

/**
 * @brief Rotate a geometry around the X axis
 * @param geometry The geometry to rotate
 * @param angle Rotation angle in radians
 */
SFCGAL_API void
rotateX(Geometry &geometry, const Kernel::FT &angle);

/**
 * @brief Rotate a geometry around the Y axis
 * @param geometry The geometry to rotate
 * @param angle Rotation angle in radians
 */
SFCGAL_API void
rotateY(Geometry &geometry, const Kernel::FT &angle);

/**
 * @brief Rotate a geometry around the Z axis
 * @param geometry The geometry to rotate
 * @param angle Rotation angle in radians
 */
SFCGAL_API void
rotateZ(Geometry &geometry, const Kernel::FT &angle);

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_ROTATE_H_
