// Copyright (c) 2012-2024, SFCGAL Contributors and Oslandia
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_MINKOWSKISUM3D_H_
#define SFCGAL_ALGORITHM_MINKOWSKISUM3D_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/config.h"
#include <CGAL/Nef_polyhedron_3.h>
#include <memory>

namespace SFCGAL::algorithm {

struct NoValidityCheck;

/**
 * @brief 3D Minkowski sum (p+q)
 * @param gA the first geometry
 * @param gB the second geometry
 * @return the 3D Minkowski sum as a unique_ptr<Geometry>
 * @pre gA and gB are valid 3D geometries
 */
SFCGAL_API std::unique_ptr<Geometry>
           minkowskiSum3D(const Geometry &gA, const Geometry &gB);

/**
 * @brief 3D Minkowski sum (p+q)
 * @param gA the first geometry
 * @param gB the second geometry
 * @return the 3D Minkowski sum as a unique_ptr<Geometry>
 * @pre gA and gB are valid 3D geometries
 * @warning No actual validity check is done.
 */
SFCGAL_API std::unique_ptr<Geometry>
minkowskiSum3D(const Geometry &gA, const Geometry &gB, NoValidityCheck);

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_MINKOWSKISUM3D_H_
