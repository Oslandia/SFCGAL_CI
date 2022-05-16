// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef _SFCGAL_ALGORITHM_POLYGONPARTITIONING_H_
#define _SFCGAL_ALGORITHM_POLYGONPARTITIONING_H_

#include <SFCGAL/config.h>

#include <SFCGAL/Geometry.h>

namespace SFCGAL {
namespace algorithm {
struct NoValidityCheck;

/**
 * Compute a partition of a surface geometry (mulitpolygon/polygon)
 * into convex polygons. Return the original geometry otherwise.
 * @pre g is a valid geometry
 * @ingroup public_api
 */
SFCGAL_API std::unique_ptr<Geometry>
           greene_approx_convex_partition(const Geometry &g);

/**
 * Compute a partition of a surface geometry (mulitpolygon/polygon)
 * into convex polygons. Return the original geometry otherwise.
 * @pre g is a valid geometry
 * @ingroup detail
 * @warning No actual validity check is done.
 */
SFCGAL_API std::unique_ptr<Geometry>
           greene_approx_convex_partition(const Geometry &g, NoValidityCheck);

} // namespace algorithm
} // namespace SFCGAL

#endif
