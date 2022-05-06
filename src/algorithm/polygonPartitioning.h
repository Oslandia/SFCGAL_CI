// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef _SFCGAL_ALGORITHM_POLYGONPARTITIONING_H_
#define _SFCGAL_ALGORITHM_POLYGONPARTITIONING_H_

#include <SFCGAL/config.h>

#include <SFCGAL/Geometry.h>

namespace SFCGAL {
namespace algorithm {

/**
 * Compute a partition of the geometry into convex polygons
 * @ingroup public_api
 */
SFCGAL_API std::unique_ptr<Geometry>
           greene_approx_convex_partition(const Geometry &g);

} // namespace algorithm
} // namespace SFCGAL

#endif
