// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef _SFCGAL_ALGORITHM_MAKEBUILDING_H_
#define _SFCGAL_ALGORITHM_MAKEBUILDING_H_

#include <SFCGAL/Geometry.h>
#include <SFCGAL/Kernel.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>

namespace SFCGAL {
namespace algorithm {
struct NoValidityCheck;

/**
 * Extrude a polygon with a fake roof top
 *
 * @ingroup public_api
 * @warning Z component is ignored.
 * @ingroup public_api
 * @pre g is a valid polygon
 */
SFCGAL_API auto
makebuilding(const Geometry &g, double height_building, double height_roof) -> std::unique_ptr<SFCGAL::PolyhedralSurface>;

/**
 * Extrude a polygon with a fake roof top
 *
 * @ingroup public_api
 * @warning Z component is ignored.
 * @ingroup public_api
 * @pre g is a valid polygon
 */
SFCGAL_API auto
makebuilding(const Geometry &g, Kernel::FT &height_building, Kernel::FT &height_roof) -> std::unique_ptr<SFCGAL::PolyhedralSurface>;

/**
 * Extrude a polygon with a fake roof top
 *
 * @ingroup detail
 * @pre g is a valid geometry
 * @warning No actual validity check is done
 */
SFCGAL_API auto
makebuilding(const Geometry &g, Kernel::FT &height_building, Kernel::FT &height_roof, NoValidityCheck) -> std::unique_ptr<SFCGAL::PolyhedralSurface>;

/**
 * Returns Compute the 2D area for a Polygon
 * @ingroup detail
 */
SFCGAL_API auto
makebuilding(const Polygon &g, Kernel::FT &height_building, Kernel::FT &height_roof) -> std::unique_ptr<SFCGAL::PolyhedralSurface>;

} // namespace algorithm
} // namespace SFCGAL
#endif 
