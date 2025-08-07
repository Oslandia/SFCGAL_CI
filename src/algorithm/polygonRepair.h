// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_POLYGONREPAIR_H_
#define SFCGAL_ALGORITHM_POLYGONREPAIR_H_

#include <SFCGAL/config.h>
#include <cstdint>
#include <memory>

namespace SFCGAL {

class Geometry;
class Polygon;
class MultiPolygon;

namespace algorithm {

/**
 * @brief Polygon repair rules available in CGAL 6.1
 * @ingroup algorithms
 */
enum class PolygonRepairRule : std::uint8_t {
  EVEN_ODD_RULE = 0, ///< Even-odd rule (default)
#if CGAL_VERSION_MAJOR == 6 && CGAL_VERSION_MINOR >= 1
  NON_ZERO_RULE     = 1, ///< Non-zero winding rule
  UNION_RULE        = 2, ///< Union of all polygons
  INTERSECTION_RULE = 3  ///< Intersection of all polygons
#endif
};

/**
 * @brief Repairs invalid polygons using CGAL's 2D Polygon Repair algorithm
 * https://doc.cgal.org/latest/Polygon_repair/index.html
 * CGAL 6.0 have only EVEN_ODD_RULE
 * Other rules have been added in 6.1(beta1)
 * //doc.cgal.org/6.1-beta1/Polygon_repair/index.html
 *
 * This function takes potentially invalid polygon geometries
 * (self-intersecting, self-touching, badly nested holes, etc.) and returns
 * valid polygon(s) using CGAL 6.1's robust polygon repair implementation.
 *
 * The algorithm handles:
 * - Self-intersecting polygons
 * - Self-touching polygons
 * - Invalid hole configurations
 * - Overlapping polygons in multipolygons
 * - Polygons touching at edges
 * - Incorrect polygon orientations
 *
 * All processing is done using CGAL's Multipolygon_with_holes_2 to uniformly
 * handle single polygons, polygons with holes, and multipolygons.
 *
 * @param geometry input geometry (Polygon or MultiPolygon)
 * @param repairRule repair strategy to use (default: EVEN_ODD_RULE)
 * @return repaired geometry as MultiPolygon
 *
 * @post returned geometry is guaranteed to be valid
 *
 * @ingroup algorithms
 */
SFCGAL_API auto
polygonRepair(const Geometry   &geometry,
              PolygonRepairRule repairRule = PolygonRepairRule::EVEN_ODD_RULE)
    -> std::unique_ptr<Geometry>;

} // namespace algorithm
} // namespace SFCGAL

#endif
