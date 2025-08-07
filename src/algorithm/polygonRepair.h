// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_POLYGONREPAIR_H_
#define SFCGAL_ALGORITHM_POLYGONREPAIR_H_

#include "SFCGAL/Geometry.h"
#include <cstdint>
#include <memory>

namespace SFCGAL {

class Geometry;
class Polygon;
class MultiPolygon;

namespace algorithm {

/**
 * @brief Polygon repair rules available in CGAL 6.1
 * Note: EVEN_ODD_RULE is available in CGAL 6.0+
 * Other rules require CGAL 6.1+
 */
enum class PolygonRepairRule : std::uint8_t {
  EVEN_ODD_RULE     = 0, ///< Even-odd rule (default, available in CGAL 6.0+)
  NON_ZERO_RULE     = 1, ///< Non-zero winding rule (requires CGAL 6.1+)
  UNION_RULE        = 2, ///< Union of all polygons (requires CGAL 6.1+)
  INTERSECTION_RULE = 3  ///< Intersection of all polygons (requires CGAL 6.1+)
};

/**
 * @brief Repairs invalid polygons using CGAL's 2D Polygon Repair algorithm
 * CGAL 6.0 have only EVEN_ODD_RULE
 * Other rules have been added in 6.1(beta1)
 * https://doc.cgal.org/6.1-beta1/Polygon_repair/
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
 * ZM values are dropped.
 *
 * @param geometry input geometry (Polygon or MultiPolygon)
 * @param repairRule repair strategy to use (default: EVEN_ODD_RULE)
 * @return repaired geometry as MultiPolygon
 *
 * @post returned geometry is guaranteed to be valid
 */
SFCGAL_API auto
polygonRepair(const Geometry   &geometry,
              PolygonRepairRule repairRule = PolygonRepairRule::EVEN_ODD_RULE)
    -> std::unique_ptr<Geometry>;

} // namespace algorithm
} // namespace SFCGAL

#endif
