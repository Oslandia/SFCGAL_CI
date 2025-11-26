// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_DETAIL_FACEFILTERS_H
#define SFCGAL_ALGORITHM_DETAIL_FACEFILTERS_H

#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"

#include <algorithm>

namespace SFCGAL::algorithm::detail {

/**
 * @brief Check if a polygon face is not at a given height
 *
 * A face is considered to be at a given height if ALL its vertices
 * have a Z coordinate equal to that height.
 *
 * @param patch The polygon to check
 * @param height The reference height
 * @return true if the face is not at the given height (at least one vertex
 * differs), false otherwise
 */
inline auto
isNotFaceAtHeight(const Polygon &patch, double height) -> bool
{
  const LineString &exterior = patch.exteriorRing();

  // Check if all points have z == height
  bool allAtHeight =
      std::all_of(exterior.begin(), exterior.end(), [height](const Point &point) {
        return (point.z() - height) == 0;
      });

  // Return true if not all points are at height (keep all faces except those
  // at height)
  return !allAtHeight;
}

/**
 * @brief Check if a polygon face is not a base face (z = 0)
 *
 * A base face is a face where all vertices have Z = 0.
 * This function returns true if at least one vertex has Z != 0.
 *
 * @param patch The polygon to check
 * @return true if the face has at least one vertex with z != 0, false
 * otherwise
 */
inline auto
isNotBaseFace(const Polygon &patch) -> bool
{
  const LineString &exterior = patch.exteriorRing();

  // Check if any point has z != 0 (not a base face)
  return std::any_of(exterior.begin(), exterior.end(),
                     [](const Point &point) -> bool { return point.z() != 0; });
}

/**
 * @brief Check if a polygon face is a sloped roof face
 *
 * A roof slope is identified by having at least one vertex with Z != 0.
 * This is an alias for isNotBaseFace with a more semantic name for roof
 * contexts.
 *
 * @param patch The polygon to check
 * @return true if the face is a roof slope (has elevated points), false
 * otherwise
 */
inline auto
isRoofSlope(const Polygon &patch) -> bool
{
  return isNotBaseFace(patch);
}

} // namespace SFCGAL::algorithm::detail

#endif // SFCGAL_ALGORITHM_DETAIL_FACEFILTERS_H
