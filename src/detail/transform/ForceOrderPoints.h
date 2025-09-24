// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TRANSFORM_FORCE_ORDER_POINTS_H_
#define SFCGAL_TRANSFORM_FORCE_ORDER_POINTS_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Kernel.h"
#include "SFCGAL/Transform.h"

namespace SFCGAL {

namespace transform {

/**
 * If the 2D surface is pointing down, reverse its points
 * @todo unittest
 * @todo move outside (it's not a coordinate transform)?
 */
class SFCGAL_API ForceOrderPoints : public Transform {
public:
  /**
   * @brief Pass the forced orientation as parameter
   * @param orientCCW Whether to force counter-clockwise orientation (default
   * true)
   */
  ForceOrderPoints(bool orientCCW = true);

  /**
   * @brief Transform a point (no-op for points)
   * @param p The point to transform
   */
  void
  transform(Point &p) override;

  /**
   * @brief Visit and force point order for a triangle
   * @param triangle The triangle to visit
   */
  void
  visit(Triangle &triangle) override;
  /**
   * @brief Visit and force point order for a polygon
   * @param polygon The polygon to visit
   */
  void
  visit(Polygon &polygon) override;

private:
  bool _orientCCW;
};

} // namespace transform
} // namespace SFCGAL

#endif
