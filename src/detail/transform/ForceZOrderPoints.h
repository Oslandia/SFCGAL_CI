// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TRANSFORM_FORCEZ_ORDER_POINTS_H_
#define SFCGAL_TRANSFORM_FORCEZ_ORDER_POINTS_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Kernel.h"
#include "SFCGAL/Transform.h"

namespace SFCGAL {

namespace transform {

/**
 * Force Z. If the 2D surface is pointing down, reverse its points
 *
 * @todo unittest
 * @todo move outside (it's not a coordinate transform)?
 */
class SFCGAL_API ForceZOrderPoints : public Transform {
public:
  /**
   * Constructor with a default Z value
   */
  ForceZOrderPoints(const Kernel::FT defaultZ = 0);

  /*
   * [SFCGAL::Transform]
   */
  void
  transform(Point &p) override;

  void
  visit(Triangle &) override;
  void
  visit(Polygon &) override;

private:
  Kernel::FT _defaultZ;
};

} // namespace transform
} // namespace SFCGAL

#endif
