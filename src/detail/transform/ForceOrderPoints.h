// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
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
   * Pass the forced orientation as parameter
   */
  ForceOrderPoints(bool orientCCW = true);

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
  bool _orientCCW;
};

} // namespace transform
} // namespace SFCGAL

#endif
