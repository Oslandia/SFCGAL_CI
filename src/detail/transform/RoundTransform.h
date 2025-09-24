// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TRANSFORM_ROUNDTRANSFORM_H_
#define SFCGAL_TRANSFORM_ROUNDTRANSFORM_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Kernel.h"
#include "SFCGAL/Transform.h"

namespace SFCGAL {
namespace transform {

/**
 * Round the coordinates of a Geometry
 */
class SFCGAL_API RoundTransform : public Transform {
public:
  /**
   * @brief Constructor with a scale factor (default is nearest integer)
   * @param scale The scale factor for rounding (default 1)
   */
  RoundTransform(const long &scale = 1);

  /**
   * @brief Transform a point by rounding its coordinates
   * @param p The point to transform
   */
  void
  transform(Point &p) override;

private:
  long _scale;
};

} // namespace transform
} // namespace SFCGAL

#endif
