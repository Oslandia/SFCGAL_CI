// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TRANSFORM_FORCEM_H_
#define SFCGAL_TRANSFORM_FORCEM_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Transform.h"

namespace SFCGAL::transform {

/**
 * Force M definitions
 */
class SFCGAL_API ForceM : public Transform {
public:
  /**
   * Constructor with a default M value
   */
  ForceM(const double &defaultM = 0);

  /*
   * [SFCGAL::Transform]
   */
  void
  transform(Point &point) override;

private:
  double _defaultM;
};

} // namespace SFCGAL::transform

#endif
