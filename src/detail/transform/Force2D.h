// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TRANSFORM_FORCE2D_H_
#define SFCGAL_TRANSFORM_FORCE2D_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Kernel.h"
#include "SFCGAL/Transform.h"

namespace SFCGAL {
namespace transform {

/**
 * Force 2D definitions
 */
class SFCGAL_API Force2D : public Transform {
public:
  /*
   * [SFCGAL::Transform]
   */
  void
  transform(Point &p) override;
};

} // namespace transform
} // namespace SFCGAL

#endif
