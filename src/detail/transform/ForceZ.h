// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TRANSFORM_FORCEZ_H_
#define SFCGAL_TRANSFORM_FORCEZ_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Kernel.h"
#include "SFCGAL/Transform.h"

namespace SFCGAL {
namespace transform {

/**
 * Force Z definitions
 */
class SFCGAL_API ForceZ : public Transform {
public:
  /**
   * Constructor with a default Z value
   */
  ForceZ(const Kernel::FT defaultZ = 0);

  /*
   * [SFCGAL::Transform]
   */
  void
  transform(Point &p) override;

private:
  Kernel::FT _defaultZ;
};

} // namespace transform
} // namespace SFCGAL

#endif
