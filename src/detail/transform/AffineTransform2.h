// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TRANSFORM_AFFINETRANSFORM2_H_
#define SFCGAL_TRANSFORM_AFFINETRANSFORM2_H_

#include "SFCGAL/Kernel.h"
#include "SFCGAL/Transform.h"
#include "SFCGAL/config.h"

#include <CGAL/Aff_transformation_2.h>

namespace SFCGAL {
namespace transform {

/**
 * Wrapper for CGAL::Aff_transform_2
 * @todo unittest
 */
class SFCGAL_API AffineTransform2 : public Transform {
public:
  /**
   * Constructor with a transform
   */
  AffineTransform2(CGAL::Aff_transformation_2<Kernel> transform);

  /*
   * [SFCGAL::Transform]
   */
  void
  transform(Point &p) override;

private:
  CGAL::Aff_transformation_2<Kernel> _transform;
};

} // namespace transform
} // namespace SFCGAL

#endif
