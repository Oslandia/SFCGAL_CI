// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/transform/AffineTransform2.h"

#include "SFCGAL/Point.h"

#include <utility>

namespace SFCGAL::transform {

AffineTransform2::AffineTransform2(CGAL::Aff_transformation_2<Kernel> transform)
    : _transform(std::move(transform))
{
}

/*
 * [SFCGAL::Transform]
 */
void
AffineTransform2::transform(Point &p)
{
  if (!p.isEmpty()) {
    // FIXME add a Point::Point( Kernel::Point_2&, double m );
    Point pt(p.toPoint_2().transform(_transform));
    if (p.isMeasured()) {
      pt.setM(p.m());
    }
    p = pt;
  }
}

} // namespace SFCGAL::transform
