// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/scale.h"
#include "SFCGAL/detail/transform/AffineTransform2.h"
#include "SFCGAL/detail/transform/AffineTransform3.h"
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Aff_transformation_3.h>

namespace SFCGAL::algorithm {

void
scale(Geometry &g, double s)
{
  scale(g, s, s, s);
}

void
scale(Geometry &g, double sx, double sy, double sz)
{
  if (g.is3D()) {
    CGAL::Aff_transformation_3<Kernel> transform(sx, 0, 0, 0, sy, 0, 0, 0, sz);
    transform::AffineTransform3        scaleTransform(transform);
    g.accept(scaleTransform);
  } else {
    CGAL::Aff_transformation_2<Kernel> transform(sx, 0, 0, sy);
    transform::AffineTransform2        scaleTransform(transform);
    g.accept(scaleTransform);
  }
}

void
scale(Geometry &g, double sx, double sy, double sz, double cx, double cy,
      double cz)
{
  if (g.is3D()) {
    CGAL::Aff_transformation_3<Kernel> toOrigin(
        CGAL::TRANSLATION, Kernel::Vector_3(-cx, -cy, -cz));
    CGAL::Aff_transformation_3<Kernel> scaleTransform(sx, 0, 0, 0, sy, 0, 0, 0,
                                                      sz);
    CGAL::Aff_transformation_3<Kernel> fromOrigin(CGAL::TRANSLATION,
                                                  Kernel::Vector_3(cx, cy, cz));

    CGAL::Aff_transformation_3<Kernel> combinedTransform =
        fromOrigin * scaleTransform * toOrigin;

    transform::AffineTransform3 affineTransform(combinedTransform);
    g.accept(affineTransform);
  } else {
    CGAL::Aff_transformation_2<Kernel> toOrigin(CGAL::TRANSLATION,
                                                Kernel::Vector_2(-cx, -cy));
    CGAL::Aff_transformation_2<Kernel> scaleTransform(sx, 0, 0, sy);
    CGAL::Aff_transformation_2<Kernel> fromOrigin(CGAL::TRANSLATION,
                                                  Kernel::Vector_2(cx, cy));

    CGAL::Aff_transformation_2<Kernel> combinedTransform =
        fromOrigin * scaleTransform * toOrigin;

    transform::AffineTransform2 affineTransform(combinedTransform);
    g.accept(affineTransform);
  }
}

} // namespace SFCGAL::algorithm
