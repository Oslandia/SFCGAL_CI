// Copyright (c) 2024-2024, SFCGAL Contributors and Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later
#include "SFCGAL/algorithm/rotate.h"
#include "SFCGAL/detail/transform/AffineTransform2.h"
#include "SFCGAL/detail/transform/AffineTransform3.h"
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Aff_transformation_3.h>
#include <cmath>

namespace SFCGAL::algorithm {

void
rotate(Geometry &g, const Kernel::FT &angle)
{
  double                       s = std::sin(CGAL::to_double(angle));
  double                       c = std::cos(CGAL::to_double(angle));
  Kernel::Aff_transformation_2 rot(c, -s, s, c);
  transform::AffineTransform2  visitor(rot);
  g.accept(visitor);
}

void
rotate(Geometry &g, const Kernel::FT &angle, const Point &origin)
{
  double                       s  = std::sin(CGAL::to_double(angle));
  double                       c  = std::cos(CGAL::to_double(angle));
  Kernel::FT                   ox = origin.x();
  Kernel::FT                   oy = origin.y();
  Kernel::Aff_transformation_2 toOrigin(CGAL::TRANSLATION,
                                        Kernel::Vector_2(-ox, -oy));
  Kernel::Aff_transformation_2 rot(c, -s, s, c);
  Kernel::Aff_transformation_2 fromOrigin(CGAL::TRANSLATION,
                                          Kernel::Vector_2(ox, oy));
  Kernel::Aff_transformation_2 combined = fromOrigin * rot * toOrigin;
  transform::AffineTransform2  visitor(combined);
  g.accept(visitor);
}

void
rotate(Geometry &g, const Kernel::FT &angle, const Kernel::Vector_3 &axis,
       const Point &origin)
{
  Kernel::Vector_3 u = axis / std::sqrt(CGAL::to_double(axis.squared_length()));
  double           theta = CGAL::to_double(angle);
  double           c     = std::cos(theta);
  double           s     = std::sin(theta);
  double           t     = 1 - c;

  double ux = CGAL::to_double(u.x());
  double uy = CGAL::to_double(u.y());
  double uz = CGAL::to_double(u.z());

  Kernel::Aff_transformation_3 rot(
      t * ux * ux + c, t * ux * uy - s * uz, t * ux * uz + s * uy,
      t * ux * uy + s * uz, t * uy * uy + c, t * uy * uz - s * ux,
      t * ux * uz - s * uy, t * uy * uz + s * ux, t * uz * uz + c);

  Kernel::FT                   ox = origin.x();
  Kernel::FT                   oy = origin.y();
  Kernel::FT                   oz = origin.z();
  Kernel::Aff_transformation_3 toOrigin(CGAL::TRANSLATION,
                                        Kernel::Vector_3(-ox, -oy, -oz));
  Kernel::Aff_transformation_3 fromOrigin(CGAL::TRANSLATION,
                                          Kernel::Vector_3(ox, oy, oz));

  Kernel::Aff_transformation_3 combined = fromOrigin * rot * toOrigin;
  transform::AffineTransform3  visitor(combined);
  g.accept(visitor);
}

void
rotateX(Geometry &g, const Kernel::FT &angle)
{
  double                       s = std::sin(CGAL::to_double(angle));
  double                       c = std::cos(CGAL::to_double(angle));
  Kernel::Aff_transformation_3 rot(1, 0, 0, 0, c, -s, 0, s, c);
  transform::AffineTransform3  visitor(rot);
  g.accept(visitor);
}

void
rotateY(Geometry &g, const Kernel::FT &angle)
{
  double                       s = std::sin(CGAL::to_double(angle));
  double                       c = std::cos(CGAL::to_double(angle));
  Kernel::Aff_transformation_3 rot(c, 0, s, 0, 1, 0, -s, 0, c);
  transform::AffineTransform3  visitor(rot);
  g.accept(visitor);
}

void
rotateZ(Geometry &g, const Kernel::FT &angle)
{
  double                       s = std::sin(CGAL::to_double(angle));
  double                       c = std::cos(CGAL::to_double(angle));
  Kernel::Aff_transformation_3 rot(c, -s, 0, s, c, 0, 0, 0, 1);
  transform::AffineTransform3  visitor(rot);
  g.accept(visitor);
}

} // namespace SFCGAL::algorithm
