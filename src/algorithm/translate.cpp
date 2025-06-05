// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/translate.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"

#include "SFCGAL/detail/transform/AffineTransform2.h"
#include "SFCGAL/detail/transform/AffineTransform3.h"

namespace SFCGAL::algorithm {

void
translate(Geometry &g, const Kernel::Vector_3 &v)
{
  transform::AffineTransform3 visitor(
      CGAL::Aff_transformation_3<Kernel>(CGAL::TRANSLATION, v));
  g.accept(visitor);
}

void
translate(Geometry &g, const Kernel::Vector_2 &v)
{
  transform::AffineTransform2 visitor(
      CGAL::Aff_transformation_2<Kernel>(CGAL::TRANSLATION, v));
  g.accept(visitor);
}

void
translate(Geometry &g, const Kernel::FT &dx, const Kernel::FT &dy,
          const Kernel::FT &dz)
{
  translate(g, Kernel::Vector_3(dx, dy, dz));
}

void
translate(Geometry &g, const double &dx, const double &dy, const double &dz)
{
  if (!std::isfinite(dx) || !std::isfinite(dy) || !std::isfinite(dz)) {
    BOOST_THROW_EXCEPTION(NonFiniteValueException(
        "trying to translate with non finite value in direction"));
  }

  return translate(g, Kernel::FT(dx), Kernel::FT(dy), Kernel::FT(dz));
}
} // namespace SFCGAL::algorithm
