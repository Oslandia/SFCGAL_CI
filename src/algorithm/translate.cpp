// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
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
translate(Geometry &geometry, const Kernel::Vector_3 &vector)
{
  transform::AffineTransform3 visitor(
      CGAL::Aff_transformation_3<Kernel>(CGAL::TRANSLATION, vector));
  geometry.accept(visitor);
}

void
translate(Geometry &geometry, const Kernel::Vector_2 &vector)
{
  transform::AffineTransform2 visitor(
      CGAL::Aff_transformation_2<Kernel>(CGAL::TRANSLATION, vector));
  geometry.accept(visitor);
}

void
translate(Geometry &geometry, const Kernel::FT &dx, const Kernel::FT &dy,
          const Kernel::FT &dz)
{
  translate(geometry, Kernel::Vector_3(dx, dy, dz));
}

void
translate(Geometry &geometry, const double &dx, const double &dy,
          const double &dz)
{
  if (!std::isfinite(dx) || !std::isfinite(dy) || !std::isfinite(dz)) {
    BOOST_THROW_EXCEPTION(NonFiniteValueException(
        "trying to translate with non finite value in direction"));
  }

  translate(geometry, Kernel::FT(dx), Kernel::FT(dy), Kernel::FT(dz));
}
} // namespace SFCGAL::algorithm
