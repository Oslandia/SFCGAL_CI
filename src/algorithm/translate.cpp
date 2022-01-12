// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/algorithm/translate.h>

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/MultiSolid.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>

#include <SFCGAL/detail/transform/AffineTransform2.h>
#include <SFCGAL/detail/transform/AffineTransform3.h>

namespace SFCGAL {
namespace algorithm {

///
///
///
void
translate(Geometry &g, const Kernel::Vector_3 &v)
{
  transform::AffineTransform3 visitor(
      CGAL::Aff_transformation_3<Kernel>(CGAL::TRANSLATION, v));
  g.accept(visitor);
}

///
///
///
void
translate(Geometry &g, const Kernel::Vector_2 &v)
{
  transform::AffineTransform2 visitor(
      CGAL::Aff_transformation_2<Kernel>(CGAL::TRANSLATION, v));
  g.accept(visitor);
}

///
///
///
void
translate(Geometry &g, Kernel::FT dx, Kernel::FT dy, Kernel::FT dz)
{
  translate(g, Kernel::Vector_3(dx, dy, dz));
}

} // namespace algorithm
} // namespace SFCGAL
