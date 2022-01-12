// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/detail/transform/AffineTransform3.h>

#include <CGAL/Aff_transformation_3.h>
#include <SFCGAL/Transform.h>

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

namespace SFCGAL {
namespace transform {

///
///
///
AffineTransform3::AffineTransform3(CGAL::Aff_transformation_3<Kernel> transform)
    : _transform(transform)
{
}

///
///
///
void
AffineTransform3::transform(Point &p)
{
  if (!p.isEmpty()) {
    Point pt(p.toPoint_3().transform(_transform));
    if (p.isMeasured())
      pt.setM(p.m());
    p = pt;
  }
}

///
///
///
void
AffineTransform3::transform(LineString &ls)
{
  for (size_t i = 0; i < ls.numPoints(); ++i) {
    transform(ls.pointN(i));
  }
}

///
///
///
void
AffineTransform3::transform(Triangle &tri)
{
  transform(tri.vertex(0));
  transform(tri.vertex(1));
  transform(tri.vertex(2));
}

///
///
///
void
AffineTransform3::transform(Polygon &poly)
{
  transform(poly.exteriorRing());

  for (size_t i = 0; i < poly.numInteriorRings(); ++i) {
    transform(poly.interiorRingN(i));
  }
}

///
///
///
void
AffineTransform3::transform(PolyhedralSurface &surf)
{
  for (size_t i = 0; i < surf.numPolygons(); ++i) {
    transform(surf.polygonN(i));
  }
}

///
///
///
void
AffineTransform3::transform(TriangulatedSurface &surf)
{
  for (size_t i = 0; i < surf.numGeometries(); ++i) {
    transform(surf.geometryN(i));
  }
}

///
///
///
void
AffineTransform3::transform(Solid &solid)
{
  transform(solid.exteriorShell());

  for (size_t i = 0; i < solid.numInteriorShells(); ++i) {
    transform(solid.interiorShellN(i));
  }
}

} // namespace transform
} // namespace SFCGAL
