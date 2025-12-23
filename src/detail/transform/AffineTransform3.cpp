// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/transform/AffineTransform3.h"

#include "SFCGAL/Transform.h"
#include <CGAL/Aff_transformation_3.h>

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

#include <utility>

namespace SFCGAL::transform {

AffineTransform3::AffineTransform3(CGAL::Aff_transformation_3<Kernel> transform)
    : _transform(std::move(transform))
{
}

void
AffineTransform3::transform(Point &point)
{
  if (!point.isEmpty()) {
    Point pt(point.toPoint_3().transform(_transform));
    if (point.isMeasured()) {
      pt.setM(point.m());
    }
    point = pt;
  }
}

void
AffineTransform3::transform(LineString &linestring)
{
  if (linestring.isEmpty()) {
    return;
  }

  for (size_t i = 0; i < linestring.numPoints(); ++i) {
    transform(linestring.pointN(i));
  }
}

void
AffineTransform3::transform(Triangle &triangle)
{
  if (triangle.isEmpty()) {
    return;
  }

  transform(triangle.vertex(0));
  transform(triangle.vertex(1));
  transform(triangle.vertex(2));
}

void
AffineTransform3::transform(Polygon &polygon)
{
  if (polygon.isEmpty()) {
    return;
  }

  transform(polygon.exteriorRing());

  for (size_t i = 0; i < polygon.numInteriorRings(); ++i) {
    transform(polygon.interiorRingN(i));
  }
}

void
AffineTransform3::transform(PolyhedralSurface &surface)
{
  for (size_t i = 0; i < surface.numPatches(); ++i) {
    transform(surface.patchN(i));
  }
}

void
AffineTransform3::transform(TriangulatedSurface &surface)
{
  for (size_t i = 0; i < surface.numPatches(); ++i) {
    transform(surface.patchN(i));
  }
}

void
AffineTransform3::transform(Solid &solid)
{
  if (solid.isEmpty()) {
    return;
  }

  transform(solid.exteriorShell());

  for (size_t i = 0; i < solid.numInteriorShells(); ++i) {
    transform(solid.interiorShellN(i));
  }
}

} // namespace SFCGAL::transform
