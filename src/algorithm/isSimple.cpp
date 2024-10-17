// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/isSimple.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"

#include "SFCGAL/Kernel.h"
#include "SFCGAL/algorithm/intersects.h"
#include "SFCGAL/algorithm/plane.h"
#include "SFCGAL/detail/algorithm/coversPoints.h"

using namespace SFCGAL::detail::algorithm;

namespace SFCGAL {

namespace algorithm {

auto
isSimple(const LineString &linestring) -> const Simplicity
{
  if (linestring.is3D() ? selfIntersects3D(linestring)
                        : selfIntersects(linestring)) {
    return Simplicity::complex(
        boost::format("linestring self intersects").str());
  }
  return Simplicity::simple();
}

auto
isSimple(const Polygon &polygon,
         const double  &toleranceAbs = 1e-9) -> const Simplicity
{
  // Polygone must be planar (all points in the same plane)
  if (polygon.is3D() && !isPlane3D<Kernel>(polygon, toleranceAbs)) {
    return Simplicity::complex("Points don't lie in the same plane.");
  }
  return Simplicity::simple();
}

auto
isSimple(const MultiPoint &multipoint) -> const Simplicity
{
  const size_t numPoint = multipoint.numGeometries();

  for (size_t l = 0; l != numPoint - 1; ++l) {
    for (size_t other_l = l + 1; other_l != numPoint; ++other_l) {
      bool const duplicated_points =
          multipoint.pointN(l) == multipoint.pointN(other_l);
      if (duplicated_points) {
        return Simplicity::complex(
            (boost::format(
                 "Points %d and %d are duplicated in the MultiPoint.") %
             l % other_l)
                .str());
      }
    }
  }

  return Simplicity::simple();
}

auto
isSimple(const Geometry &g, const double &toleranceAbs) -> const Simplicity
{
  switch (g.geometryTypeId()) {
  case TYPE_POINT:
    return Simplicity::simple();

  case TYPE_LINESTRING: // no self-intersecting (excepted at ending point)
    return isSimple(g.as<LineString>());

  case TYPE_POLYGON: // are the rings simple? (3D) -> check the coplanarity
    return isSimple(g.as<Polygon>(), toleranceAbs);

  case TYPE_TRIANGLE:
    return Simplicity::simple();

  case TYPE_SOLID: // every phs is simple
    return isSimple(g.as<Solid>());

  case TYPE_MULTIPOINT: // no equal points
    return isSimple(g.as<MultiPoint>());

  case TYPE_MULTILINESTRING: // every ls is simple, and only intersections are
                             // at boundaries
    return isSimple(g.as<MultiLineString>());

  case TYPE_MULTIPOLYGON: // every polygon is simple
    return isSimple(g.as<MultiPolygon>());

  case TYPE_MULTISOLID: // every solid is simple
    return isSimple(g.as<MultiSolid>());

  case TYPE_GEOMETRYCOLLECTION: // every geometry is simple
    return isSimple(g.as<GeometryCollection>());

  case TYPE_TRIANGULATEDSURFACE: // every triangle is simple
    return isSimple(g.as<TriangulatedSurface>());

  case TYPE_POLYHEDRALSURFACE: // every polygon is simple
    return isSimple(g.as<PolyhedralSurface>());
  }

  BOOST_THROW_EXCEPTION(Exception(
      (boost::format("isSimple( %s ) is not defined") % g.geometryType())
          .str()));
  return Simplicity::complex(
      (boost::format("isSimple( %s ) is not defined") % g.geometryType())
          .str()); // to avoid warning
}

} // namespace algorithm
} // namespace SFCGAL
