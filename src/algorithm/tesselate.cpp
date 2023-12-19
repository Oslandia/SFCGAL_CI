// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/TriangulatedSurface.h>
#include <SFCGAL/algorithm/isValid.h>
#include <SFCGAL/algorithm/tesselate.h>
#include <SFCGAL/triangulate/triangulatePolygon.h>

namespace SFCGAL::algorithm {

///
///
///
auto
tesselate(const Geometry &g, NoValidityCheck /*unused*/)
    -> std::unique_ptr<Geometry>
{
  switch (g.geometryTypeId()) {
  case TYPE_POINT:
  case TYPE_LINESTRING:
  case TYPE_TRIANGLE:
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
    return std::unique_ptr<Geometry>(g.clone());

  case TYPE_POLYGON:
  case TYPE_POLYHEDRALSURFACE: {
    auto *triSurf = new TriangulatedSurface();
    triangulate::triangulatePolygon3D(g, *triSurf);
    return std::unique_ptr<Geometry>(triSurf);
  }

  case TYPE_SOLID: {
    std::unique_ptr<GeometryCollection> ret(new GeometryCollection);

    for (size_t i = 0; i < g.as<Solid>().numShells(); ++i) {
      const PolyhedralSurface &shellN = g.as<Solid>().shellN(i);

      if (!shellN.isEmpty()) {
        ret->addGeometry(tesselate(shellN).release());
      }
    }

    return std::unique_ptr<Geometry>(ret.release());
  }

  // multipolygon and multisolid return a geometrycollection
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION: {
    std::unique_ptr<GeometryCollection> ret(new GeometryCollection);

    for (size_t i = 0; i < g.numGeometries(); ++i) {
      ret->addGeometry(tesselate(g.geometryN(i)).release());
    }

    return std::unique_ptr<Geometry>(ret.release());
  }

  default:
    break;
  }

  return std::unique_ptr<Geometry>(g.clone());
}

auto
tesselate(const Geometry &g) -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY(g);

  return tesselate(g, NoValidityCheck());
}

} // namespace SFCGAL::algorithm
