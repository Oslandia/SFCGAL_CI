// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/volume.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/algorithm/tesselate.h"

namespace SFCGAL::algorithm {

auto
volume(const Solid &solid, NoValidityCheck /*unused*/) -> Kernel::FT
{
  Kernel::FT                  vol = 0;
  const CGAL::Point_3<Kernel> origin(0, 0, 0);
  const size_t                numShells = solid.numShells();

  for (size_t i = 0; i < numShells; i++) {
    std::unique_ptr<Geometry>  t(tesselate(solid.shellN(i), NoValidityCheck()));
    const TriangulatedSurface &tin        = t->as<TriangulatedSurface>();
    const size_t               numPatches = tin.numPatches();

    for (size_t j = 0; j < numPatches; j++) {
      const Triangle &tri = tin.patchN(j);
      vol = vol + CGAL::volume(origin, tri.vertex(0).toPoint_3(),
                               tri.vertex(1).toPoint_3(),
                               tri.vertex(2).toPoint_3());
    }
  }

  return vol;
}

auto
volume(const Geometry &g) -> Kernel::FT
{
  if (g.isEmpty()) {
    return 0;
  }

  SFCGAL_ASSERT_GEOMETRY_VALIDITY(g);

  switch (g.geometryTypeId()) {
  case TYPE_POINT:
  case TYPE_LINESTRING:
  case TYPE_POLYGON:
  case TYPE_TRIANGLE:
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
    return 0;

  case TYPE_SOLID:
    return volume(g.as<Solid>(), NoValidityCheck());

  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION:
    Kernel::FT  v = 0;
    const auto &c = g.as<GeometryCollection>();

    for (size_t i = 0; i < c.numGeometries(); i++) {
      if (c.geometryN(i).is<Solid>()) {
        v = v + volume(c.geometryN(i).as<Solid>(), NoValidityCheck());
      }
    }

    return v;
  }

  BOOST_THROW_EXCEPTION(Exception(
      (boost::format("volume( %s ) is not defined") % g.geometryType()).str()));
  return 0; // to avoid warning
}

} // namespace SFCGAL::algorithm
