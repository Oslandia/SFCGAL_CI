// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/volume.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/algorithm/isValid.h"

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

namespace SFCGAL::algorithm {

auto
volume(const PolyhedralSurface &surface, NoValidityCheck /*unused*/)
    -> Kernel::FT
{
  // Convert to CGAL Polyhedron
  auto cgal_polyhedron =
      surface.toPolyhedron_3<SFCGAL::Kernel, detail::MarkedPolyhedron>();

  // Ensure correct orientation
  // TODO: Could be an option or a dedicated algorithm
  // if (!CGAL::Polygon_mesh_processing::is_outward_oriented(*cgal_polyhedron))
  // {
  //   CGAL::Polygon_mesh_processing::orient_to_bound_a_volume(*cgal_polyhedron);
  // }

  // Calculate and return the volume
  return CGAL::Polygon_mesh_processing::volume(*cgal_polyhedron);
}

auto
volume(const Solid &solid, NoValidityCheck /*unused*/) -> Kernel::FT
{
  Kernel::FT   vol       = 0;
  const size_t numShells = solid.numShells();

  // A Solid is composed of PolyhedralSurfaces (shells)
  // The outer shell contributes positively to the volume
  // Inner shells (holes) contribute negatively
  for (size_t i = 0; i < numShells; i++) {
    const PolyhedralSurface &shell = solid.shellN(i);

    // Calculate volume of this shell
    Kernel::FT shell_volume = volume(shell, NoValidityCheck());

    // First shell (outer) is positive, others (holes) are negative
    if (i == 0) {
      vol = shell_volume;
    } else {
      // Inner shells represent cavities, so subtract their volume
      vol = vol - CGAL::abs(shell_volume);
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
  case TYPE_NURBSCURVE:
  case TYPE_POLYGON:
  case TYPE_TRIANGLE:
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_TRIANGULATEDSURFACE:
    return 0;

  case TYPE_POLYHEDRALSURFACE:
    return volume(g.as<PolyhedralSurface>(), NoValidityCheck());

  case TYPE_SOLID:
    return volume(g.as<Solid>(), NoValidityCheck());

  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION: {
    Kernel::FT  v = 0;
    const auto &c = g.as<GeometryCollection>();

    for (size_t i = 0; i < c.numGeometries(); i++) {
      if (c.geometryN(i).is<Solid>()) {
        v = v + volume(c.geometryN(i).as<Solid>(), NoValidityCheck());
      } else if (c.geometryN(i).is<PolyhedralSurface>()) {
        // Also handle PolyhedralSurfaces in collections
        const auto &surf = c.geometryN(i).as<PolyhedralSurface>();
        v                = v + volume(surf, NoValidityCheck());
      }
    }

    return v;
  }
  }

  BOOST_THROW_EXCEPTION(Exception(
      (boost::format("volume( %s ) is not defined") % g.geometryType()).str()));
  return 0; // to avoid warning
}

} // namespace SFCGAL::algorithm
