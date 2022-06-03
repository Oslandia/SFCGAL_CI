// Copyright (c) 2MultiLineString-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "force3D.h"
#include <SFCGAL/algorithm/extrude.h>
#include <SFCGAL/algorithm/force3D.h>
#include <SFCGAL/algorithm/makebuilding.h>
#include <SFCGAL/algorithm/straightSkeleton.h>
#include <SFCGAL/algorithm/tesselate.h>
#include <SFCGAL/algorithm/union.h>
#include <SFCGAL/triangulate/triangulate2DZ.h>

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>

#include <SFCGAL/algorithm/isValid.h>
#include <SFCGAL/algorithm/plane.h>

#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>

#include <SFCGAL/Exception.h>
#include <boost/format.hpp>
#include <type_traits>

namespace SFCGAL {
namespace algorithm {

using Point_2    = CGAL::Point_2<SFCGAL::Kernel>;
using Triangle_2 = CGAL::Triangle_2<SFCGAL::Kernel>;
using Polygon_2  = CGAL::Polygon_2<SFCGAL::Kernel>;

///
///
///
auto
makebuilding(const Geometry &g, Kernel::FT &height_building,
             Kernel::FT &height_roof, NoValidityCheck)
    -> std::unique_ptr<SFCGAL::PolyhedralSurface>
{
  switch (g.geometryTypeId()) {
  case TYPE_POINT:
  case TYPE_LINESTRING:
  case TYPE_TRIANGLE:
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
  case TYPE_SOLID:
  case TYPE_MULTISOLID:
  case TYPE_MULTIPOLYGON:
    return 0;

  case TYPE_POLYGON:
    return makebuilding(g.as<Polygon>(), height_building, height_roof);
  }

  BOOST_THROW_EXCEPTION(Exception(
      (boost::format(
           "Unexpected geometry type (%s) in SFCGAL::algorithm::makebuilding") %
       g.geometryType())
          .str()));
}

auto
makebuilding(const Geometry &g, Kernel::FT &height_building,
             Kernel::FT &height_roof)
    -> std::unique_ptr<SFCGAL::PolyhedralSurface>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(g);
  return makebuilding(g, height_building, height_roof, NoValidityCheck());
}

auto
makebuilding(const Geometry &g, double height_building, double height_roof)
    -> std::unique_ptr<SFCGAL::PolyhedralSurface>
{

  Kernel::FT hbuild(height_building);
  Kernel::FT hroof(height_roof);
  return makebuilding(g, hbuild, hroof);
}

auto
makebuilding(const Polygon &polygon, Kernel::FT &height_building,
             Kernel::FT &height_roof)
    -> std::unique_ptr<SFCGAL::PolyhedralSurface>
{
  auto sk{straightSkeleton(polygon)};
  for (const auto &line : polygon) {
    sk->addGeometry(line);
  }
  MultiLineString multilinestring;
  for (int i = 0; i < sk->numGeometries(); ++i) {
    const LineString line = sk->lineStringN(i);
    LineString       out;
    out.reserve(line.numPoints());
    for (int v = 0; v < line.numPoints(); ++v) {
      Point      p{line.pointN(v)};
      Kernel::FT z{height_building};
      if (polygon.toPolygon_2().bounded_side(p.toPoint_2()) ==
          CGAL::ON_BOUNDED_SIDE) {
        z += height_roof;
      }
      out.addPoint(Point(p.x(), p.y(), z, p.m()));
    }
    multilinestring.addGeometry(out);
  }

  auto *surf          = new SFCGAL::TriangulatedSurface;
  auto *building_surf = new SFCGAL::TriangulatedSurface;
  SFCGAL::triangulate::ConstraintDelaunayTriangulation cdt;
  SFCGAL::triangulate::triangulate2DZ(multilinestring, cdt);
  cdt.getTriangles(*surf);

  auto extruded{extrude(polygon.exteriorRing(), 0, 0, height_building)};
  // TODO: with hole
  // SOLIDZ+TINZ ?

  SFCGAL::triangulate::ConstraintDelaunayTriangulation building_cdt;
  SFCGAL::triangulate::triangulate2DZ(*extruded, building_cdt);
  building_cdt.getTriangles(*building_surf);

  for (int i = 0; i < surf->numGeometries(); ++i) {
    extruded->as<PolyhedralSurface>().addPolygon(surf->triangleN(i));
  }

  return std::make_unique<PolyhedralSurface>(
      *(extruded->as<PolyhedralSurface>().clone()));
}

} // namespace algorithm
} // namespace SFCGAL
