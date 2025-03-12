// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/triangulate/triangulate2DZ.h"

#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/Triangle.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/algorithm/isValid.h"

namespace SFCGAL::triangulate {

using Vertex_handle = ConstraintDelaunayTriangulation::Vertex_handle;

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

///
///
///
void
triangulate2DZ(const Point &g, ConstraintDelaunayTriangulation &triangulation)
{
  triangulation.addVertex(g.coordinate());
}
///
///
///
void
triangulate2DZ(const LineString                &g,
               ConstraintDelaunayTriangulation &triangulation)
{
  Vertex_handle last;

  for (size_t i = 0; i < g.numPoints(); i++) {
    Vertex_handle const vertex =
        triangulation.addVertex(g.pointN(i).coordinate());

    if (i != 0) {
      triangulation.addConstraint(last, vertex);
    }

    last = vertex;
  }
}
///
///
///
void
triangulate2DZ(const Polygon &g, ConstraintDelaunayTriangulation &triangulation)
{
  for (size_t i = 0; i < g.numRings(); i++) {
    triangulate2DZ(g.ringN(i), triangulation);
  }
}
///
///
///
void
triangulate2DZ(const Triangle                  &g,
               ConstraintDelaunayTriangulation &triangulation)
{
  Vertex_handle last;

  for (size_t i = 0; i < 4; i++) {
    Vertex_handle const vertex =
        triangulation.addVertex(g.vertex(i).coordinate());

    if (i != 0) {
      triangulation.addConstraint(last, vertex);
    }

    last = vertex;
  }
}

/// @} end of private section

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

///
///
///
void
triangulateCollection2DZ(const Geometry                  &g,
                         ConstraintDelaunayTriangulation &triangulation)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    triangulate2DZ(g.geometryN(i), triangulation);
  }
}

///
///
///
void
triangulate2DZ(const Geometry                  &g,
               ConstraintDelaunayTriangulation &triangulation)
{

  if (g.isEmpty()) {
    return;
  }

  if (triangulation.hasProjectionPlane()) {
    SFCGAL_ASSERT_GEOMETRY_VALIDITY_ON_PLANE(g);
  } else {
    SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(g);
  }

  switch (g.geometryTypeId()) {
  case TYPE_POINT:
    triangulate2DZ(g.as<Point>(), triangulation);
    return;

  case TYPE_LINESTRING:
    triangulate2DZ(g.as<LineString>(), triangulation);
    return;

  case TYPE_POLYGON:
    triangulate2DZ(g.as<Polygon>(), triangulation);
    return;

  case TYPE_TRIANGLE:
    triangulate2DZ(g.as<Triangle>(), triangulation);
    return;

  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_POLYHEDRALSURFACE:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_GEOMETRYCOLLECTION:
    triangulateCollection2DZ(g, triangulation);
    return;

  case TYPE_SOLID:
  case TYPE_MULTISOLID:
    // note: we can't have a valid geom in 2D that comes from a solid
    // since a solid closed and thus self-intersects once projected
    BOOST_THROW_EXCEPTION(InappropriateGeometryException(
        (boost::format("can't process 2DZ triangulation for type '%1%'") %
         g.geometryType())
            .str()));
  }
}

///
///
///
auto
triangulate2DZ(const Geometry &g) -> ConstraintDelaunayTriangulation
{
  ConstraintDelaunayTriangulation triangulation;
  triangulate2DZ(g, triangulation);
  return triangulation;
}

} // namespace SFCGAL::triangulate
