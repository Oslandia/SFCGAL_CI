// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/triangulate/triangulatePolygon.h>

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/Kernel.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Solid.h>

#include <SFCGAL/detail/triangulate/ConstraintDelaunayTriangulation.h>

#include <SFCGAL/algorithm/isValid.h>
#include <SFCGAL/algorithm/normal.h>
#include <SFCGAL/algorithm/plane.h>

#include <iostream>

namespace SFCGAL {
namespace triangulate {

using Vertex_handle = ConstraintDelaunayTriangulation::Vertex_handle;

///
///
///
void
triangulatePolygon3D(const Geometry      &g,
                     TriangulatedSurface &triangulatedSurface)
{
  if (g.isEmpty()) {
    return;
  }

  SFCGAL_ASSERT_GEOMETRY_VALIDITY(g);

  switch (g.geometryTypeId()) {
  case TYPE_TRIANGLE:
    return triangulatePolygon3D(g.as<Triangle>(), triangulatedSurface);

  case TYPE_POLYGON:
    return triangulatePolygon3D(g.as<Polygon>(), triangulatedSurface);

  case TYPE_TRIANGULATEDSURFACE:
    return triangulatePolygon3D(g.as<TriangulatedSurface>(),
                                triangulatedSurface);

  case TYPE_POLYHEDRALSURFACE:
    return triangulatePolygon3D(g.as<PolyhedralSurface>(), triangulatedSurface);

  case TYPE_SOLID:
    return triangulatePolygon3D(g.as<Solid>(), triangulatedSurface);

  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION: {
    for (size_t i = 0; i < g.numGeometries(); i++) {
      triangulatePolygon3D(g.geometryN(i), triangulatedSurface);
    }

    return;
  }

  default:
    BOOST_THROW_EXCEPTION(InappropriateGeometryException(
        (boost::format("can't triangulate 3d polygons for type '%1%'") %
         g.geometryType())
            .str()));
  }
}

///
///
///
void
triangulatePolygon3D(const Triangle      &triangle,
                     TriangulatedSurface &triangulatedSurface)
{
  /*
   * filter empty polygon
   */
  if (triangle.isEmpty()) {
    return;
  }

  triangulatedSurface.addTriangle(triangle);
}

///
///
///
void
triangulatePolygon3D(const TriangulatedSurface &g,
                     TriangulatedSurface       &triangulatedSurface)
{
  triangulatedSurface.addTriangles(g);
}

///
///
///
void
triangulatePolygon3D(const GeometryCollection &g,
                     TriangulatedSurface      &triangulatedSurface)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    triangulatePolygon3D(g.geometryN(i), triangulatedSurface);
  }
}

///
///
///
void
triangulatePolygon3D(const Polygon       &polygon,
                     TriangulatedSurface &triangulatedSurface)
{
  /*
   * filter empty polygon
   */
  if (polygon.isEmpty()) {
    return;
  }

  /*
   * Make sure polygon is not degenerate
   */
  if (!algorithm::hasPlane3D<Kernel>(polygon)) {
    BOOST_THROW_EXCEPTION(Exception(
        (boost::format("can't find plane for polygon %s") % polygon.asText())
            .str()));
  }

  /*
   * Prepare a Constraint Delaunay Triangulation
   */
  ConstraintDelaunayTriangulation cdt;

  /*
   * find polygon plane
   */
  Kernel::Plane_3 polygonPlane =
      algorithm::plane3D<Kernel>(polygon, algorithm::Plane3DInexactUnsafe());

  cdt.setProjectionPlane(polygonPlane);

  /*
   * insert each ring in the triangulation
   */
  for (size_t i = 0; i < polygon.numRings(); i++) {
    const LineString &ring = polygon.ringN(i);

    /*
     * note, we do not include the last point, since it's equal to the last and
     * that
     */
    if (!ring.numPoints()) {
      continue;
    }

    Vertex_handle v_prev = cdt.addVertex(ring.pointN(0).coordinate());
    Vertex_handle v_0    = v_prev;

    for (size_t j = 1; j < ring.numPoints() - 1; j++) {
      Vertex_handle vh = cdt.addVertex(ring.pointN(j).coordinate());
      cdt.addConstraint(v_prev, vh);
      v_prev = vh;
    }

    cdt.addConstraint(v_prev, v_0);
  }

  /*
   * Mark facets that are inside the domain bounded by the polygon
   */
  cdt.markDomains();
  cdt.getTriangles(triangulatedSurface, true);
}

///
///
///
void
triangulatePolygon3D(const PolyhedralSurface &g,
                     TriangulatedSurface     &triangulatedSurface)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    triangulatePolygon3D(g.polygonN(i), triangulatedSurface);
  }
}

///
///
///
void
triangulatePolygon3D(const Solid &g, TriangulatedSurface &triangulatedSurface)
{
  for (size_t i = 0; i < g.numShells(); i++) {
    triangulatePolygon3D(g.shellN(i), triangulatedSurface);
  }
}

} // namespace triangulate
} // namespace SFCGAL
