// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/isClosed.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/algorithm/distance.h"

#include <limits>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/helpers.h>

namespace SFCGAL::algorithm {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------

namespace detail {

auto
isClosedLineString(const LineString &linestring) -> Closure
{
  if (linestring.isEmpty()) {
    return Closure::closed();
  }

  if (linestring.numPoints() < 2) {
    return Closure::open("LineString has less than 2 points");
  }

  // Check if first and last points are identical
  if (linestring.startPoint() != linestring.endPoint()) {
    return Closure::open("First and last points are not identical");
  }

  return Closure::closed();
}

auto
isClosedPolyhedralSurface(const PolyhedralSurface &surface) -> Closure
{
  if (surface.isEmpty()) {
    return Closure::closed();
  }

  try {
    // Convert to CGAL Polyhedron to check closure
    auto cgal_poly =
        surface
            .toPolyhedron_3<SFCGAL::Kernel, SFCGAL::detail::MarkedPolyhedron>();

    // Check if the polyhedron is closed (no boundary edges)
    if (!cgal_poly->is_closed()) {
      // Count boundary edges for diagnostic
      int boundary_count = 0;
      for (auto it = cgal_poly->halfedges_begin();
           it != cgal_poly->halfedges_end(); ++it) {
        if (it->is_border()) {
          boundary_count++;
        }
      }

      return Closure::open(
          (boost::format("Surface has %d boundary halfedges") % boundary_count)
              .str());
    }

    return Closure::closed();

  } catch (const std::exception &e) {
    return Closure::open(
        (boost::format("Failed to convert to polyhedron: %s") % e.what())
            .str());
  }
}

auto
isClosedTriangulatedSurface(const TriangulatedSurface &surface) -> Closure
{
  if (surface.isEmpty()) {
    return Closure::closed();
  }

  try {
    // Convert to CGAL Polyhedron to check closure
    auto cgal_poly =
        surface
            .toPolyhedron_3<SFCGAL::Kernel, SFCGAL::detail::MarkedPolyhedron>();

    // Check if the polyhedron is closed
    if (!cgal_poly->is_closed()) {
      // Count boundary edges for diagnostic
      int boundary_count = 0;
      for (auto it = cgal_poly->halfedges_begin();
           it != cgal_poly->halfedges_end(); ++it) {
        if (it->is_border()) {
          boundary_count++;
        }
      }

      return Closure::open(
          (boost::format("Surface has %d boundary halfedges") % boundary_count)
              .str());
    }

    return Closure::closed();

  } catch (const std::exception &e) {
    return Closure::open(
        (boost::format("Failed to convert to polyhedron: %s") % e.what())
            .str());
  }
}

auto
isClosedSolid(const Solid &solid) -> Closure
{
  // Solids are closed by definition
  // However, we could check if all shells are closed
  const size_t numShells = solid.numShells();

  for (size_t i = 0; i < numShells; ++i) {
    Closure shellClosure = isClosedPolyhedralSurface(solid.shellN(i));
    if (!shellClosure) {
      return Closure::open((boost::format("Shell %d is not closed: %s") % i %
                            shellClosure.reason())
                               .str());
    }
  }

  return Closure::closed();
}

auto
isClosedMultiLineString(const MultiLineString &mls) -> Closure
{
  if (mls.isEmpty()) {
    return Closure::closed();
  }

  const size_t numGeom = mls.numGeometries();

  for (size_t i = 0; i < numGeom; ++i) {
    Closure ls_closure = isClosedLineString(mls.lineStringN(i));
    if (!ls_closure) {
      return Closure::open((boost::format("LineString %d is not closed: %s") %
                            i % ls_closure.reason())
                               .str());
    }
  }

  return Closure::closed();
}

auto
isClosedMultiSolid(const MultiSolid &msolid) -> Closure
{
  // MultiSolids are always closed (all solids are closed)
  // But we could verify each solid
  const size_t numGeom = msolid.numGeometries();

  for (size_t i = 0; i < numGeom; ++i) {
    Closure solid_closure = isClosedSolid(msolid.solidN(i));
    if (!solid_closure) {
      return Closure::open((boost::format("Solid %d is not closed: %s") % i %
                            solid_closure.reason())
                               .str());
    }
  }

  return Closure::closed();
}

auto
isClosedNURBSCurve(const NURBSCurve &nurbsCurve) -> Closure
{
  if (nurbsCurve.isEmpty()) {
    return Closure::closed();
  }

  if (nurbsCurve.numControlPoints() < 2) {
    return Closure::open("NURBSCurve has less than 2 points");
  }

  // Use efficient endpoint accessors and squared distance comparison
  Point startPoint = nurbsCurve.startPoint();
  Point endPoint   = nurbsCurve.endPoint();

  // Compute squared distance to avoid sqrt
  auto dx = endPoint.x() - startPoint.x();
  auto dy = endPoint.y() - startPoint.y();
  auto dz = endPoint.is3D() && startPoint.is3D() ? endPoint.z() - startPoint.z()
                                                 : NURBSCurve::FT(0);
  auto squaredDistance = dx * dx + dy * dy + dz * dz;

  // Use relative tolerance based on curve length or numeric epsilon
  auto curveLength = nurbsCurve.length();
  auto tolerance =
      curveLength > NURBSCurve::FT(0)
          ? curveLength * std::numeric_limits<double>::epsilon() * 100
          : NURBSCurve::FT(std::numeric_limits<double>::epsilon());
  auto squaredTolerance = tolerance * tolerance;

  if (squaredDistance <= squaredTolerance) {
    return Closure::closed();
  }

  // Return actual distance (not squared) in error message
  auto actualDistance = CGAL::sqrt(CGAL::to_double(squaredDistance));
  return Closure::open(
      (boost::format(
           "NURBS curve endpoints are not coincident (distance: %g)") %
       actualDistance)
          .str());
}

} // namespace detail

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------

auto
isClosed(const Geometry &g) -> Closure
{
  switch (g.geometryTypeId()) {
  case TYPE_POINT:
  case TYPE_POLYGON:
  case TYPE_TRIANGLE:
  case TYPE_MULTIPOINT:
  case TYPE_MULTIPOLYGON:
    return Closure::closed();

  case TYPE_LINESTRING:
    return detail::isClosedLineString(g.as<LineString>());

  case TYPE_POLYHEDRALSURFACE:
    return detail::isClosedPolyhedralSurface(g.as<PolyhedralSurface>());

  case TYPE_TRIANGULATEDSURFACE:
    return detail::isClosedTriangulatedSurface(g.as<TriangulatedSurface>());

  case TYPE_SOLID:
    return detail::isClosedSolid(g.as<Solid>());

  case TYPE_MULTILINESTRING:
    return detail::isClosedMultiLineString(g.as<MultiLineString>());

  case TYPE_MULTISOLID:
    return detail::isClosedMultiSolid(g.as<MultiSolid>());

  case TYPE_NURBSCURVE:
    return detail::isClosedNURBSCurve(g.as<NURBSCurve>());

  case TYPE_GEOMETRYCOLLECTION: {
    const auto &collection = g.as<GeometryCollection>();

    if (collection.isEmpty()) {
      return Closure::closed();
    }

    const size_t numGeom = collection.numGeometries();

    for (size_t i = 0; i < numGeom; ++i) {
      Closure geom_closure = isClosed(collection.geometryN(i));
      if (!geom_closure) {
        return Closure::open(
            (boost::format("Geometry %d (%s) is not closed: %s") % i %
             collection.geometryN(i).geometryType() % geom_closure.reason())
                .str());
      }
    }

    return Closure::closed();
  }
  }

  BOOST_THROW_EXCEPTION(Exception(
      (boost::format("isClosed( %s ) is not defined") % g.geometryType())
          .str()));

  return Closure::open("Unknown geometry type");
}

} // namespace SFCGAL::algorithm
