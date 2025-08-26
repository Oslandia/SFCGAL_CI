// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/isSimple.h"

#include "SFCGAL/Curve.h"
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

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

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
isSimple(const Polygon &polygon, const double &toleranceAbs = 1e-9)
    -> const Simplicity
{
  // Polygon must be planar (all points in the same plane)
  if (polygon.is3D() && !isPlane3D<Kernel>(polygon, toleranceAbs)) {
    return Simplicity::complex("Points don't lie in the same plane.");
  }
  return Simplicity::simple();
}

auto
isSimple(const PolyhedralSurface &phs, const double &toleranceAbs)
    -> const Simplicity
{
  if (phs.isEmpty()) {
    return Simplicity::simple();
  }
  const size_t numPatches = phs.numPatches();

  for (size_t g = 0; g != numPatches; ++g) {
    Simplicity const s = isSimple(phs.patchN(g), toleranceAbs);

    if (!s) {
      return Simplicity::complex(
          (boost::format("Polygon %d is complex: %s") % g % s.reason()).str());
    }
  }

  return Simplicity::simple();
}

auto
isSimple(const Solid &solid, const double &toleranceAbs) -> const Simplicity
{
  if (solid.isEmpty()) {
    return Simplicity::simple();
  }
  const size_t numGeom = solid.numGeometries();

  for (size_t g = 0; g != numGeom; ++g) {
    Simplicity const s = isSimple(solid.shellN(g), toleranceAbs);

    if (!s) {
      return Simplicity::complex(
          (boost::format("Polygon %d is complex: %s") % g % s.reason()).str());
    }
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
isSimple(const MultiLineString &mls) -> const Simplicity
{
  if (mls.isEmpty()) {
    return Simplicity::simple();
  }
  const size_t numGeom = mls.numGeometries();

  for (size_t g = 0; g != numGeom; ++g) {
    Simplicity const s = isSimple(mls.lineStringN(g));

    if (!s) {
      return Simplicity::complex(
          (boost::format("LineString %d is complex: %s") % g % s.reason())
              .str());
    }
  }

  return Simplicity::simple();
}

auto
isSimple(const MultiPolygon &mpoly, const double &toleranceAbs)
    -> const Simplicity
{
  if (mpoly.isEmpty()) {
    return Simplicity::simple();
  }
  const size_t numGeom = mpoly.numGeometries();
  for (size_t g = 0; g != numGeom; ++g) {
    Simplicity const s = isSimple(mpoly.polygonN(g), toleranceAbs);
    if (!s) {
      return Simplicity::complex(
          (boost::format("Polygon %d is complex: %s") % g % s.reason()).str());
    }
  }

  return Simplicity::simple();
}

auto
isSimple(const MultiSolid &msolid, const double &toleranceAbs)
    -> const Simplicity
{
  if (msolid.isEmpty()) {
    return Simplicity::simple();
  }
  const size_t numGeom = msolid.numGeometries();

  for (size_t g = 0; g != numGeom; ++g) {
    Simplicity const s = isSimple(msolid.solidN(g), toleranceAbs);

    if (!s) {
      return Simplicity::complex(
          (boost::format("Solid %d is complex: %s") % g % s.reason()).str());
    }
  }

  return Simplicity::simple();
}

auto
isSimple(const GeometryCollection &collection, const double &toleranceAbs)
    -> const Simplicity
{
  if (collection.isEmpty()) {
    return Simplicity::simple();
  }
  const size_t numGeom = collection.numGeometries();

  for (size_t g = 0; g != numGeom; ++g) {
    Simplicity const s = isSimple(collection.geometryN(g), toleranceAbs);

    if (!s) {
      return Simplicity::complex(
          (boost::format("%s at index %d is complex: %s") %
           collection.geometryN(g).geometryType() % g % s.reason())
              .str());
    }
  }

  return Simplicity::simple();
}

/// @} end of private section

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

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

  case TYPE_POLYHEDRALSURFACE: // every polygon is simple
    return isSimple(g.as<PolyhedralSurface>(), toleranceAbs);

  case TYPE_TRIANGLE:
    return Simplicity::simple();

  case TYPE_TRIANGULATEDSURFACE: // every triangle is simple
    return Simplicity::simple();

  case TYPE_SOLID: // every phs is simple
    return isSimple(g.as<Solid>(), toleranceAbs);

  case TYPE_MULTIPOINT: // no equal points
    return isSimple(g.as<MultiPoint>());

  case TYPE_MULTILINESTRING: // every ls is simple, and only intersections are
                             // at boundaries
    return isSimple(g.as<MultiLineString>());

  case TYPE_MULTIPOLYGON: // every polygon is simple
    return isSimple(g.as<MultiPolygon>(), toleranceAbs);

  case TYPE_MULTISOLID: // every solid is simple
    return isSimple(g.as<MultiSolid>(), toleranceAbs);

  case TYPE_NURBSCURVE: {
    // Tessellate NURBS and check LineString simpleness to detect
    // self-intersections
    auto lineString = g.as<Curve>().toLineStringAdaptive();
    if (!lineString || lineString->isEmpty()) {
      lineString = g.as<Curve>().toLineString(256);
    }
    if (!lineString || lineString->isEmpty()) {
      return Simplicity::simple(); // Empty curves are simple
    }
    return isSimple(*lineString, toleranceAbs);
  }

  case TYPE_GEOMETRYCOLLECTION: // every geometry is simple
    return isSimple(g.as<GeometryCollection>(), toleranceAbs);
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
