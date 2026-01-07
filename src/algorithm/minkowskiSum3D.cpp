// Copyright (c) 2012-2024, SFCGAL Contributors and Oslandia
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

/**
 * @file minkowskiSum3D.cpp
 * @brief Implementation of 3D Minkowski sum algorithm using CGAL Nef polyhedra.
 *
 * This implementation wraps CGAL::minkowski_sum_3 to compute the Minkowski sum
 * of two 3D geometries. The algorithm supports all SFCGAL geometry types:
 *
 * - **Points**: Converted to singular vertices in Nef polyhedron
 * - **LineStrings**: Converted to polylines using CGAL's native Polylines_tag
 * - **NURBSCurves**: Converted to LineString first, then to polylines
 * - **Triangles**: Converted to singular facets
 * - **Polygons**: Triangulated and converted to facets (handles non-convex and
 * holes)
 * - **TriangulatedSurface/PolyhedralSurface**: Converted directly to Nef
 * polyhedra
 * - **Solids**: Exterior shell converted, interior shells subtracted
 * - **Multi* and GeometryCollection**: Combined via union
 *
 * According to CGAL documentation, input polyhedra may consist of:
 * 1. Singular vertices
 * 2. Singular edges (polylines)
 * 3. Singular convex facets without holes
 * 4. Surfaces with convex facets without holes
 * 5. Three-dimensional features with coplanar facets having common selection
 * marks
 *
 * @see https://doc.cgal.org/latest/Minkowski_sum_3/index.html
 */

#include "SFCGAL/algorithm/minkowskiSum3D.h"
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
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/detail/GeometrySet.h"
#include "SFCGAL/triangulate/triangulatePolygon.h"

#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/minkowski_sum_3.h>

#include <list>
#include <vector>

namespace SFCGAL::algorithm {

// ----------------------------------------------------------------------------------
// -- Type aliases
// ----------------------------------------------------------------------------------
namespace {

/// @brief Nef polyhedron type used throughout this implementation
using Nef_polyhedron_3 = CGAL::Nef_polyhedron_3<Kernel>;

/// @brief Point iterator type for polyline handling
using PointIterator = Kernel::Point_3 *;

/// @brief Point range as iterator pair for CGAL polyline constructor
using PointRange = std::pair<PointIterator, PointIterator>;

/// @brief List of polyline ranges for CGAL Nef_polyhedron_3 Polylines_tag
/// constructor
using PolylineList = std::list<PointRange>;

// ----------------------------------------------------------------------------------
// -- Empty result helpers
// ----------------------------------------------------------------------------------

/**
 * @brief Create an empty Nef polyhedron
 * @return An empty Nef_polyhedron_3 instance
 */
inline auto
emptyNef() -> Nef_polyhedron_3
{
  return {};
}

/**
 * @brief Create an empty geometry result (GeometryCollection)
 * @return A unique_ptr to an empty GeometryCollection
 */
inline auto
emptyResult() -> std::unique_ptr<Geometry>
{
  return std::make_unique<GeometryCollection>();
}

// ----------------------------------------------------------------------------------
// -- Nef union accumulator helper
// ----------------------------------------------------------------------------------

/**
 * @brief Accumulate Nef polyhedra via union operation
 *
 * This template helper reduces duplication between polygonToNef() and
 * geometryCollectionToNef() by providing a generic way to combine
 * multiple Nef polyhedra into a single result via union.
 *
 * @tparam Generator A callable type that takes a size_t index and returns
 *                   a Nef_polyhedron_3
 * @param count Number of elements to process
 * @param gen Generator function that produces Nef_polyhedron_3 for each index
 * @return Combined Nef_polyhedron_3 (union of all non-empty results)
 */
template <typename Generator>
auto
accumulateNefUnion(size_t count, Generator &&gen) -> Nef_polyhedron_3
{
  Nef_polyhedron_3 result;
  for (size_t i = 0; i < count; ++i) {
    Nef_polyhedron_3 nef = gen(i);
    if (nef.is_empty()) {
      continue;
    }
    if (result.is_empty()) {
      result = std::move(nef);
    } else {
      result += nef;
    }
  }
  return result;
}

} // anonymous namespace

// ----------------------------------------------------------------------------------
// -- Private converter functions
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

/**
 * @brief Convert a Point to Nef_polyhedron_3
 *
 * CGAL Nef_polyhedron_3 natively supports singular vertices.
 *
 * @param point The Point to convert
 * @return Nef_polyhedron_3 containing a singular vertex
 */
[[nodiscard]] auto
pointToNef(const Point &point) -> Nef_polyhedron_3
{
  return Nef_polyhedron_3(point.toPoint_3());
}

/**
 * @brief Convert a LineString to Nef_polyhedron_3 using CGAL's native polyline
 * support
 *
 * CGAL Nef_polyhedron_3 natively supports singular edges (polylines) via the
 * Polylines_tag constructor. This is the correct way to represent 1D geometry
 * for Minkowski sums, avoiding the previous incorrect prism approximation.
 *
 * Note: The points vector is created locally and must remain valid until
 * the Nef_polyhedron_3 is constructed. The Nef constructor copies the data.
 *
 * @param linestring The LineString to convert
 * @return Nef_polyhedron_3 containing the polyline, or empty if < 2 points
 */
[[nodiscard]] auto
lineStringToNef(const LineString &linestring) -> Nef_polyhedron_3
{
  if (linestring.numPoints() < 2) {
    return emptyNef();
  }

  // Create a local copy of the points
  // The vector must remain in scope until Nef_polyhedron_3 is constructed
  std::vector<Kernel::Point_3> points;
  points.reserve(linestring.numPoints());

  for (size_t i = 0; i < linestring.numPoints(); ++i) {
    points.push_back(linestring.pointN(i).toPoint_3());
  }

  // Create polyline list with pointer range into local storage
  PolylineList polylines;
  polylines.emplace_back(&points.front(), &points.back() + 1);

  // Use CGAL's native polyline constructor
  // The constructor copies the data, so the local vector can be destroyed after
  return {polylines.begin(), polylines.end(),
          Nef_polyhedron_3::Polylines_tag()};
}

/**
 * @brief Convert a Triangle to Nef_polyhedron_3
 *
 * Triangles are (surface) primitives. CGAL Nef_polyhedron_3 supports
 * singular convex facets without holes.
 *
 * @param triangle The Triangle to convert
 * @return Nef_polyhedron_3 containing the triangle facet, or empty if triangle
 * is empty
 */
[[nodiscard]] auto
triangleToNef(const Triangle &triangle) -> Nef_polyhedron_3
{
  if (triangle.isEmpty()) {
    return emptyNef();
  }

  Polyhedron_3 poly;
  poly.make_triangle(triangle.vertex(0).toPoint_3(),
                     triangle.vertex(1).toPoint_3(),
                     triangle.vertex(2).toPoint_3());
  return Nef_polyhedron_3(poly);
}

/**
 * @brief Convert a Polygon to Nef_polyhedron_3
 *
 * Properly handles both convex and non-convex polygons by triangulating them
 * first. Also properly handles polygons with interior rings (holes).
 *
 * According to CGAL documentation, Minkowski sum input polyhedra may consist
 * of "singular convex facets without holes". For non-convex polygons or
 * polygons with holes, we triangulate them to get a set of convex triangles.
 *
 * @param polygon The Polygon to convert
 * @return Nef_polyhedron_3 containing the triangulated polygon facets, or empty
 * if polygon is empty
 */
[[nodiscard]] auto
polygonToNef(const Polygon &polygon) -> Nef_polyhedron_3
{
  if (polygon.isEmpty()) {
    return emptyNef();
  }

  // Triangulate the polygon to handle non-convex shapes and interior rings
  TriangulatedSurface triangulatedSurface;
  try {
    triangulate::triangulatePolygon3D(polygon, triangulatedSurface);
  } catch (...) {
    // If triangulation fails, try to use convex hull as fallback
    // This maintains backward compatibility for simple convex polygons
    // But, may considred to be removed in the next version
    std::vector<Kernel::Point_3> points;
    points.reserve(polygon.exteriorRing().numPoints() - 1);
    for (size_t i = 0; i < polygon.exteriorRing().numPoints() - 1; ++i) {
      points.push_back(polygon.exteriorRing().pointN(i).toPoint_3());
    }
    Polyhedron_3 cgalPoly;
    CGAL::convex_hull_3(points.begin(), points.end(), cgalPoly);
    if (!cgalPoly.is_empty()) {
      return Nef_polyhedron_3(cgalPoly);
    }
    return emptyNef();
  }

  if (triangulatedSurface.numPatches() == 0) {
    return emptyNef();
  }

  // Build Nef from triangulated surface using accumulator helper
  return accumulateNefUnion(
      triangulatedSurface.numPatches(), [&triangulatedSurface](size_t index) {
        const Triangle &tri = triangulatedSurface.patchN(index);
        Polyhedron_3    triPoly;
        triPoly.make_triangle(tri.vertex(0).toPoint_3(),
                              tri.vertex(1).toPoint_3(),
                              tri.vertex(2).toPoint_3());
        return Nef_polyhedron_3(triPoly);
      });
}

/**
 * @brief Convert a TriangulatedSurface to Nef_polyhedron_3
 *
 * @param surface The TriangulatedSurface to convert
 * @return Nef_polyhedron_3 containing the surface, or empty if surface is empty
 */
[[nodiscard]] auto
triangulatedSurfaceToNef(const TriangulatedSurface &surface) -> Nef_polyhedron_3
{
  if (surface.isEmpty()) {
    return emptyNef();
  }

  std::unique_ptr<Polyhedron_3> polyPtr =
      surface.toPolyhedron_3<Polyhedron_3>();
  if (polyPtr && !polyPtr->is_empty()) {
    return Nef_polyhedron_3(*polyPtr);
  }
  return emptyNef();
}

/**
 * @brief Convert a PolyhedralSurface to Nef_polyhedron_3
 *
 * @param surface The PolyhedralSurface to convert
 * @return Nef_polyhedron_3 containing the surface, or empty if surface is empty
 */
[[nodiscard]] auto
polyhedralSurfaceToNef(const PolyhedralSurface &surface) -> Nef_polyhedron_3
{
  if (surface.isEmpty()) {
    return emptyNef();
  }

  std::unique_ptr<Polyhedron_3> polyPtr =
      surface.toPolyhedron_3<Polyhedron_3>();
  if (polyPtr && !polyPtr->is_empty()) {
    return Nef_polyhedron_3(*polyPtr);
  }
  return emptyNef();
}

/**
 * @brief Convert a Solid to Nef_polyhedron_3
 *
 * Handles solids with interior shells (voids) by computing the exterior shell
 * and subtracting interior shells.
 *
 * @param solid The Solid to convert
 * @return Nef_polyhedron_3 containing the solid volume, or empty if solid is
 * empty
 */
[[nodiscard]] auto
solidToNef(const Solid &solid) -> Nef_polyhedron_3
{
  if (solid.isEmpty()) {
    return emptyNef();
  }

  // Convert exterior shell
  Nef_polyhedron_3 result = polyhedralSurfaceToNef(solid.exteriorShell());

  // Handle interior shells (voids) - subtract them from the exterior
  for (size_t i = 0; i < solid.numInteriorShells(); ++i) {
    Nef_polyhedron_3 interior = polyhedralSurfaceToNef(solid.interiorShellN(i));
    if (!interior.is_empty()) {
      result -= interior;
    }
  }

  return result;
}

// Forward declaration for recursive geometry collection handling
[[nodiscard]] auto
geometryToNef(const Geometry &g) -> Nef_polyhedron_3;

/**
 * @brief Convert a GeometryCollection (or Multi* type) to Nef_polyhedron_3
 *
 * Combines all geometries in the collection via union.
 *
 * @param collection The GeometryCollection to convert
 * @return Nef_polyhedron_3 containing the union of all geometries, or empty if
 * collection is empty
 */
[[nodiscard]] auto
geometryCollectionToNef(const GeometryCollection &collection)
    -> Nef_polyhedron_3
{
  if (collection.isEmpty()) {
    return emptyNef();
  }

  // Use accumulator helper to combine all geometries via union
  return accumulateNefUnion(collection.numGeometries(),
                            [&collection](size_t index) {
                              return geometryToNef(collection.geometryN(index));
                            });
}

/**
 * @brief Convert any SFCGAL::Geometry to Nef_polyhedron_3
 *
 * This is the main conversion function that dispatches to type-specific
 * converters.
 *
 * @param g The geometry to convert
 * @return Nef_polyhedron_3 representation of the geometry
 * @throws GeometryInvalidityException for unsupported geometry types
 */
// NOLINTBEGIN(readability-function-cognitive-complexity)
[[nodiscard]] auto
geometryToNef(const Geometry &g) -> Nef_polyhedron_3
{
  if (g.isEmpty()) {
    return emptyNef();
  }

  switch (g.geometryTypeId()) {
  case TYPE_POINT:
    return pointToNef(g.as<Point>());

  case TYPE_LINESTRING:
    return lineStringToNef(g.as<LineString>());

  case TYPE_NURBSCURVE: {
    // Convert NURBS curve to LineString using default parameters, then to Nef
    // This matches the pattern used in minkowskiSum (2D)
    const auto &nurbs      = g.as<NURBSCurve>();
    auto        lineString = nurbs.toLineString(); // default discretization
    if (!lineString || lineString->isEmpty()) {
      return emptyNef();
    }
    return lineStringToNef(*lineString);
  }

  case TYPE_TRIANGLE:
    return triangleToNef(g.as<Triangle>());

  case TYPE_POLYGON:
    return polygonToNef(g.as<Polygon>());

  case TYPE_TRIANGULATEDSURFACE:
    return triangulatedSurfaceToNef(g.as<TriangulatedSurface>());

  case TYPE_POLYHEDRALSURFACE:
    return polyhedralSurfaceToNef(g.as<PolyhedralSurface>());

  case TYPE_SOLID:
    return solidToNef(g.as<Solid>());

  case TYPE_MULTIPOINT:
    return geometryCollectionToNef(g.as<MultiPoint>());

  case TYPE_MULTILINESTRING:
    return geometryCollectionToNef(g.as<MultiLineString>());

  case TYPE_MULTIPOLYGON:
    return geometryCollectionToNef(g.as<MultiPolygon>());

  case TYPE_MULTISOLID:
    return geometryCollectionToNef(g.as<MultiSolid>());

  case TYPE_GEOMETRYCOLLECTION:
    return geometryCollectionToNef(g.as<GeometryCollection>());

  default:
    BOOST_THROW_EXCEPTION(GeometryInvalidityException(
        "Unsupported geometry type for Minkowski sum 3D: " + g.geometryType() +
        " (type id: " + std::to_string(static_cast<int>(g.geometryTypeId())) +
        ")"));
  }
}
// NOLINTEND(readability-function-cognitive-complexity)

/**
 * @brief Convert Nef_polyhedron_3 to SFCGAL::Geometry
 *
 * @param nef The Nef polyhedron to convert
 * @return unique_ptr<Geometry> containing the result (PolyhedralSurface or
 * empty GeometryCollection)
 */
[[nodiscard]] auto
nefToGeometry(const Nef_polyhedron_3 &nef) -> std::unique_ptr<Geometry>
{
  if (nef.is_empty()) {
    return emptyResult();
  }

  Polyhedron_3 poly;
  nef.convert_to_polyhedron(poly);

  if (poly.is_empty()) {
    return emptyResult();
  }

  return std::make_unique<PolyhedralSurface>(poly);
}

auto
minkowskiSum3D(const Geometry &gA, const Geometry &gB,
               NoValidityCheck /*unused*/) -> std::unique_ptr<Geometry>
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return emptyResult();
  }

  Nef_polyhedron_3 nefA = geometryToNef(gA);
  Nef_polyhedron_3 nefB = geometryToNef(gB);

  if (nefA.is_empty() || nefB.is_empty()) {
    return emptyResult();
  }

  Nef_polyhedron_3 result = CGAL::minkowski_sum_3(nefA, nefB);

  if (result.is_empty()) {
    return emptyResult();
  }

  return nefToGeometry(result);
}

auto
minkowskiSum3D(const Geometry &gA, const Geometry &gB)
    -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(gA);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(gB);

  std::unique_ptr<Geometry> result(minkowskiSum3D(gA, gB, NoValidityCheck()));
  propagateValidityFlag(*result, true);
  return result;
}

/// @} end of private section
} // namespace SFCGAL::algorithm
