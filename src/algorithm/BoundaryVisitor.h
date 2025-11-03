// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_BOUNDARYVISITOR_H_
#define SFCGAL_ALGORITHM_BOUNDARYVISITOR_H_

#include "SFCGAL/GeometryVisitor.h"

#include "SFCGAL/detail/graph/GeometryGraph.h"
#include "SFCGAL/detail/graph/GeometryGraphBuilder.h"

/**
 * Compute the topological boundary of a Geometry.
 *
 * This visitor computes the boundary for various geometry types:
 * - Point, MultiPoint: empty geometry collection.
 * - LineString, Triangle: empty if closed, otherwise a MultiPoint of endpoints.
 * - Polygon: the polygon's exterior and interior rings as a LineString or
 *   MultiLineString (one per ring).
 * - MultiLineString: either empty (if all component lines are closed) or the
 *   collection of unmatched endpoints.
 * - MultiPolygon, PolyhedralSurface, TriangulatedSurface, Solid, MultiSolid:
 *   handled where supported; see TODOs for incomplete cases.
 * - NURBSCurve: processed for boundary calculation like other curve types.
 *
 * GeometryCollection handling is limited and not supported in the general,
 * heterogeneous case.
 *
 * The computed boundary is stored internally and returned via
 * releaseBoundary().
 *
 * @warning GeometryCollection is not supported in the general case.
 * @todo Improve support for Solid, MultiPolygon, PolyhedralSurface,
 *       TriangulatedSurface and heterogeneous GeometryCollection.
 */
namespace SFCGAL::algorithm {

/**
 * Compute the boundary for a Geometry
 *
 * boundary( Point )      : GEOMETRYCOLLECTION EMPTY
 * boundary( LineString ) : either GEOMETRYCOLLECTION EMPTY is the LineString is
 * closed, or MULTIPOINT(2) boundary( Polygon )    : LINESTRING |
 * MULTILINESTRING (polygon rings) boundary( Triangle )   : either
 * GEOMETRYCOLLECTION EMPTY is the LineString is closed, or MULTIPOINT(2)
 *
 * boundary( MultiPoint )      : GEOMETRYCOLLECTION EMPTY
 * boundary( MultiLineString ) : either GEOMETRYCOLLECTION EMPTY or single
 * occurance points
 *
 * @warning GeometryCollection are not supported in the general case
 *
 * @ŧodo Solid
 *
 * @todo MultiPolygon, PolyhedralSurface, TriangulatedSurface (same graph
 * algorithm, edges without parallel or opposite)
 *
 *
 * @todo GeometryCollection : complex for heterogeneous collection (not
 * supported in GEOS)
 * @todo MultiSolid : faced elimination
 *
 */
class SFCGAL_API BoundaryVisitor : public ConstGeometryVisitor {
public:
  /// @brief Visit a Point geometry to compute its boundary
  /// @param g The Point geometry to visit
  void
  visit(const Point &g) override;
  /// @brief Visit a LineString geometry to compute its boundary
  /// @param g The LineString geometry to visit
  void
  visit(const LineString &g) override;
  /// @brief Visit a Polygon geometry to compute its boundary
  /// @param g The Polygon geometry to visit
  void
  visit(const Polygon &g) override;
  /// @brief Visit a Triangle geometry to compute its boundary
  /// @param g The Triangle geometry to visit
  void
  visit(const Triangle &g) override;
  /// @brief Visit a Solid geometry to compute its boundary
  /// @param g The Solid geometry to visit
  void
  visit(const Solid &g) override;
  /// @brief Visit a MultiPoint geometry to compute its boundary
  /// @param g The MultiPoint geometry to visit
  void
  visit(const MultiPoint &g) override;
  /// @brief Visit a MultiLineString geometry to compute its boundary
  /// @param g The MultiLineString geometry to visit
  void
  visit(const MultiLineString &g) override;
  /// @brief Visit a MultiPolygon geometry to compute its boundary
  /// @param g The MultiPolygon geometry to visit
  void
  visit(const MultiPolygon &g) override;
  /// @brief Visit a MultiSolid geometry to compute its boundary
  /// @param g The MultiSolid geometry to visit
  void
  visit(const MultiSolid &g) override;
  /// @brief Visit a GeometryCollection to compute its boundary
  /// @param g The GeometryCollection to visit
  void
  visit(const GeometryCollection &g) override;
  /// @brief Visit a PolyhedralSurface geometry to compute its boundary
  /// @param g The PolyhedralSurface geometry to visit
  void
  visit(const PolyhedralSurface &g) override;
  /// @brief Visit a TriangulatedSurface geometry to compute its boundary
  /// @param g The TriangulatedSurface geometry to visit
  void
  visit(const TriangulatedSurface &g) override;
  /// @brief Process NURBSCurve for boundary calculation
  /// @param g The NURBSCurve geometry to process
  void
  visit(const NURBSCurve &g) override;

  /**
   * get the boundary
   */
  Geometry *
  releaseBoundary();

protected:
  /**
   * get the boundary vertices for a set of LineString in a GeometryGraph
   * @param g The GeometryGraph containing LineStrings
   */
  void
  getBoundaryFromLineStrings(const graph::GeometryGraph &g);
  /**
   * get the boundary edges for a set of Polygons in a GeometryGraph
   * @param g The GeometryGraph containing Polygons
   * @warning not optimal (edges could be counted using complex<
   * vertex_descriptor >)
   * @todo merge resulting edges
   */
  void
  getBoundaryFromPolygons(const graph::GeometryGraph &g);

private:
  std::unique_ptr<Geometry> _boundary;
};

} // namespace SFCGAL::algorithm

#endif
