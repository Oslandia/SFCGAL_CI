// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_DETAIL_ENVELOPEVISITOR_H_
#define SFCGAL_DETAIL_ENVELOPEVISITOR_H_

#include "SFCGAL/Envelope.h"
#include "SFCGAL/GeometryVisitor.h"
#include "SFCGAL/config.h"

#include <vector>

/**
 * Visitor that extracts point coordinates from geometries to expand an
 * Envelope.
 *
 * Visits concrete Geometry types and updates the referenced `envelope` with all
 * point coordinates encountered (vertices, control points, etc.), so the
 * envelope bounds encompass the geometry. Supports points, line strings,
 * polygons, triangles, solids, multis, collections, surfaces, triangulated
 * surfaces and NURBS curves.
 *
 * The visitor holds a reference to the Envelope to populate (public member
 * `envelope`). Use an instance of this visitor to traverse a Geometry with
 * the goal of computing or extending its bounding envelope.
 */
namespace SFCGAL::detail {

/**
 * Get the list of points from a Geometry
 *
 * @todo ConstPointVisitor
 */
class SFCGAL_API EnvelopeVisitor : public ConstGeometryVisitor {
public:
  /**
   * @brief Constructor for envelope visitor
   * @param envelope_ Reference to envelope to update
   */
  EnvelopeVisitor(Envelope &envelope_);

  /**
   * @brief Visit a Point and update envelope
   * @param g The point to visit
   */
  void
  visit(const Point &g) override;
  /**
   * @brief Visit a LineString and update envelope
   * @param g The linestring to visit
   */
  void
  visit(const LineString &g) override;
  /**
   * @brief Visit a Polygon and update envelope
   * @param g The polygon to visit
   */
  void
  visit(const Polygon &g) override;
  /**
   * @brief Visit a Triangle and update envelope
   * @param g The triangle to visit
   */
  void
  visit(const Triangle &g) override;
  /**
   * @brief Visit a Solid and update envelope
   * @param g The solid to visit
   */
  void
  visit(const Solid &g) override;
  /**
   * @brief Visit a MultiPoint and update envelope
   * @param g The multipoint to visit
   */
  void
  visit(const MultiPoint &g) override;
  /**
   * @brief Visit a MultiLineString and update envelope
   * @param g The multilinestring to visit
   */
  void
  visit(const MultiLineString &g) override;
  /**
   * @brief Visit a MultiPolygon and update envelope
   * @param g The multipolygon to visit
   */
  void
  visit(const MultiPolygon &g) override;
  /**
   * @brief Visit a MultiSolid and update envelope
   * @param g The multisolid to visit
   */
  void
  visit(const MultiSolid &g) override;
  /**
   * @brief Visit a GeometryCollection and update envelope
   * @param g The geometry collection to visit
   */
  void
  visit(const GeometryCollection &g) override;
  /**
   * @brief Visit a PolyhedralSurface and update envelope
   * @param g The polyhedral surface to visit
   */
  void
  visit(const PolyhedralSurface &g) override;
  /**
   * @brief Visit a TriangulatedSurface and update envelope
   * @param g The triangulated surface to visit
   */
  void
  visit(const TriangulatedSurface &g) override;
  /// @brief Process NURBSCurve to update envelope bounds
  /// @param g The NURBSCurve geometry to process
  void
  visit(const NURBSCurve &g) override;

public:
  Envelope &envelope; ///< Reference to the envelope being computed
};

} // namespace SFCGAL::detail

#endif
