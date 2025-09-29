// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_GEOMETRYVISITOR_H_
#define SFCGAL_GEOMETRYVISITOR_H_

#include "SFCGAL/Geometry.h"

namespace SFCGAL {

class NURBSCurve;

/**
 * Visitor interface for non-const Geometry types.
 *
 * Implementations should provide behavior for each concrete geometry
 * type via the type-specific visit overloads. The destructor is virtual
 * to allow safe subclass destruction.
 */
class SFCGAL_API GeometryVisitor {
public:
  virtual ~GeometryVisitor();

  /**
   * Apply this visitor to a generic Geometry reference.
   *
   * Concrete visitors can override this to provide a common fallback for
   * any Geometry not handled by the type-specific overloads.
   *
   * @param g The geometry to visit.
   */
  virtual void
  visit(Geometry &g);

  /// Visit a Point geometry
  virtual void
  visit(Point &g) = 0;
  /// Visit a LineString geometry
  virtual void
  visit(LineString &g) = 0;
  /// Visit a Polygon geometry
  virtual void
  visit(Polygon &g) = 0;
  /// Visit a Triangle geometry
  virtual void
  visit(Triangle &g) = 0;
  /// Visit a Solid geometry
  virtual void
  visit(Solid &g) = 0;
  /// Visit a MultiPoint geometry
  virtual void
  visit(MultiPoint &g) = 0;
  /// Visit a MultiLineString geometry
  virtual void
  visit(MultiLineString &g) = 0;
  /// Visit a MultiPolygon geometry
  virtual void
  visit(MultiPolygon &g) = 0;
  /// Visit a MultiSolid geometry
  virtual void
  visit(MultiSolid &g) = 0;
  /// Visit a GeometryCollection geometry
  virtual void
  visit(GeometryCollection &g) = 0;
  /// Visit a PolyhedralSurface geometry
  virtual void
  visit(PolyhedralSurface &g) = 0;
  /// Visit a TriangulatedSurface geometry
  virtual void
  visit(TriangulatedSurface &g) = 0;
  /**
   * Visit a NURBSCurve geometry.
   *
   * Default implementation provides no-op; override to handle NURBSCurve.
   *
   * @param g The NURBSCurve to visit.
   */
  virtual void
  visit(NURBSCurve &g);
};

/**
 * Visitor interface for const Geometry types.
 *
 * Implementations should provide behavior for each concrete const
 * geometry type via the type-specific visit overloads. The destructor
 * is virtual to allow safe subclass destruction.
 */
class SFCGAL_API ConstGeometryVisitor {
public:
  virtual ~ConstGeometryVisitor();

  /**
   * Apply this visitor to a generic const Geometry reference.
   *
   * Concrete visitors can override this to provide a common fallback for
   * any const Geometry not handled by the type-specific overloads.
   *
   * @param g The const geometry to visit.
   */
  virtual void
  visit(const Geometry &g);

  /// Visit a Point geometry
  virtual void
  visit(const Point &g) = 0;
  /// Visit a LineString geometry
  virtual void
  visit(const LineString &g) = 0;
  /// Visit a Polygon geometry
  virtual void
  visit(const Polygon &g) = 0;
  /// Visit a Triangle geometry
  virtual void
  visit(const Triangle &g) = 0;
  /// Visit a Solid geometry
  virtual void
  visit(const Solid &g) = 0;
  /// Visit a MultiPoint geometry
  virtual void
  visit(const MultiPoint &g) = 0;
  /// Visit a MultiLineString geometry
  virtual void
  visit(const MultiLineString &g) = 0;
  /// Visit a MultiPolygon geometry
  virtual void
  visit(const MultiPolygon &g) = 0;
  /// Visit a MultiSolid geometry
  virtual void
  visit(const MultiSolid &g) = 0;
  /// Visit a GeometryCollection geometry
  virtual void
  visit(const GeometryCollection &g) = 0;
  /// Visit a PolyhedralSurface geometry
  virtual void
  visit(const PolyhedralSurface &g) = 0;
  /// Visit a TriangulatedSurface geometry
  virtual void
  visit(const TriangulatedSurface &g) = 0;
  /**
   * Visit a const NURBSCurve geometry.
   *
   * Default implementation provides no-op; override to handle const NURBSCurve.
   *
   * @param g The const NURBSCurve to visit.
   */
  virtual void
  visit(const NURBSCurve &g);
};

} // namespace SFCGAL

#endif
