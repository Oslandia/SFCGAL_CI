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
   * @param geometry The geometry to visit.
   */
  virtual void
  visit(Geometry &geometry);

  /**
   * Visit a Point geometry
   * @param geometry The Point geometry to visit
   */
  virtual void
  visit(Point &geometry) = 0;
  /** Visit a LineString geometry
   * @param geometry The LineString geometry to visit */
  virtual void
  visit(LineString &geometry) = 0;
  /** Visit a Polygon geometry
   * @param geometry The Polygon geometry to visit */
  virtual void
  visit(Polygon &geometry) = 0;
  /** Visit a Triangle geometry
   * @param geometry The Triangle geometry to visit */
  virtual void
  visit(Triangle &geometry) = 0;
  /** Visit a Solid geometry
   * @param geometry The Solid geometry to visit */
  virtual void
  visit(Solid &geometry) = 0;
  /** Visit a MultiPoint geometry
   * @param geometry The MultiPoint geometry to visit */
  virtual void
  visit(MultiPoint &geometry) = 0;
  /** Visit a MultiLineString geometry
   * @param geometry The MultiLineString geometry to visit */
  virtual void
  visit(MultiLineString &geometry) = 0;
  /** Visit a MultiPolygon geometry
   * @param geometry The MultiPolygon geometry to visit */
  virtual void
  visit(MultiPolygon &geometry) = 0;
  /** Visit a MultiSolid geometry
   * @param geometry The MultiSolid geometry to visit */
  virtual void
  visit(MultiSolid &geometry) = 0;
  /** Visit a GeometryCollection geometry
   * @param geometry The GeometryCollection geometry to visit */
  virtual void
  visit(GeometryCollection &geometry) = 0;
  /** Visit a PolyhedralSurface geometry
   * @param geometry The PolyhedralSurface geometry to visit */
  virtual void
  visit(PolyhedralSurface &geometry) = 0;
  /** Visit a TriangulatedSurface geometry
   * @param geometry The TriangulatedSurface geometry to visit */
  virtual void
  visit(TriangulatedSurface &geometry) = 0;
  /**
   * Visit a NURBSCurve geometry.
   *
   * Default implementation provides no-op; override to handle NURBSCurve.
   *
   * @param geometry The NURBSCurve to visit.
   */
  virtual void
  visit(NURBSCurve &geometry);
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
   * @param geometry The const geometry to visit.
   */
  virtual void
  visit(const Geometry &geometry);

  /** Visit a const Point geometry
   * @param geometry The const Point geometry to visit */
  virtual void
  visit(const Point &geometry) = 0;
  /** Visit a const LineString geometry
   * @param geometry The const LineString geometry to visit */
  virtual void
  visit(const LineString &geometry) = 0;
  /** Visit a const Polygon geometry
   * @param geometry The const Polygon geometry to visit */
  virtual void
  visit(const Polygon &geometry) = 0;
  /** Visit a const Triangle geometry
   * @param geometry The const Triangle geometry to visit */
  virtual void
  visit(const Triangle &geometry) = 0;
  /** Visit a const Solid geometry
   * @param geometry The const Solid geometry to visit */
  virtual void
  visit(const Solid &geometry) = 0;
  /** Visit a const MultiPoint geometry
   * @param geometry The const MultiPoint geometry to visit */
  virtual void
  visit(const MultiPoint &geometry) = 0;
  /** Visit a const MultiLineString geometry
   * @param geometry The const MultiLineString geometry to visit */
  virtual void
  visit(const MultiLineString &geometry) = 0;
  /** Visit a const MultiPolygon geometry
   * @param geometry The const MultiPolygon geometry to visit */
  virtual void
  visit(const MultiPolygon &geometry) = 0;
  /** Visit a const MultiSolid geometry
   * @param geometry The const MultiSolid geometry to visit */
  virtual void
  visit(const MultiSolid &geometry) = 0;
  /** Visit a const GeometryCollection geometry
   * @param geometry The const GeometryCollection geometry to visit */
  virtual void
  visit(const GeometryCollection &geometry) = 0;
  /** Visit a const PolyhedralSurface geometry
   * @param geometry The const PolyhedralSurface geometry to visit */
  virtual void
  visit(const PolyhedralSurface &geometry) = 0;
  /** Visit a const TriangulatedSurface geometry
   * @param geometry The const TriangulatedSurface geometry to visit */
  virtual void
  visit(const TriangulatedSurface &geometry) = 0;
  /**
   * Visit a const NURBSCurve geometry.
   *
   * Default implementation provides no-op; override to handle const NURBSCurve.
   *
   * @param geometry The const NURBSCurve to visit.
   */
  virtual void
  visit(const NURBSCurve &geometry);
};

} // namespace SFCGAL

#endif
