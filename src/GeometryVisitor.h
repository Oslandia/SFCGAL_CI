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

  /**
   * Visit a Point geometry
   * @param g The Point geometry to visit
   */
  virtual void
  visit(Point &g) = 0;
  /** Visit a LineString geometry
   * @param g The LineString geometry to visit */
  virtual void
  visit(LineString &g) = 0;
  /** Visit a Polygon geometry
   * @param g The Polygon geometry to visit */
  virtual void
  visit(Polygon &g) = 0;
  /** Visit a Triangle geometry
   * @param g The Triangle geometry to visit */
  virtual void
  visit(Triangle &g) = 0;
  /** Visit a Solid geometry
   * @param g The Solid geometry to visit */
  virtual void
  visit(Solid &g) = 0;
  /** Visit a MultiPoint geometry
   * @param g The MultiPoint geometry to visit */
  virtual void
  visit(MultiPoint &g) = 0;
  /** Visit a MultiLineString geometry
   * @param g The MultiLineString geometry to visit */
  virtual void
  visit(MultiLineString &g) = 0;
  /** Visit a MultiPolygon geometry
   * @param g The MultiPolygon geometry to visit */
  virtual void
  visit(MultiPolygon &g) = 0;
  /** Visit a MultiSolid geometry
   * @param g The MultiSolid geometry to visit */
  virtual void
  visit(MultiSolid &g) = 0;
  /** Visit a GeometryCollection geometry
   * @param g The GeometryCollection geometry to visit */
  virtual void
  visit(GeometryCollection &g) = 0;
  /** Visit a PolyhedralSurface geometry
   * @param g The PolyhedralSurface geometry to visit */
  virtual void
  visit(PolyhedralSurface &g) = 0;
  /** Visit a TriangulatedSurface geometry
   * @param g The TriangulatedSurface geometry to visit */
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

  /** Visit a const Point geometry
   * @param g The const Point geometry to visit */
  virtual void
  visit(const Point &g) = 0;
  /** Visit a const LineString geometry
   * @param g The const LineString geometry to visit */
  virtual void
  visit(const LineString &g) = 0;
  /** Visit a const Polygon geometry
   * @param g The const Polygon geometry to visit */
  virtual void
  visit(const Polygon &g) = 0;
  /** Visit a const Triangle geometry
   * @param g The const Triangle geometry to visit */
  virtual void
  visit(const Triangle &g) = 0;
  /** Visit a const Solid geometry
   * @param g The const Solid geometry to visit */
  virtual void
  visit(const Solid &g) = 0;
  /** Visit a const MultiPoint geometry
   * @param g The const MultiPoint geometry to visit */
  virtual void
  visit(const MultiPoint &g) = 0;
  /** Visit a const MultiLineString geometry
   * @param g The const MultiLineString geometry to visit */
  virtual void
  visit(const MultiLineString &g) = 0;
  /** Visit a const MultiPolygon geometry
   * @param g The const MultiPolygon geometry to visit */
  virtual void
  visit(const MultiPolygon &g) = 0;
  /** Visit a const MultiSolid geometry
   * @param g The const MultiSolid geometry to visit */
  virtual void
  visit(const MultiSolid &g) = 0;
  /** Visit a const GeometryCollection geometry
   * @param g The const GeometryCollection geometry to visit */
  virtual void
  visit(const GeometryCollection &g) = 0;
  /** Visit a const PolyhedralSurface geometry
   * @param g The const PolyhedralSurface geometry to visit */
  virtual void
  visit(const PolyhedralSurface &g) = 0;
  /** Visit a const TriangulatedSurface geometry
   * @param g The const TriangulatedSurface geometry to visit */
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
