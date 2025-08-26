// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TRANSFORM_H_
#define SFCGAL_TRANSFORM_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryVisitor.h"

namespace SFCGAL {

/**
 * Base interface for applying a coordinate transformation to geometry objects.
 *
 * Implementations must provide transform(Point&) to define how a single point's
 * coordinates are mapped; all other geometry visitor methods apply that point-
 * level transform to the coordinates contained in their respective geometry
 * types in place.
 */
class SFCGAL_API Transform : public GeometryVisitor {
public:
  /**
   * Virtual destructor to allow safe cleanup of derived Transform
   * implementations.
   */
  virtual ~Transform();
  /**
   * Apply the coordinate transformation to a single Point in place.
   *
   * Implementations should update the point's coordinate values according to
   * the desired transform.
   *
   * @param p Point to transform (modified in place).
   */
  virtual void
  transform(Point &p) = 0;

  /**
   * Apply the transformation to a Point geometry (visitor override).
   *
   * The Point is modified in place.
   *
   * @param g Point to transform (modified in place).
   */
  void
  visit(Point &g) override;
  /**
   * Apply the transformation to a LineString geometry (visitor override).
   *
   * Each vertex of the LineString is transformed in place.
   *
   * @param g LineString to transform (modified in place).
   */
  void
  visit(LineString &g) override;
  /**
   * Apply the transformation to a Polygon geometry (visitor override).
   *
   * All rings and their vertices are transformed in place.
   *
   * @param g Polygon to transform (modified in place).
   */
  void
  visit(Polygon &g) override;
  /**
   * Apply the transformation to a Triangle geometry (visitor override).
   *
   * Each vertex of the Triangle is transformed in place.
   *
   * @param g Triangle to transform (modified in place).
   */
  void
  visit(Triangle &g) override;
  /**
   * Apply the transformation to a Solid geometry (visitor override).
   *
   * All constituent faces and vertices of the Solid are transformed in place.
   *
   * @param g Solid to transform (modified in place).
   */
  void
  visit(Solid &g) override;
  /**
   * Apply the transformation to a MultiPoint geometry (visitor override).
   *
   * Every point in the collection is transformed in place.
   *
   * @param g MultiPoint to transform (modified in place).
   */
  void
  visit(MultiPoint &g) override;
  /**
   * Apply the transformation to a MultiLineString geometry (visitor override).
   *
   * All lines and their vertices are transformed in place.
   *
   * @param g MultiLineString to transform (modified in place).
   */
  void
  visit(MultiLineString &g) override;
  /**
   * Apply the transformation to a MultiPolygon geometry (visitor override).
   *
   * All polygons, rings, and vertices are transformed in place.
   *
   * @param g MultiPolygon to transform (modified in place).
   */
  void
  visit(MultiPolygon &g) override;
  /**
   * Apply the transformation to a MultiSolid geometry (visitor override).
   *
   * All solids and their constituent geometry are transformed in place.
   *
   * @param g MultiSolid to transform (modified in place).
   */
  void
  visit(MultiSolid &g) override;
  /**
   * Apply the transformation to a GeometryCollection (visitor override).
   *
   * Each element in the collection is visited and transformed in place.
   *
   * @param g GeometryCollection to transform (modified in place).
   */
  void
  visit(GeometryCollection &g) override;
  /**
   * Apply the transformation to a PolyhedralSurface geometry (visitor
   * override).
   *
   * All faces and their vertices are transformed in place.
   *
   * @param g PolyhedralSurface to transform (modified in place).
   */
  void
  visit(PolyhedralSurface &g) override;
  /**
   * Apply the transformation to a TriangulatedSurface geometry (visitor
   * override).
   *
   * All triangles and their vertices are transformed in place.
   *
   * @param g TriangulatedSurface to transform (modified in place).
   */
  void
  visit(TriangulatedSurface &g) override;
  /**
   * Apply the transformation to a NURBSCurve geometry (visitor override).
   *
   * Control points (and any relevant coordinate data) of the NURBSCurve are
   * transformed in place.
   *
   * @param g NURBSCurve to transform (modified in place).
   */
  void
  visit(NURBSCurve &g) override;
};

} // namespace SFCGAL

#endif
