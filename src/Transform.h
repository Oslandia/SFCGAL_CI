// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TRANSFORM_H_
#define SFCGAL_TRANSFORM_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryVisitor.h"

namespace SFCGAL {

/**
 * Represents a coordinate transform
 */
class SFCGAL_API Transform : public GeometryVisitor {
public:
  virtual ~Transform();
  /**
   * apply transform to a geometry
   */
  virtual void
  transform(Point &p) = 0;

  void
  visit(Point &g) override;
  void
  visit(LineString &g) override;
  void
  visit(Polygon &g) override;
  void
  visit(Triangle &g) override;
  void
  visit(Solid &g) override;
  void
  visit(MultiPoint &g) override;
  void
  visit(MultiLineString &g) override;
  void
  visit(MultiPolygon &g) override;
  void
  visit(MultiSolid &g) override;
  void
  visit(GeometryCollection &g) override;
  void
  visit(PolyhedralSurface &g) override;
  void
  visit(TriangulatedSurface &g) override;
};

} // namespace SFCGAL

#endif
