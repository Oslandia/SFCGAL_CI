// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_DETAIL_GETPOINTSVISITOR_H_
#define SFCGAL_DETAIL_GETPOINTSVISITOR_H_

#include "SFCGAL/config.h"

#include "SFCGAL/GeometryVisitor.h"
#include <vector>

namespace SFCGAL {
namespace detail {

/**
 * Get the list of points from a Geometry
 */
class SFCGAL_API GetPointsVisitor : public ConstGeometryVisitor {
public:
  void
  visit(const Point &g) override;
  void
  visit(const LineString &g) override;
  void
  visit(const Polygon &g) override;
  void
  visit(const Triangle &g) override;
  void
  visit(const Solid &g) override;
  void
  visit(const MultiPoint &g) override;
  void
  visit(const MultiLineString &g) override;
  void
  visit(const MultiPolygon &g) override;
  void
  visit(const MultiSolid &g) override;
  void
  visit(const GeometryCollection &g) override;
  void
  visit(const PolyhedralSurface &g) override;
  void
  visit(const TriangulatedSurface &g) override;
  void
  visit(const BezierCurve &g) override;

public:
  typedef std::vector<const Point *>::const_iterator const_iterator;
  std::vector<const Point *>                         points;
};

} // namespace detail
} // namespace SFCGAL

#endif
