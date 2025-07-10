// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_DETAIL_FORCEVALIDITY_VISITOR_H_
#define SFCGAL_DETAIL_FORCEVALIDITY_VISITOR_H_

#include "SFCGAL/config.h"

#include "SFCGAL/GeometryVisitor.h"

namespace SFCGAL {
namespace detail {

class SFCGAL_API ForceValidityVisitor : public GeometryVisitor {
public:
  ForceValidityVisitor(bool valid) : valid_(valid) {}
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

private:
  bool valid_;
};

} // namespace detail
} // namespace SFCGAL

#endif
