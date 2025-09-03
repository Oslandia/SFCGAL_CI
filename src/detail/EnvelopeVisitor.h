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

namespace SFCGAL {
namespace detail {

/**
 * Get the list of points from a Geometry
 *
 * @todo ConstPointVisitor
 */
class SFCGAL_API EnvelopeVisitor : public ConstGeometryVisitor {
public:
  EnvelopeVisitor(Envelope &envelope_);

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
  void
  visit(const BSplineCurve &g) override;
  void
  visit(const NURBSCurve &g) override;

public:
  Envelope &envelope;
};

} // namespace detail
} // namespace SFCGAL

#endif
