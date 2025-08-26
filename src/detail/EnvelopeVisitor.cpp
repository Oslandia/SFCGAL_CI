// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/EnvelopeVisitor.h"
#include "SFCGAL/algorithm/convexHull.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"

namespace SFCGAL::detail {

EnvelopeVisitor::EnvelopeVisitor(Envelope &envelope_) : envelope(envelope_) {}

void
EnvelopeVisitor::visit(const Point &g)
{
  envelope.expandToInclude(g.coordinate());
}

void
EnvelopeVisitor::visit(const LineString &g)
{
  for (size_t i = 0; i < g.numPoints(); i++) {
    visit(g.pointN(i));
  }
}

void
EnvelopeVisitor::visit(const Polygon &g)
{
  for (size_t i = 0; i < g.numRings(); i++) {
    visit(g.ringN(i));
  }
}

void
EnvelopeVisitor::visit(const Triangle &g)
{
  for (size_t i = 0; i < 3; i++) {
    visit(g.vertex(i));
  }
}

void
EnvelopeVisitor::visit(const Solid &g)
{
  for (size_t i = 0; i < g.numShells(); i++) {
    visit(g.shellN(i));
  }
}

void
EnvelopeVisitor::visit(const MultiPoint &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.pointN(i));
  }
}

void
EnvelopeVisitor::visit(const MultiLineString &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.lineStringN(i));
  }
}

void
EnvelopeVisitor::visit(const MultiPolygon &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.polygonN(i));
  }
}

void
EnvelopeVisitor::visit(const MultiSolid &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.solidN(i));
  }
}

void
EnvelopeVisitor::visit(const GeometryCollection &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    g.geometryN(i).accept(*this);
  }
}

void
EnvelopeVisitor::visit(const PolyhedralSurface &g)
{
  for (size_t i = 0; i < g.numPatches(); i++) {
    visit(g.patchN(i));
  }
}

void
EnvelopeVisitor::visit(const TriangulatedSurface &g)
{
  for (size_t i = 0; i < g.numPatches(); i++) {
    visit(g.patchN(i));
  }
}

void
EnvelopeVisitor::visit(const NURBSCurve &g)
{
  // Note: For precise bounds, we should evaluate the curve at critical points
  // but control points provide a conservative bound
  for (size_t i = 0; i < g.numControlPoints(); i++) {
    visit(g.controlPointN(i));
  }

  // Optionally, for better precision, sample curve at additional points
  if (!g.isEmpty()) {
    auto      bounds      = g.parameterBounds();
    const int sampleCount = 20; // Reasonable sampling for envelope

    for (int i = 0; i <= sampleCount; ++i) {
      auto  t = bounds.first + (bounds.second - bounds.first) * i / sampleCount;
      Point samplePoint = g.evaluate(t);
      visit(samplePoint);
    }
  }
}

} // namespace SFCGAL::detail
