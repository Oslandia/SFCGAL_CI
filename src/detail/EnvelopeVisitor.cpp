// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/detail/EnvelopeVisitor.h>

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/MultiSolid.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>

namespace SFCGAL {
namespace detail {

///
///
///
EnvelopeVisitor::EnvelopeVisitor(Envelope &envelope_) : envelope(envelope_) {}

///
///
///
void
EnvelopeVisitor::visit(const Point &g)
{
  envelope.expandToInclude(g.coordinate());
}

///
///
///
void
EnvelopeVisitor::visit(const LineString &g)
{
  for (size_t i = 0; i < g.numPoints(); i++) {
    visit(g.pointN(i));
  }
}

///
///
///
void
EnvelopeVisitor::visit(const Polygon &g)
{
  for (size_t i = 0; i < g.numRings(); i++) {
    visit(g.ringN(i));
  }
}

///
///
///
void
EnvelopeVisitor::visit(const Triangle &g)
{
  for (size_t i = 0; i < 3; i++) {
    visit(g.vertex(i));
  }
}

///
///
///
void
EnvelopeVisitor::visit(const Solid &g)
{
  for (size_t i = 0; i < g.numShells(); i++) {
    visit(g.shellN(i));
  }
}

///
///
///
void
EnvelopeVisitor::visit(const MultiPoint &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.pointN(i));
  }
}

///
///
///
void
EnvelopeVisitor::visit(const MultiLineString &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.lineStringN(i));
  }
}

///
///
///
void
EnvelopeVisitor::visit(const MultiPolygon &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.polygonN(i));
  }
}

///
///
///
void
EnvelopeVisitor::visit(const MultiSolid &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.solidN(i));
  }
}

///
///
///
void
EnvelopeVisitor::visit(const GeometryCollection &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    g.geometryN(i).accept(*this);
  }
}

///
///
///
void
EnvelopeVisitor::visit(const PolyhedralSurface &g)
{
  for (size_t i = 0; i < g.numPolygons(); i++) {
    visit(g.polygonN(i));
  }
}

///
///
///
void
EnvelopeVisitor::visit(const TriangulatedSurface &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.geometryN(i));
  }
}

} // namespace detail
} // namespace SFCGAL
