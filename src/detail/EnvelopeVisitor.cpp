/**
 *   SFCGAL
 *
 *   Copyright (C) 2012-2013 Oslandia <infos@oslandia.com>
 *   Copyright (C) 2012-2013 IGN (http://www.ign.fr)
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Library General Public
 *   License as published by the Free Software Foundation; either
 *   version 2 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Library General Public License for more details.

 *   You should have received a copy of the GNU Library General Public
 *   License along with this library; if not, see
 <http://www.gnu.org/licenses/>.
 */

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
