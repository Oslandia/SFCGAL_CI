// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/detail/GetPointsVisitor.h>

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

namespace SFCGAL::detail {

///
///
///
void
GetPointsVisitor::visit(const Point &g)
{
  points.push_back(&g);
}

///
///
///
void
GetPointsVisitor::visit(const LineString &g)
{
  for (size_t i = 0; i < g.numPoints(); i++) {
    visit(g.pointN(i));
  }
}

///
///
///
void
GetPointsVisitor::visit(const Polygon &g)
{
  for (size_t i = 0; i < g.numRings(); i++) {
    visit(g.ringN(i));
  }
}

///
///
///
void
GetPointsVisitor::visit(const Triangle &g)
{
  visit(g.vertex(0));
  visit(g.vertex(1));
  visit(g.vertex(2));
}

///
///
///
void
GetPointsVisitor::visit(const Solid &g)
{
  for (size_t i = 0; i < g.numShells(); i++) {
    visit(g.shellN(i));
  }
}

///
///
///
void
GetPointsVisitor::visit(const MultiPoint &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.pointN(i));
  }
}

///
///
///
void
GetPointsVisitor::visit(const MultiLineString &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.lineStringN(i));
  }
}

///
///
///
void
GetPointsVisitor::visit(const MultiPolygon &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.polygonN(i));
  }
}

///
///
///
void
GetPointsVisitor::visit(const MultiSolid &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.solidN(i));
  }
}

///
///
///
void
GetPointsVisitor::visit(const GeometryCollection &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    g.geometryN(i).accept(*this);
  }
}

///
///
///
void
GetPointsVisitor::visit(const PolyhedralSurface &g)
{
  for (size_t i = 0; i < g.numPolygons(); i++) {
    visit(g.polygonN(i));
  }
}

///
///
///
void
GetPointsVisitor::visit(const TriangulatedSurface &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.geometryN(i));
  }
}

} // namespace SFCGAL
