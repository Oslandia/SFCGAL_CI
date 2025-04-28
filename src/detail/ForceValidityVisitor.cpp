// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/ForceValidityVisitor.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"

namespace SFCGAL::detail {

void
ForceValidityVisitor::visit(Point &g)
{
  g.forceValidityFlag(valid_);
}

void
ForceValidityVisitor::visit(LineString &g)
{
  g.forceValidityFlag(valid_);
  for (size_t i = 0; i < g.numPoints(); i++) {
    visit(g.pointN(i));
  }
}

void
ForceValidityVisitor::visit(Polygon &g)
{
  g.forceValidityFlag(valid_);
  for (size_t i = 0; i < g.numRings(); i++) {
    visit(g.ringN(i));
  }
}

void
ForceValidityVisitor::visit(Triangle &g)
{
  g.forceValidityFlag(valid_);
  visit(g.vertex(0));
  visit(g.vertex(1));
  visit(g.vertex(2));
}

void
ForceValidityVisitor::visit(Solid &g)
{
  g.forceValidityFlag(valid_);
  for (size_t i = 0; i < g.numShells(); i++) {
    visit(g.shellN(i));
  }
}

void
ForceValidityVisitor::visit(MultiPoint &g)
{
  g.forceValidityFlag(valid_);
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.pointN(i));
  }
}

void
ForceValidityVisitor::visit(MultiLineString &g)
{
  g.forceValidityFlag(valid_);
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.lineStringN(i));
  }
}

void
ForceValidityVisitor::visit(MultiPolygon &g)
{
  g.forceValidityFlag(valid_);
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.polygonN(i));
  }
}

void
ForceValidityVisitor::visit(MultiSolid &g)
{
  g.forceValidityFlag(valid_);
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.solidN(i));
  }
}

void
ForceValidityVisitor::visit(GeometryCollection &g)
{
  g.forceValidityFlag(valid_);
  for (size_t i = 0; i < g.numGeometries(); i++) {
    g.geometryN(i).accept(*this);
  }
}

void
ForceValidityVisitor::visit(PolyhedralSurface &g)
{
  g.forceValidityFlag(valid_);
  for (size_t i = 0; i < g.numPatches(); i++) {
    visit(g.patchN(i));
  }
}

void
ForceValidityVisitor::visit(TriangulatedSurface &g)
{
  g.forceValidityFlag(valid_);
  for (size_t i = 0; i < g.numPatches(); i++) {
    visit(g.patchN(i));
  }
}

} // namespace SFCGAL::detail
