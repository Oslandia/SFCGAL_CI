// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/Transform.h"

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

namespace SFCGAL {

Transform::~Transform() = default;

void
Transform::visit(Point &g)
{
  transform(g);
}

void
Transform::visit(LineString &g)
{
  for (size_t i = 0; i < g.numPoints(); i++) {
    visit(g.pointN(i));
  }
}

void
Transform::visit(Polygon &g)
{
  for (size_t i = 0; i < g.numRings(); i++) {
    visit(g.ringN(i));
  }
}

void
Transform::visit(Triangle &g)
{
  visit(g.vertex(0));
  visit(g.vertex(1));
  visit(g.vertex(2));
}

void
Transform::visit(Solid &g)
{
  for (size_t i = 0; i < g.numShells(); i++) {
    visit(g.shellN(i));
  }
}

void
Transform::visit(MultiPoint &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.pointN(i));
  }
}

void
Transform::visit(MultiLineString &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.lineStringN(i));
  }
}

void
Transform::visit(MultiPolygon &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.polygonN(i));
  }
}

void
Transform::visit(MultiSolid &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    visit(g.solidN(i));
  }
}

void
Transform::visit(GeometryCollection &g)
{
  for (size_t i = 0; i < g.numGeometries(); i++) {
    GeometryVisitor::visit(g.geometryN(i));
  }
}

void
Transform::visit(PolyhedralSurface &g)
{
  for (size_t i = 0; i < g.numPatches(); i++) {
    visit(g.patchN(i));
  }
}

/**
 * @brief Apply the transformation to every patch of a TriangulatedSurface.
 *
 * Iterates all patches in the given TriangulatedSurface and recursively
 * visits each patch so the transformation is applied to their constituent
 * geometries (control vertices/triangles).
 */
void
Transform::visit(TriangulatedSurface &g)
{
  for (size_t i = 0; i < g.numPatches(); i++) {
    visit(g.patchN(i));
  }
}

/**
 * @brief Visit a NURBS curve and apply the transform to each control point.
 *
 * Iterates over the curve's control points and recursively visits each one
 * so the transformation is applied to all control vertices of the NURBS curve.
 *
 * @param g NURBSCurve to traverse; its control points are visited in index
 * order.
 */
void
Transform::visit(NURBSCurve &g)
{
  for (size_t i = 0; i < g.numControlPoints(); i++) {
    visit(g.controlPointN(i));
  }
}

} // namespace SFCGAL
