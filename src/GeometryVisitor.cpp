// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/GeometryVisitor.h"

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

#include <stdexcept>

namespace SFCGAL {

/**
 * @brief Default destructor for GeometryVisitor.
 *
 * Defaulted destructor provided for proper destruction of derived visitor
 * instances.
 */
GeometryVisitor::~GeometryVisitor() = default;

void
GeometryVisitor::visit(Geometry &geometry)
{
  geometry.accept(*this);
}
//
/////
/////
/////
// void GeometryVisitor::visit( MultiPoint & g )
//{
//	for ( size_t i = 0; i < g.numGeometries(); i++ ){
//		visit( g.geometryN(i).as< Point >() );
//	}
//}
//
/////
/////
/////
// void GeometryVisitor::visit( MultiLineString & g )
//{
//	for ( size_t i = 0; i < g.numGeometries(); i++ ){
//		visit( g.geometryN(i).as< LineString >() );
//	}
//}
//
/////
/////
/////
// void GeometryVisitor::visit( MultiPolygon & g )
//{
//	for ( size_t i = 0; i < g.numGeometries(); i++ ){
//		visit( g.geometryN(i).as< Polygon >() );
//	}
//}
//
//
/////
/////
/////
// void GeometryVisitor::visit( GeometryCollection & g )
//{
//	for ( size_t i = 0; i < g.numGeometries(); i++ ){
//		visit( g.geometryN(i) );
//	}
//}
//
/////
/////
/////
// void GeometryVisitor::visit( PolyhedralSurface & g )
//{
//	for ( size_t i = 0; i < g.numPatches(); i++ ){
//		visit( g.patchN(i) );
//	}
//}
//
/////
/////
/////
// void GeometryVisitor::visit( TriangulatedSurface & g )
//{
//	for ( size_t i = 0; i < g.numPatches(); i++ ){
//		visit( g.patchN(i) );
//	}
//}

//---------------- ConstGeometryVisitor

/**
 * @brief Default destructor for ConstGeometryVisitor.
 *
 * Defaulted destructor provided for proper destruction of derived visitor
 * instances.
 */
ConstGeometryVisitor::~ConstGeometryVisitor() = default;

void
ConstGeometryVisitor::visit(const Geometry &geometry)
{
  geometry.accept(*this);
}
//
/////
/////
/////
// void ConstGeometryVisitor::visit( const MultiPoint & g )
//{
//	for ( size_t i = 0; i < g.numGeometries(); i++ ){
//		visit( g.geometryN(i).as< Point >() );
//	}
//}
//
/////
/////
/////
// void ConstGeometryVisitor::visit( const MultiLineString & g )
//{
//	for ( size_t i = 0; i < g.numGeometries(); i++ ){
//		visit( g.geometryN(i).as< LineString >() );
//	}
//}
//
/////
/////
/////
// void ConstGeometryVisitor::visit( const MultiPolygon & g )
//{
//	for ( size_t i = 0; i < g.numGeometries(); i++ ){
//		visit( g.geometryN(i).as< Polygon >() );
//	}
//}
//
//
/////
/////
/////
// void ConstGeometryVisitor::visit( const GeometryCollection & g )
//{
//	for ( size_t i = 0; i < g.numGeometries(); i++ ){
//		visit( g.geometryN(i) );
//	}
//}
//
/////
/////
/////
// void ConstGeometryVisitor::visit( const PolyhedralSurface & g )
//{
//	for ( size_t i = 0; i < g.numPatches(); i++ ){
//		visit( g.patchN(i) );
//	}
//}
//
/////
/////
/////
// void ConstGeometryVisitor::visit( const TriangulatedSurface & g )
//{
//	for ( size_t i = 0; i < g.numPatches(); i++ ){
//		visit( g.patchN(i) );
//	}
//}
//

/**
 * @brief ABI-compatible default visitor for NURBSCurve.
 *
 * This default implementation does not process the NURBSCurve and exists only
 * to preserve ABI compatibility. It throws std::runtime_error; derived
 * visitors should override this method to provide proper NURBS handling
 * (for example, by converting the NURBSCurve to a LineString and visiting it).
 *
 * @throws std::runtime_error Always thrown to indicate the NURBSCurve visitor
 * is not implemented.
 */
void
GeometryVisitor::visit(NURBSCurve & /*geometry*/)
{
  // Default implementation: convert to LineString and visit that
  // Derived classes should override this method for proper NURBS handling
  throw std::runtime_error("NURBSCurve visitor not implemented");
}

/**
 * @brief Default const visitor for NURBSCurve.
 *
 * This ABI-compatible default implementation does not handle NURBSCurve
 * directly. It throws std::runtime_error to force derived visitors to provide
 * proper NURBS handling (the intended fallback behavior would be to convert the
 * curve to a LineString and visit that).
 *
 * @throws std::runtime_error Always thrown with message "NURBSCurve const
 * visitor not implemented".
 */
void
ConstGeometryVisitor::visit(const NURBSCurve & /*geometry*/)
{
  // Default implementation: convert to LineString and visit that
  // Derived classes should override this method for proper NURBS handling
  throw std::runtime_error("NURBSCurve const visitor not implemented");
}

} // namespace SFCGAL
