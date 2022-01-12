// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/GeometryVisitor.h>

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

///
///
///
GeometryVisitor::~GeometryVisitor() = default;

///
///
///
void
GeometryVisitor::visit(Geometry &g)
{
  g.accept(*this);
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
//	for ( size_t i = 0; i < g.numPolygons(); i++ ){
//		visit( g.polygonN(i) );
//	}
//}
//
/////
/////
/////
// void GeometryVisitor::visit( TriangulatedSurface & g )
//{
//	for ( size_t i = 0; i < g.numTriangles(); i++ ){
//		visit( g.triangleN(i) );
//	}
//}

//---------------- ConstGeometryVisitor

///
///
///
ConstGeometryVisitor::~ConstGeometryVisitor() = default;

///
///
///
void
ConstGeometryVisitor::visit(const Geometry &g)
{
  g.accept(*this);
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
//	for ( size_t i = 0; i < g.numPolygons(); i++ ){
//		visit( g.polygonN(i) );
//	}
//}
//
/////
/////
/////
// void ConstGeometryVisitor::visit( const TriangulatedSurface & g )
//{
//	for ( size_t i = 0; i < g.numTriangles(); i++ ){
//		visit( g.triangleN(i) );
//	}
//}
//

} // namespace SFCGAL
