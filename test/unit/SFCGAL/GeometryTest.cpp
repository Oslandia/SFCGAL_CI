// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

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

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_GeometryTest)

// virtual ~Geometry();

// MAY BE TESTED IN SUBCLASSES

// virtual Geometry *   Geometry::clone() const = 0 ;
// virtual Geometry*    Geometry::boundary() const ;
// Envelope             Geometry::envelope() const ;
// std::string          Geometry::asText( const int & numDecimals = -1 ) const ;
// virtual std::string  Geometry::geometryType() const = 0 ;
// virtual GeometryType Geometry::geometryTypeId() const = 0 ;
// virtual int          Geometry::dimension() const = 0 ;
// virtual int          Geometry::coordinateDimension() const = 0 ;
// virtual bool         Geometry::isEmpty() const = 0 ;
// virtual bool         Geometry::is3D() const = 0 ;
// virtual bool         Geometry::isMeasured() const = 0 ;
// virtual bool         Geometry::isSimple() const = 0 ;
// template < typename Derived > inline bool Geometry::is() const
// template < typename Derived > inline const Derived &  Geometry::as() const
// template < typename Derived > inline Derived &        Geometry::as()

// TODO
// void round( const long & scale = 1 ) ;

// TODO
// virtual size_t             numGeometries() const ;
// virtual const Geometry  &  geometryN( size_t const& n ) const ;
// virtual Geometry &         geometryN( size_t const& n ) ;

//-- tested in distance
// double distance( const Geometry & other ) const ;
// double distance3D( const Geometry & other ) const ;

//-- tested in geometry visitor test
// virtual void accept( GeometryVisitor & visitor ) = 0 ;
// virtual void accept( ConstGeometryVisitor & visitor ) const = 0 ;

// TODO
// template <class Archive> void serialize( Archive& ar, const unsigned int
// version )

BOOST_AUTO_TEST_CASE(getCoordinateType) {}

BOOST_AUTO_TEST_SUITE_END()
