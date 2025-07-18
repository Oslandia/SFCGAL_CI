// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/io/wkt.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_PolygonTest)

// Polygon() ;
BOOST_AUTO_TEST_CASE(defaultConstructor)
{
  Polygon const g;
  BOOST_CHECK(g.isEmpty());
  BOOST_CHECK(!g.is3D());
  BOOST_CHECK_EQUAL(g.numInteriorRings(), 0U);
}

// Polygon( const std::vector< LineString > & rings ) ;
// Polygon( const LineString & exteriorRing ) ;
BOOST_AUTO_TEST_CASE(exteriorRingConstructor)
{
  LineString exteriorRing;
  exteriorRing.addPoint(Point(0.0, 0.0));
  exteriorRing.addPoint(Point(1.0, 0.0));
  exteriorRing.addPoint(Point(1.0, 1.0));
  exteriorRing.addPoint(Point(0.0, 1.0));
  exteriorRing.addPoint(Point(0.0, 0.0));

  Polygon g(exteriorRing);
  BOOST_CHECK(!g.isEmpty());
  BOOST_CHECK(!g.is3D());
  BOOST_CHECK_EQUAL(g.numInteriorRings(), 0U);
  BOOST_CHECK_EQUAL(g.exteriorRing().numPoints(), 5U);
}
BOOST_AUTO_TEST_CASE(exteriorRingConstructor3D)
{
  Polygon g;
  g.exteriorRing().addPoint(Point(0.0, 0.0, 2.0));
  g.exteriorRing().addPoint(Point(1.0, 0.0, 2.0));
  g.exteriorRing().addPoint(Point(1.0, 1.0, 2.0));
  g.exteriorRing().addPoint(Point(0.0, 1.0, 2.0));
  g.exteriorRing().addPoint(Point(0.0, 0.0, 2.0));

  BOOST_CHECK(!g.isEmpty());
  BOOST_CHECK(g.is3D());
  BOOST_CHECK_EQUAL(g.numInteriorRings(), 0U);
}

// Polygon( LineString * exteriorRing ) ;
// Polygon( const Triangle & triangle ) ;
BOOST_AUTO_TEST_CASE(testConstructorTriangle)
{
  Polygon g(Triangle(Point(0.0, 0.0), Point(0.0, 1.0), Point(1.0, 1.0)));
  BOOST_CHECK_EQUAL(g.numRings(), 1U);
  BOOST_CHECK_EQUAL(g.exteriorRing().numPoints(), 4U);
}

// Polygon( Polygon const& other ) ;

// Polygon( const CGAL::Polygon_2< Kernel >& other );
// Polygon( const CGAL::Polygon_with_holes_2< Kernel >& poly );

// Polygon& operator = ( const Polygon & other ) ;
//~Polygon() ;

// bool isCounterClockWiseOriented() const;

// void reverse() ;
BOOST_AUTO_TEST_CASE(testReverse)
{
  Polygon g(Triangle(Point(0.0, 0.0), Point(1.0, 0.0), Point(1.0, 1.0)));
  BOOST_CHECK(g.isCounterClockWiseOriented());
  g.reverse();
  BOOST_CHECK(!g.isCounterClockWiseOriented());
}
// TODO same with holes

// inline const LineString &    exteriorRing() const
// inline LineString &          exteriorRing()
// inline void  setExteriorRing( const LineString& ring )
// inline void  setExteriorRing( LineString* ring )
// inline bool                  hasInteriorRings() const
// inline size_t                numInteriorRings() const
// inline const LineString &    interiorRingN( const size_t & n ) const
// inline LineString &          interiorRingN( const size_t & n )
// inline size_t  numRings() const
// inline const LineString &    ringN( const size_t & n ) const
// inline LineString &          ringN( const size_t & n )
// inline void            addInteriorRing( const LineString & ls )
// inline void            addInteriorRing( LineString* ls )
// inline void            addRing( const LineString & ls )
// inline void            addRing( LineString* ls )
// inline iterator       begin() { return _rings.begin() ; }
// inline const_iterator begin() const { return _rings.begin() ; }
// inline iterator       end() { return _rings.end() ; }
// inline const_iterator end() const { return _rings.end() ; }

// CGAL::Polygon_2<Kernel> toPolygon_2( bool fixOrientation = true ) const;
// CGAL::Polygon_with_holes_2<Kernel> toPolygon_with_holes_2( bool
// fixOrientation = true ) const;

//-- Geometry

// virtual Geometry *   Geometry::clone() const = 0 ;
BOOST_AUTO_TEST_CASE(testClone)
{
  LineString exteriorRing;
  exteriorRing.addPoint(Point(0.0, 0.0));
  exteriorRing.addPoint(Point(1.0, 0.0));
  exteriorRing.addPoint(Point(1.0, 1.0));
  exteriorRing.addPoint(Point(0.0, 1.0));
  exteriorRing.addPoint(Point(0.0, 0.0));

  Polygon const            g(exteriorRing);
  std::unique_ptr<Polygon> copy(g.clone());

  BOOST_CHECK(!copy->isEmpty());
  BOOST_CHECK(!copy->is3D());
  BOOST_CHECK_EQUAL(copy->numInteriorRings(), 0U);
  BOOST_CHECK_EQUAL(copy->exteriorRing().numPoints(), 5U);
}

// virtual Geometry*    Geometry::boundary() const ;
BOOST_AUTO_TEST_CASE(testBoundaryEmpty)
{
  std::unique_ptr<Geometry> boundary(Polygon().boundary());
  BOOST_CHECK(boundary->isEmpty());
  BOOST_CHECK(boundary->is<GeometryCollection>());
}
BOOST_AUTO_TEST_CASE(testBoundaryWithoutHoles)
{
  std::string const         wkt("POLYGON ((0 0,0 1,1 1,0 0))");
  std::unique_ptr<Geometry> boundary(io::readWkt(wkt)->boundary());
  BOOST_CHECK(!boundary->isEmpty());
  BOOST_CHECK_EQUAL(boundary->asText(0), "LINESTRING (0 0,0 1,1 1,0 0)");
}
BOOST_AUTO_TEST_CASE(testBoundaryWithHoles)
{
  std::string const wkt("POLYGON ((0 0,0 5,5 5,0 5,0 0),(1 1,2 1,2 2,1 1))");
  std::unique_ptr<Geometry> boundary(io::readWkt(wkt)->boundary());
  BOOST_CHECK(!boundary->isEmpty());
  BOOST_CHECK_EQUAL(
      boundary->asText(0),
      "MULTILINESTRING ((0 0,0 5,5 5,0 5,0 0),(1 1,2 1,2 2,1 1))");
}

// Envelope             Geometry::envelope() const ;

// std::string          Geometry::asText( const int & numDecimals = -1 ) const ;
BOOST_AUTO_TEST_CASE(asTextEmpty)
{
  Polygon const g;
  BOOST_CHECK_EQUAL(g.asText(1), "POLYGON EMPTY");
}
BOOST_AUTO_TEST_CASE(asText2d)
{
  Polygon g;
  g.exteriorRing().addPoint(Point(0.0, 0.0));
  g.exteriorRing().addPoint(Point(1.0, 0.0));
  g.exteriorRing().addPoint(Point(1.0, 1.0));
  g.exteriorRing().addPoint(Point(0.0, 1.0));
  g.exteriorRing().addPoint(Point(0.0, 0.0));

  BOOST_CHECK_EQUAL(g.asText(1),
                    "POLYGON ((0.0 0.0,1.0 0.0,1.0 1.0,0.0 1.0,0.0 0.0))");
}
BOOST_AUTO_TEST_CASE(asText3d)
{
  Polygon g;
  g.exteriorRing().addPoint(Point(0.0, 0.0, 2.0));
  g.exteriorRing().addPoint(Point(1.0, 0.0, 2.0));
  g.exteriorRing().addPoint(Point(1.0, 1.0, 2.0));
  g.exteriorRing().addPoint(Point(0.0, 1.0, 2.0));
  g.exteriorRing().addPoint(Point(0.0, 0.0, 2.0));

  BOOST_CHECK_EQUAL(g.asText(1), "POLYGON Z ((0.0 0.0 2.0,1.0 0.0 2.0,1.0 1.0 "
                                 "2.0,0.0 1.0 2.0,0.0 0.0 2.0))");
}
// virtual std::string  Geometry::geometryType() const = 0 ;
BOOST_AUTO_TEST_CASE(testGeometryType)
{
  Polygon const g;
  BOOST_CHECK_EQUAL(g.geometryType(), "Polygon");
}
// virtual GeometryType Geometry::geometryTypeId() const = 0 ;
BOOST_AUTO_TEST_CASE(testGeometryTypeId)
{
  Polygon const g;
  BOOST_CHECK_EQUAL(g.geometryTypeId(), TYPE_POLYGON);
}
// virtual int          Geometry::dimension() const = 0 ;
BOOST_AUTO_TEST_CASE(testDimension)
{
  Polygon const g;
  BOOST_CHECK_EQUAL(g.dimension(), 2);
}
// virtual int          Geometry::coordinateDimension() const = 0 ;
// virtual bool         Geometry::isEmpty() const = 0 ;
// virtual bool         Geometry::is3D() const = 0 ;
// virtual bool         Geometry::isMeasured() const = 0 ;
// virtual bool         Geometry::isSimple() const = 0 ;

BOOST_AUTO_TEST_CASE(testDropZM)
{
  Polygon emptyPolygon;
  BOOST_CHECK(!emptyPolygon.is3D());
  BOOST_CHECK(!emptyPolygon.isMeasured());
  BOOST_CHECK(!emptyPolygon.dropZ());
  BOOST_CHECK(!emptyPolygon.dropM());

  Polygon polygon2D;
  polygon2D.exteriorRing().addPoint(Point(0.0, 0.0));
  polygon2D.exteriorRing().addPoint(Point(1.0, 0.0));
  polygon2D.exteriorRing().addPoint(Point(1.0, 1.0));
  BOOST_CHECK(!polygon2D.is3D());
  BOOST_CHECK(!polygon2D.isMeasured());
  BOOST_CHECK(!polygon2D.dropZ());
  BOOST_CHECK(!polygon2D.dropM());

  Polygon polygon3D;
  polygon3D.exteriorRing().addPoint(Point(0.0, 0.0, 2.0));
  polygon3D.exteriorRing().addPoint(Point(1.0, 0.0, 2.0));
  polygon3D.exteriorRing().addPoint(Point(1.0, 1.0, 2.0));
  polygon3D.exteriorRing().addPoint(Point(0.0, 1.0, 2.0));
  polygon3D.exteriorRing().addPoint(Point(0.0, 0.0, 2.0));
  BOOST_CHECK(polygon3D.is3D());
  BOOST_CHECK(!polygon3D.isMeasured());
  BOOST_CHECK(!polygon3D.dropM());
  BOOST_CHECK(polygon3D.dropZ());

  BOOST_CHECK_EQUAL(polygon3D.asText(1),
                    "POLYGON ((0.0 0.0,1.0 0.0,1.0 1.0,0.0 1.0,0.0 0.0))");

  BOOST_CHECK(!polygon3D.is3D());
  BOOST_CHECK(!polygon3D.dropZ());
  BOOST_CHECK(!polygon3D.isMeasured());
  BOOST_CHECK(!polygon3D.dropM());

  std::unique_ptr<Geometry> polygonM(
      io::readWkt("POLYGON M ((0 0 4, 0 3 5, 3 3 6, 3 0 7, 0 0 4))").release());
  BOOST_CHECK(!polygonM->is3D());
  BOOST_CHECK(polygonM->isMeasured());
  BOOST_CHECK(!polygonM->dropZ());
  BOOST_CHECK(polygonM->dropM());
  BOOST_CHECK_EQUAL(polygonM->asText(0), "POLYGON ((0 0,0 3,3 3,3 0,0 0))");
  BOOST_CHECK(!polygonM->is3D());
  BOOST_CHECK(!polygonM->isMeasured());
  BOOST_CHECK(!polygonM->dropZ());
  BOOST_CHECK(!polygonM->dropM());

  Polygon polygonZM;
  polygonZM.exteriorRing().addPoint(Point(0.0, 0.0, 2.0, 1.0));
  polygonZM.exteriorRing().addPoint(Point(1.0, 0.0, 2.0, 1.0));
  polygonZM.exteriorRing().addPoint(Point(1.0, 1.0, 2.0, 1.0));
  polygonZM.exteriorRing().addPoint(Point(0.0, 1.0, 2.0, 1.0));
  polygonZM.exteriorRing().addPoint(Point(0.0, 0.0, 2.0, 1.0));
  BOOST_CHECK(polygonZM.is3D());
  BOOST_CHECK(polygonZM.isMeasured());

  BOOST_CHECK(polygonZM.dropM());
  BOOST_CHECK(polygonZM.is3D());
  BOOST_CHECK(!polygonZM.isMeasured());
  BOOST_CHECK_EQUAL(polygonZM.asText(0),
                    "POLYGON Z ((0 0 2,1 0 2,1 1 2,0 1 2,0 0 2))");
  BOOST_CHECK(!polygonZM.dropM());

  BOOST_CHECK(polygonZM.dropZ());
  BOOST_CHECK(!polygonZM.is3D());
  BOOST_CHECK(!polygonZM.isMeasured());
  BOOST_CHECK_EQUAL(polygonZM.asText(0), "POLYGON ((0 0,1 0,1 1,0 1,0 0))");
  BOOST_CHECK(!polygonZM.dropZ());
  BOOST_CHECK(!polygonZM.dropM());
}

BOOST_AUTO_TEST_CASE(testSwapXY)
{
  Polygon emptyPolygon;
  BOOST_CHECK(emptyPolygon.isEmpty());
  emptyPolygon.swapXY();
  BOOST_CHECK(emptyPolygon.isEmpty());

  Polygon polygon2D;
  polygon2D.exteriorRing().addPoint(Point(3.0, 0.0));
  polygon2D.exteriorRing().addPoint(Point(1.0, 0.0));
  polygon2D.exteriorRing().addPoint(Point(1.0, 5.0));
  polygon2D.swapXY();
  BOOST_CHECK_EQUAL(polygon2D.asText(0), "POLYGON ((0 3,0 1,5 1))");

  Polygon polygon3D;
  polygon3D.exteriorRing().addPoint(Point(3.0, 0.0, 2.0));
  polygon3D.exteriorRing().addPoint(Point(1.0, 0.0, 2.0));
  polygon3D.exteriorRing().addPoint(Point(1.0, 1.0, 2.0));
  polygon3D.exteriorRing().addPoint(Point(0.0, 1.0, 2.0));
  polygon3D.exteriorRing().addPoint(Point(7.0, 2.0, 2.0));
  polygon3D.swapXY();
  BOOST_CHECK_EQUAL(polygon3D.asText(1),
                    "POLYGON Z ((0.0 3.0 2.0,0.0 1.0 2.0,1.0 1.0 2.0,1.0 0.0 "
                    "2.0,2.0 7.0 2.0))");

  std::unique_ptr<Geometry> polygonM(
      io::readWkt("POLYGON M ((0 0 4, 0 3 5, 3 3 6, 3 0 7, 0 0 4))").release());
  polygonM->swapXY();
  BOOST_CHECK_EQUAL(polygonM->asText(0),
                    "POLYGON M ((0 0 4,3 0 5,3 3 6,0 3 7,0 0 4))");

  Polygon polygonZM;
  polygonZM.exteriorRing().addPoint(Point(0.0, 0.0, 2.0, 1.0));
  polygonZM.exteriorRing().addPoint(Point(1.0, 0.0, 2.0, 1.0));
  polygonZM.exteriorRing().addPoint(Point(1.0, 1.0, 2.0, 1.0));
  polygonZM.exteriorRing().addPoint(Point(0.0, 1.0, 2.0, 1.0));
  polygonZM.exteriorRing().addPoint(Point(0.0, 0.0, 2.0, 1.0));
  polygonZM.swapXY();
  BOOST_CHECK_EQUAL(polygonZM.asText(0),
                    "POLYGON ZM ((0 0 2 1,0 1 2 1,1 1 2 1,1 0 2 1,0 0 2 1))");
}

// template < typename Derived > inline bool Geometry::is() const
BOOST_AUTO_TEST_CASE(isPolygon)
{
  Polygon const g;
  BOOST_CHECK(g.is<Polygon>());
}

BOOST_AUTO_TEST_CASE(getCoordinateType)
{
  BOOST_CHECK_EQUAL(
      io::readWkt("POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))")->getCoordinateType(),
      CoordinateType::COORDINATE_XY);
  BOOST_CHECK_EQUAL(
      io::readWkt("POLYGON Z((0 0 1, 1 0 1, 1 1 1, 0 1 1, 0 0 1))")
          ->getCoordinateType(),
      CoordinateType::COORDINATE_XYZ);
  BOOST_CHECK_EQUAL(
      io::readWkt("POLYGON M((0 0 2, 1 0 2, 1 1 2, 0 1 2, 0 0 2))")
          ->getCoordinateType(),
      CoordinateType::COORDINATE_XYM);
  BOOST_CHECK_EQUAL(
      io::readWkt("POLYGON ZM((0 0 1 2, 1 0 1 2, 1 1 1 2, 0 1 1 2, 0 0 1 2))")
          ->getCoordinateType(),
      CoordinateType::COORDINATE_XYZM);
}

// template < typename Derived > inline const Derived &  Geometry::as() const
// template < typename Derived > inline Derived &        Geometry::as()

BOOST_AUTO_TEST_SUITE_END()
