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
#include <boost/test/unit_test.hpp>

#include <exception>

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
#include "SFCGAL/io/wkt.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_GeometryCollectionTest)

// GeometryCollection() ;
BOOST_AUTO_TEST_CASE(defaultConstructor)
{
  GeometryCollection const g;
  BOOST_CHECK(g.isEmpty());
}

// GeometryCollection( GeometryCollection const& other ) ;
// GeometryCollection& operator = ( const GeometryCollection & other ) ;
// virtual ~GeometryCollection() ;

// virtual size_t              numGeometries() const ;
// virtual const Geometry  &  geometryN( size_t const& n ) const ;
// virtual Geometry &          geometryN( size_t const& n ) ;
// void                      addGeometry( Geometry * geometry ) ;
// void                      addGeometry( Geometry const& geometry ) ;
BOOST_AUTO_TEST_CASE(testAccessors)
{
  GeometryCollection g;
  BOOST_CHECK_EQUAL(g.numGeometries(), 0U);
  BOOST_CHECK_THROW(g.geometryN(0), Exception);

  g.addGeometry(new Point(2.0, 3.0));
  BOOST_CHECK_EQUAL(g.numGeometries(), 1U);
  g.addGeometry(new LineString(Point(0.0, 0.0), Point(1.0, 1.0)));
  BOOST_CHECK_EQUAL(g.numGeometries(), 2U);
  g.addGeometry(
      new Triangle(Point(0.0, 0.0), Point(1.0, 0.0), Point(1.0, 1.0)));
  BOOST_CHECK_EQUAL(g.numGeometries(), 3U);

  BOOST_CHECK_EQUAL(g.geometryN(0).asText(0), "POINT (2 3)");
  BOOST_CHECK_EQUAL(g.geometryN(1).asText(0), "LINESTRING (0 0,1 1)");
  BOOST_CHECK_EQUAL(g.geometryN(2).asText(0), "TRIANGLE ((0 0,1 0,1 1,0 0))");
  BOOST_CHECK_THROW(g.geometryN(3), Exception);

  g.setGeometryN(Point(0.0, 0.0), 1);
  BOOST_CHECK_EQUAL(g.numGeometries(), 3U);
  BOOST_CHECK_EQUAL(g.geometryN(0).asText(0), "POINT (2 3)");
  BOOST_CHECK_EQUAL(g.geometryN(1).asText(0), "POINT (0 0)");
  BOOST_CHECK_EQUAL(g.geometryN(2).asText(0), "TRIANGLE ((0 0,1 0,1 1,0 0))");

  g.setGeometryN(
      new Triangle(Point(3.0, 0.0), Point(4.0, 0.0), Point(4.0, 1.0)), 2);
  BOOST_CHECK_EQUAL(g.numGeometries(), 3U);
  BOOST_CHECK_EQUAL(g.geometryN(0).asText(0), "POINT (2 3)");
  BOOST_CHECK_EQUAL(g.geometryN(1).asText(0), "POINT (0 0)");
  BOOST_CHECK_EQUAL(g.geometryN(2).asText(0), "TRIANGLE ((3 0,4 0,4 1,3 0))");
}

//-- iterators

// inline iterator       begin()
// inline const_iterator begin() const
// inline iterator       end()
// inline const_iterator end() const
BOOST_AUTO_TEST_CASE(testIterators)
{
  GeometryCollection g;
  g.addGeometry(Point(0.0, 0.0));
  g.addGeometry(Point(1.0, 1.0));

  GeometryCollection::const_iterator it = g.begin();

  BOOST_CHECK_EQUAL(it->asText(0), "POINT (0 0)");
  ++it;
  BOOST_CHECK_EQUAL(it->asText(0), "POINT (1 1)");
  ++it;
  BOOST_CHECK(it == g.end());
}

//-- Geometry

// virtual Geometry *   Geometry::clone() const = 0 ;
// virtual Geometry*    Geometry::boundary() const ;
// Envelope             Geometry::envelope() const ;

// std::string          Geometry::asText( const int & numDecimals = -1 ) const ;
BOOST_AUTO_TEST_CASE(asTextEmpty)
{
  GeometryCollection const g;
  BOOST_CHECK_EQUAL(g.asText(1), "GEOMETRYCOLLECTION EMPTY");
}
BOOST_AUTO_TEST_CASE(asText2d)
{
  GeometryCollection g;
  g.addGeometry(Point(2.0, 3.0));
  g.addGeometry(Triangle(Point(0.0, 0.0), Point(1.0, 0.0), Point(1.0, 1.0)));
  BOOST_CHECK_EQUAL(
      g.asText(1), "GEOMETRYCOLLECTION (POINT (2.0 3.0),TRIANGLE ((0.0 0.0,1.0 "
                   "0.0,1.0 1.0,0.0 0.0)))");
}
BOOST_AUTO_TEST_CASE(asText3d)
{
  GeometryCollection g;
  g.addGeometry(Point(2.0, 3.0, 5.0));
  g.addGeometry(Triangle(Point(0.0, 0.0, 6.0), Point(1.0, 0.0, 6.0),
                         Point(1.0, 1.0, 6.0)));
  BOOST_CHECK_EQUAL(
      g.asText(1),
      "GEOMETRYCOLLECTION Z (POINT Z (2.0 3.0 5.0),TRIANGLE Z ((0.0 "
      "0.0 6.0,1.0 0.0 6.0,1.0 1.0 6.0,0.0 0.0 6.0)))");
}

// virtual std::string  Geometry::geometryType() const = 0 ;
BOOST_AUTO_TEST_CASE(testGeometryType)
{
  GeometryCollection const g;
  BOOST_CHECK_EQUAL(g.geometryType(), "GeometryCollection");
}
// virtual GeometryType Geometry::geometryTypeId() const = 0 ;
BOOST_AUTO_TEST_CASE(testGeometryTypeId)
{
  GeometryCollection const g;
  BOOST_CHECK_EQUAL(g.geometryTypeId(), TYPE_GEOMETRYCOLLECTION);
}

// virtual int          Geometry::dimension() const = 0 ;
// virtual int          Geometry::coordinateDimension() const = 0 ;
// virtual bool         Geometry::isEmpty() const = 0 ;
// virtual bool         Geometry::is3D() const = 0 ;
// virtual bool         Geometry::isMeasured() const = 0 ;
// virtual bool         Geometry::isSimple() const = 0 ;

// template < typename Derived > inline bool Geometry::is() const
BOOST_AUTO_TEST_CASE(testIsGeometryCollection)
{
  BOOST_CHECK(GeometryCollection().is<GeometryCollection>());
  BOOST_CHECK(MultiPoint().is<GeometryCollection>());
  BOOST_CHECK(MultiLineString().is<GeometryCollection>());
  BOOST_CHECK(MultiPolygon().is<GeometryCollection>());
  BOOST_CHECK(MultiSolid().is<GeometryCollection>());
}

BOOST_AUTO_TEST_CASE(testDropZM)
{
  GeometryCollection geomEmpty;
  BOOST_CHECK(geomEmpty.isEmpty());
  BOOST_CHECK(!geomEmpty.dropZ());
  BOOST_CHECK(!geomEmpty.dropM());

  GeometryCollection geom2D;
  geom2D.addGeometry(Point(2.0, 3.0));
  geom2D.addGeometry(
      Triangle(Point(0.0, 0.0), Point(1.0, 0.0), Point(1.0, 1.0)));
  BOOST_CHECK(!geom2D.is3D());
  BOOST_CHECK(!geom2D.isMeasured());
  BOOST_CHECK(!geom2D.dropZ());
  BOOST_CHECK(!geom2D.dropM());

  GeometryCollection geom3D;
  geom3D.addGeometry(Point(2.0, 3.0, 5.0));
  geom3D.addGeometry(Triangle(Point(0.0, 0.0, 6.0), Point(1.0, 0.0, 6.0),
                              Point(1.0, 1.0, 6.0)));
  BOOST_CHECK(geom3D.is3D());
  BOOST_CHECK(!geom3D.isMeasured());
  BOOST_CHECK(!geom3D.dropM());
  BOOST_CHECK(geom3D.dropZ());
  BOOST_CHECK_EQUAL(geom3D.asText(1),
                    "GEOMETRYCOLLECTION (POINT (2.0 3.0),TRIANGLE ((0.0 "
                    "0.0,1.0 0.0,1.0 1.0,0.0 0.0)))");
  BOOST_CHECK(!geom3D.dropZ());
  BOOST_CHECK(!geom3D.dropM());

  GeometryCollection geomM;
  geomM.addGeometry(io::readWkt("POINT M (2 3 4)").release());
  geomM.addGeometry(
      io::readWkt("TRIANGLE M ((0 0 1, 5 5 5, 0 5 2, 0 0 1))").release());
  BOOST_REQUIRE(geomM.is<GeometryCollection>());
  BOOST_CHECK(!geomM.is3D());
  BOOST_CHECK(geomM.isMeasured());
  BOOST_CHECK(!geomM.dropZ());
  BOOST_CHECK(geomM.dropM());
  BOOST_CHECK_EQUAL(geomM.asText(1),
                    "GEOMETRYCOLLECTION (POINT (2.0 3.0),"
                    "TRIANGLE ((0.0 0.0,5.0 5.0,0.0 5.0,0.0 0.0)))");
  BOOST_CHECK(!geomM.dropZ());
  BOOST_CHECK(!geomM.dropM());
  BOOST_CHECK(!geomM.is3D());
  BOOST_CHECK(!geomM.isMeasured());

  GeometryCollection geomZM;
  geomZM.addGeometry(Point(2.0, 3.0, 5.0, 4.0));
  geomZM.addGeometry(Triangle(Point(0.0, 0.0, 6.0, 2.0),
                              Point(1.0, 0.0, 6.0, 2.0),
                              Point(1.0, 1.0, 6.0, 2.0)));
  BOOST_CHECK(geomZM.is3D());
  BOOST_CHECK(geomZM.isMeasured());
  BOOST_CHECK(geomZM.dropM());
  BOOST_CHECK(geomZM.is3D());
  BOOST_CHECK(!geomZM.isMeasured());
  BOOST_CHECK_EQUAL(geomZM.asText(0),
                    "GEOMETRYCOLLECTION Z (POINT Z (2 3 5),"
                    "TRIANGLE Z ((0 0 6,1 0 6,1 1 6,0 0 6)))");
  BOOST_CHECK(geomZM.dropZ());
  BOOST_CHECK(!geomZM.is3D());
  BOOST_CHECK(!geomZM.isMeasured());
  BOOST_CHECK_EQUAL(geomZM.asText(0), "GEOMETRYCOLLECTION (POINT (2 3),"
                                      "TRIANGLE ((0 0,1 0,1 1,0 0)))");
  BOOST_CHECK(!geomZM.dropZ());
  BOOST_CHECK(!geomZM.dropM());
}

BOOST_AUTO_TEST_CASE(testSwapXY)
{
  GeometryCollection geomEmpty;
  BOOST_CHECK(geomEmpty.isEmpty());
  geomEmpty.swapXY();
  BOOST_CHECK(geomEmpty.isEmpty());

  GeometryCollection geom2D;
  geom2D.addGeometry(Point(2.0, 3.0));
  geom2D.addGeometry(
      Triangle(Point(0.0, 0.0), Point(1.0, 0.0), Point(1.0, 1.0)));
  geom2D.swapXY();
  BOOST_CHECK_EQUAL(geom2D.asText(1),
                    "GEOMETRYCOLLECTION "
                    "(POINT (3.0 2.0),"
                    "TRIANGLE ((0.0 0.0,0.0 1.0,1.0 1.0,0.0 0.0)))");

  GeometryCollection geom3D;
  geom3D.addGeometry(Point(2.0, 3.0, 5.0));
  geom3D.addGeometry(Triangle(Point(0.0, 0.0, 6.0), Point(1.0, 0.0, 6.0),
                              Point(1.0, 1.0, 6.0)));
  geom3D.swapXY();
  BOOST_CHECK_EQUAL(
      geom3D.asText(1),
      "GEOMETRYCOLLECTION Z "
      "(POINT Z (3.0 2.0 5.0),"
      "TRIANGLE Z ((0.0 0.0 6.0,0.0 1.0 6.0,1.0 1.0 6.0,0.0 0.0 6.0)))");

  GeometryCollection geomM;
  geomM.addGeometry(io::readWkt("POINT M (2 3 4)").release());
  geomM.addGeometry(
      io::readWkt("TRIANGLE M ((0 0 1, 5 5 5, 0 5 2, 0 0 1))").release());
  BOOST_REQUIRE(geomM.is<GeometryCollection>());
  geomM.swapXY();
  BOOST_CHECK_EQUAL(
      geomM.asText(1),
      "GEOMETRYCOLLECTION M "
      "(POINT M (3.0 2.0 4.0),"
      "TRIANGLE M ((0.0 0.0 1.0,5.0 5.0 5.0,5.0 0.0 2.0,0.0 0.0 1.0)))");

  GeometryCollection geomZM;
  geomZM.addGeometry(Point(2.0, 3.0, 5.0, 4.0));
  geomZM.addGeometry(Triangle(Point(0.0, 0.0, 6.0, 2.0),
                              Point(1.0, 0.0, 6.0, 2.0),
                              Point(1.0, 1.0, 6.0, 2.0)));
  geomZM.swapXY();
  BOOST_CHECK_EQUAL(geomZM.asText(0),
                    "GEOMETRYCOLLECTION ZM "
                    "(POINT ZM (3 2 5 4),"
                    "TRIANGLE ZM ((0 0 6 2,0 1 6 2,1 1 6 2,0 0 6 2)))");
}

// template < typename Derived > inline const Derived &  Geometry::as() const
// template < typename Derived > inline Derived &        Geometry::as()

BOOST_AUTO_TEST_SUITE_END()
