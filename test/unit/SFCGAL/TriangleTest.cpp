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
#include "SFCGAL/Kernel.h"

#include "SFCGAL/Envelope.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Triangle.h"
#include <SFCGAL/io/wkt.h>

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_TriangleTest)

// Triangle();
BOOST_AUTO_TEST_CASE(testDefaultConstructor)
{
  Triangle const g;
  BOOST_CHECK(g.isEmpty());
  BOOST_CHECK_EQUAL(g.numGeometries(), 0U);
}
// Triangle( const Kernel::Triangle_2 & triangle ) ;
BOOST_AUTO_TEST_CASE(testConstructorTriangle_2)
{
  Kernel::Point_2 const a(0.0, 0.0);
  Kernel::Point_2 const b(1.0, 0.0);
  Kernel::Point_2 const c(1.0, 1.0);

  Triangle g(Kernel::Triangle_2(a, b, c));
  BOOST_CHECK(!g.isEmpty());
  BOOST_CHECK(!g.is3D());
  BOOST_CHECK_EQUAL(g.numGeometries(), 1U);

  BOOST_CHECK_EQUAL(g.vertex(0).toPoint_2(), a);
  BOOST_CHECK_EQUAL(g.vertex(1).toPoint_2(), b);
  BOOST_CHECK_EQUAL(g.vertex(2).toPoint_2(), c);
}
// Triangle( const Kernel::Triangle_3 & triangle ) ;
BOOST_AUTO_TEST_CASE(testConstructorTriangle_3)
{
  Kernel::Point_3 const a(0.0, 0.0, 1.0);
  Kernel::Point_3 const b(1.0, 0.0, 2.0);
  Kernel::Point_3 const c(1.0, 1.0, 3.0);

  Triangle g(Kernel::Triangle_3(a, b, c));
  BOOST_CHECK(!g.isEmpty());
  BOOST_CHECK(g.is3D());

  BOOST_CHECK_EQUAL(g.vertex(0).toPoint_3(), a);
  BOOST_CHECK_EQUAL(g.vertex(1).toPoint_3(), b);
  BOOST_CHECK_EQUAL(g.vertex(2).toPoint_3(), c);
}
// Triangle( const Point & p, const Point & q, const Point & r );
// Triangle( const Triangle & other );
// Triangle& operator = ( const Triangle & other );
//~Triangle();

// void                   reverse() ;
BOOST_AUTO_TEST_CASE(testReverse)
{
  Kernel::Point_3 const a(0.0, 0.0, 1.0);
  Kernel::Point_3 const b(1.0, 0.0, 2.0);
  Kernel::Point_3 const c(1.0, 1.0, 3.0);

  Triangle g(Kernel::Triangle_3(a, b, c));
  g.reverse();

  BOOST_CHECK_EQUAL(g.vertex(0).toPoint_3(), a);
  BOOST_CHECK_EQUAL(g.vertex(1).toPoint_3(), c);
  BOOST_CHECK_EQUAL(g.vertex(2).toPoint_3(), b);
}

// Polygon               toPolygon() const ;

// inline const Point & vertex( const int & i ) const
// inline Point &       vertex( const int & i )
// inline Kernel::Triangle_2 toTriangle_2() const
// inline Kernel::Triangle_3 toTriangle_3() const

// template < int D > inline typename TypeForDimension<D>::Triangle
// toTriangle_d() const

//-- Geometry tests

// virtual Geometry *   Geometry::clone() const = 0 ;
BOOST_AUTO_TEST_CASE(testClone)
{
  Kernel::Point_3 const a(0.0, 0.0, 1.0);
  Kernel::Point_3 const b(1.0, 0.0, 2.0);
  Kernel::Point_3 const c(1.0, 1.0, 3.0);

  Triangle const g(Kernel::Triangle_3(a, b, c));

  std::unique_ptr<Geometry> copy(g.clone());
  BOOST_REQUIRE(copy->is<Triangle>());
  BOOST_CHECK_EQUAL(copy->asText(0), "TRIANGLE Z ((0 0 1,1 0 2,1 1 3,0 0 1))");
}

// virtual Geometry*    Geometry::boundary() const ;
BOOST_AUTO_TEST_CASE(testBoundary)
{
  Kernel::Point_2 const a(0.0, 0.0);
  Kernel::Point_2 const b(1.0, 0.0);
  Kernel::Point_2 const c(1.0, 1.0);

  Triangle const g(Kernel::Triangle_2(a, b, c));

  std::unique_ptr<Geometry> boundary(g.boundary());
  BOOST_CHECK_EQUAL(boundary->asText(0), "LINESTRING (0 0,1 0,1 1,0 0)");
}

// Envelope             Geometry::envelope() const ;
BOOST_AUTO_TEST_CASE(testEnvelope_empty)
{
  Triangle const g;

  Envelope const bbox = g.envelope();
  BOOST_CHECK(bbox.isEmpty());
}
BOOST_AUTO_TEST_CASE(testEnvelope_2d)
{
  Kernel::Point_2 const a(0.0, -1.0);
  Kernel::Point_2 const b(0.5, 0.2);
  Kernel::Point_2 const c(1.0, 1.5);

  Triangle const g(Kernel::Triangle_2(a, b, c));

  Envelope const bbox = g.envelope();
  BOOST_CHECK(!bbox.isEmpty());
  BOOST_CHECK(!bbox.is3D());

  BOOST_CHECK_EQUAL(bbox.xMin(), 0.0);
  BOOST_CHECK_EQUAL(bbox.xMax(), 1.0);
  BOOST_CHECK_EQUAL(bbox.yMin(), -1.0);
  BOOST_CHECK_EQUAL(bbox.yMax(), 1.5);
}
BOOST_AUTO_TEST_CASE(testEnvelope_3d)
{
  Kernel::Point_3 const a(0.0, -1.0, 2.0);
  Kernel::Point_3 const b(0.5, 0.2, 4.0);
  Kernel::Point_3 const c(1.0, 1.5, 8.0);

  Triangle const g(Kernel::Triangle_3(a, b, c));

  Envelope const bbox = g.envelope();
  BOOST_CHECK(!bbox.isEmpty());
  BOOST_CHECK(bbox.is3D());
  BOOST_CHECK_EQUAL(bbox.xMin(), 0.0);
  BOOST_CHECK_EQUAL(bbox.xMax(), 1.0);
  BOOST_CHECK_EQUAL(bbox.yMin(), -1.0);
  BOOST_CHECK_EQUAL(bbox.yMax(), 1.5);
  BOOST_CHECK_EQUAL(bbox.zMin(), 2.0);
  BOOST_CHECK_EQUAL(bbox.zMax(), 8.0);
}

// std::string          Geometry::asText( const int & numDecimals = -1 ) const ;
BOOST_AUTO_TEST_CASE(asTextEmpty)
{
  Triangle const g;
  BOOST_CHECK_EQUAL(g.asText(1), "TRIANGLE EMPTY");
}
BOOST_AUTO_TEST_CASE(asText2d)
{
  Triangle const g(Point(0.0, 0.0), Point(1.0, 0.0), Point(1.0, 1.0));
  BOOST_CHECK_EQUAL(g.asText(1),
                    "TRIANGLE ((0.0 0.0,1.0 0.0,1.0 1.0,0.0 0.0))");
}
BOOST_AUTO_TEST_CASE(asText3d)
{
  Triangle const g(Point(0.0, 0.0, 2.0), Point(1.0, 0.0, 3.0),
                   Point(1.0, 1.0, 4.0));
  BOOST_CHECK_EQUAL(
      g.asText(1),
      "TRIANGLE Z ((0.0 0.0 2.0,1.0 0.0 3.0,1.0 1.0 4.0,0.0 0.0 2.0))");
}

// virtual std::string  Geometry::geometryType() const = 0 ;
BOOST_AUTO_TEST_CASE(testGeometryType)
{
  BOOST_CHECK_EQUAL(Triangle().geometryType(), "Triangle");
}
// virtual GeometryType Geometry::geometryTypeId() const = 0 ;
BOOST_AUTO_TEST_CASE(testGeometryTypeId)
{
  BOOST_CHECK_EQUAL(Triangle().geometryTypeId(), TYPE_TRIANGLE);
}

// virtual int          Geometry::dimension() const = 0 ;
// virtual int          Geometry::coordinateDimension() const = 0 ;
// virtual bool         Geometry::isEmpty() const = 0 ;
BOOST_AUTO_TEST_CASE(testIsEmpty)
{
  BOOST_CHECK(Point().isEmpty());
  BOOST_CHECK(!Triangle(Point(0.0, 0.0, 2.0), Point(1.0, 0.0, 3.0),
                        Point(1.0, 1.0, 4.0))
                   .isEmpty());
}
// virtual bool         Geometry::is3D() const = 0 ;
// virtual bool         Geometry::isMeasured() const = 0 ;
// virtual bool         Geometry::isSimple() const = 0 ;

BOOST_AUTO_TEST_CASE(testDropZM)
{
  Triangle triangleEmpty;
  BOOST_CHECK(!triangleEmpty.dropZ());
  BOOST_CHECK(!triangleEmpty.dropM());

  Triangle triangle2D(Point(2.0, 3.0), Point(4.0, 5.0), Point(7.0, -1.0));
  BOOST_CHECK(!triangle2D.dropZ());
  BOOST_CHECK(!triangle2D.dropM());

  Triangle triangle3D(Point(2.0, 3.0, 1.0), Point(4.0, 5.0, 1.0), Point(7.0, -1.0, 1.0));
  BOOST_CHECK(!triangle3D.dropM());
  BOOST_CHECK(triangle3D.dropZ());
  BOOST_CHECK_EQUAL(triangle3D.asText(1), "TRIANGLE ((2.0 3.0,4.0 5.0,7.0 -1.0,2.0 3.0))");
  BOOST_CHECK(!triangle3D.dropZ());

  std::unique_ptr<Geometry> pt1M(io::readWkt("POINT M (5 3 4)"));
  std::unique_ptr<Geometry> pt2M(io::readWkt("POINT M (2 5 4)"));
  std::unique_ptr<Geometry> pt3M(io::readWkt("POINT M (9 -1 4)"));
  Triangle triangleM(pt1M->as<Point>(), pt2M->as<Point>(), pt3M->as<Point>());
  BOOST_REQUIRE(triangleM.is<Triangle>());
  BOOST_CHECK(triangleM.isMeasured());
  BOOST_CHECK(!triangleM.is3D());
  BOOST_CHECK(!triangleM.dropZ());
  BOOST_CHECK(triangleM.dropM());
  BOOST_CHECK_EQUAL(triangleM.asText(0), "TRIANGLE ((5 3,2 5,9 -1,5 3))");
  BOOST_CHECK(!triangleM.dropM());
  BOOST_CHECK(!triangleM.dropZ());

  Triangle triangleZM(Point(2.0, 3.0, 1.0, 2.0), Point(4.0, 5.0, 1.0, 5.0), Point(7.0, -1.0, 1.0, 7.0));
  BOOST_CHECK(triangleZM.isMeasured());
  BOOST_CHECK(triangleZM.is3D());

  BOOST_CHECK(triangleZM.dropM());
  BOOST_CHECK(triangleZM.is3D());
  BOOST_CHECK(!triangleZM.isMeasured());
  BOOST_CHECK_EQUAL(triangleZM.asText(0), "TRIANGLE Z ((2 3 1,4 5 1,7 -1 1,2 3 1))");
  BOOST_CHECK(!triangleZM.dropM());

  BOOST_CHECK(triangleZM.dropZ());
  BOOST_CHECK_EQUAL(triangleZM.asText(0), "TRIANGLE ((2 3,4 5,7 -1,2 3))");
  BOOST_CHECK(!triangleZM.dropZ());
  BOOST_CHECK(!triangleZM.dropM());
  BOOST_CHECK(!triangleZM.isMeasured());
  BOOST_CHECK(!triangleZM.is3D());
}

BOOST_AUTO_TEST_CASE(testSwapXY)
{
  Triangle triangleEmpty;
  BOOST_CHECK(triangleEmpty.isEmpty());
  triangleEmpty.swapXY();
  BOOST_CHECK(triangleEmpty.isEmpty());

  Triangle triangle2D(Point(2.0, 3.0), Point(4.0, 5.0), Point(7.0, -1.0));
  triangle2D.swapXY();
  BOOST_CHECK_EQUAL(triangle2D.asText(0), "TRIANGLE ((3 2,5 4,-1 7,3 2))");

  Triangle triangle3D(Point(2.0, 3.0, 1.0), Point(4.0, 5.0, 1.0), Point(7.0, -1.0, 1.0));
  triangle3D.swapXY();
  BOOST_CHECK_EQUAL(triangle3D.asText(0), "TRIANGLE Z ((3 2 1,5 4 1,-1 7 1,3 2 1))");

  std::unique_ptr<Geometry> pt1M(io::readWkt("POINT M (2 3 4)"));
  std::unique_ptr<Geometry> pt2M(io::readWkt("POINT M (4 5 4)"));
  std::unique_ptr<Geometry> pt3M(io::readWkt("POINT M (7 -1 4)"));
  Triangle triangleM(pt1M->as<Point>(), pt2M->as<Point>(), pt3M->as<Point>());
  triangleM.swapXY();
  BOOST_CHECK_EQUAL(triangleM.asText(0), "TRIANGLE M ((3 2 4,5 4 4,-1 7 4,3 2 4))");

  Triangle triangleZM(Point(2.0, 3.0, 1.0, 2.0), Point(4.0, 5.0, 1.0, 5.0), Point(7.0, -1.0, 1.0, 7.0));
  triangleZM.swapXY();
  BOOST_CHECK_EQUAL(triangleZM.asText(0), "TRIANGLE ZM ((3 2 1 2,5 4 1 5,-1 7 1 7,3 2 1 2))");
}

// template < typename Derived > inline bool Geometry::is() const
BOOST_AUTO_TEST_CASE(isTriangle)
{
  Triangle const g;
  BOOST_CHECK(g.is<Triangle>());
}

// template < typename Derived > inline const Derived &  Geometry::as() const
// template < typename Derived > inline Derived &        Geometry::as()
BOOST_AUTO_TEST_CASE(asTriangle)
{
  std::unique_ptr<Geometry> g(new Triangle());
  BOOST_CHECK(g->as<Triangle>().isEmpty());
}

BOOST_AUTO_TEST_SUITE_END()
