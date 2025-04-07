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

#include "SFCGAL/Coordinate.h"
#include "SFCGAL/Exception.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_CoordinateTest)

/// Coordinate() ;
BOOST_AUTO_TEST_CASE(testDefaultConstructor)
{
  Coordinate const g;
  BOOST_CHECK(g.isEmpty());
  BOOST_CHECK_THROW(g.x(), Exception);
  BOOST_CHECK_THROW(g.y(), Exception);
  BOOST_CHECK_THROW(g.z(), Exception);
}

/// Coordinate( const Kernel::FT & x, const Kernel::FT & y ) ;
BOOST_AUTO_TEST_CASE(testXYConstructor)
{
  Coordinate const g(3, 4);
  BOOST_CHECK(!g.isEmpty());
  BOOST_CHECK(!g.is3D());
  BOOST_CHECK_EQUAL(g.x(), 3);
  BOOST_CHECK_EQUAL(g.y(), 4);
}

/// Coordinate( const Kernel::FT & x, const Kernel::FT & y, const Kernel::FT & z
/// ) ;
BOOST_AUTO_TEST_CASE(testXYZConstructor)
{
  Coordinate const g(3, 4, 5);
  BOOST_CHECK(!g.isEmpty());
  BOOST_CHECK(g.is3D());
  BOOST_CHECK_EQUAL(g.x(), 3);
  BOOST_CHECK_EQUAL(g.y(), 4);
  BOOST_CHECK_EQUAL(g.z(), 5);
}

/// Coordinate( const double & x, const double & y ) ;
BOOST_AUTO_TEST_CASE(testXYConstructorDouble)
{
  Coordinate const g(3.0, 4.0);
  BOOST_CHECK(!g.isEmpty());
  BOOST_CHECK(!g.is3D());
  BOOST_CHECK_EQUAL(g.x(), 3);
  BOOST_CHECK_EQUAL(g.y(), 4);
}
/// Coordinate( const double & x, const double & y, const double& z ) ;
BOOST_AUTO_TEST_CASE(testXYZConstructorDouble)
{
  Coordinate const g(3.0, 4.0, 5.0);
  BOOST_CHECK(!g.isEmpty());
  BOOST_CHECK(g.is3D());
  BOOST_CHECK_EQUAL(g.x(), 3);
  BOOST_CHECK_EQUAL(g.y(), 4);
  BOOST_CHECK_EQUAL(g.z(), 5);
  BOOST_CHECK_THROW(Coordinate(std::numeric_limits<double>::infinity(), 0, 0),
                    NonFiniteValueException);
  BOOST_CHECK_THROW(Coordinate(0, std::numeric_limits<double>::infinity(), 0),
                    NonFiniteValueException);
  BOOST_CHECK_THROW(Coordinate(0, 0, std::numeric_limits<double>::infinity()),
                    NonFiniteValueException);
}

/// Coordinate( const CGAL::Point_2< K > & other ):
/// Coordinate( const CGAL::Point_3< K > & other ):
/// Coordinate( const Coordinate & other ) ;
BOOST_AUTO_TEST_CASE(testCopyConstructorEmpty)
{
  Coordinate const  g;
  const Coordinate &copy(g);
  BOOST_CHECK(copy.isEmpty());
}
BOOST_AUTO_TEST_CASE(testCopyConstructorXY)
{
  Coordinate const  g(3, 4);
  const Coordinate &copy(g);
  BOOST_CHECK_EQUAL(copy.x(), 3);
  BOOST_CHECK_EQUAL(copy.y(), 4);
}

/// Coordinate& operator = ( const Coordinate & other ) ;
/// ~Coordinate() ;
/// int          coordinateDimension() const ;
BOOST_AUTO_TEST_CASE(testCoordinateDimensionEmpty)
{
  Coordinate const g;
  BOOST_CHECK_EQUAL(g.coordinateDimension(), 0);
}
BOOST_AUTO_TEST_CASE(testCoordinateDimensionXY)
{
  Coordinate const g(3, 4);
  BOOST_CHECK_EQUAL(g.coordinateDimension(), 2);
}
BOOST_AUTO_TEST_CASE(testCoordinateDimensionXYZ)
{
  Coordinate const g(3, 4, 5);
  BOOST_CHECK_EQUAL(g.coordinateDimension(), 3);
}

/// bool         isEmpty() const ;
/// bool         is3D() const ;
/// Kernel::FT x() const;
/// Kernel::FT y() const;
/// Kernel::FT z() const;

/// Coordinate& round( const Kernel::FT& scaleFactor ) ;
BOOST_AUTO_TEST_CASE(testRoundInteger)
{
  Coordinate g(0.5, 1.5);
  g.round();
  BOOST_CHECK_EQUAL(g.x(), 1);
  BOOST_CHECK_EQUAL(g.y(), 2);
}
BOOST_AUTO_TEST_CASE(testRoundOneDecimal)
{
  Coordinate g(0.52, 1.57);
  g.round(10);
  BOOST_CHECK_CLOSE(g.x(), 0.5, 0.1);
  BOOST_CHECK_CLOSE(g.y(), 1.6, 0.1);

  std::ostringstream oss;
  oss << CGAL::exact(g.x()) << " " << CGAL::exact(g.y());
  BOOST_CHECK_EQUAL(oss.str(), "1/2 8/5"); // 16/10
}

///--- comparators

/// bool operator < ( const Coordinate & other ) const ;
BOOST_AUTO_TEST_CASE(testLessEmpty)
{
  Coordinate const gA;
  Coordinate const gB;
  BOOST_CHECK_THROW((void)(gA < gB), Exception);
}
BOOST_AUTO_TEST_CASE(testLessXY_XY)
{
  BOOST_CHECK(!(Coordinate(0, 0) < Coordinate(0, 0)));
  BOOST_CHECK((Coordinate(0, 0) < Coordinate(1, 0)));
  BOOST_CHECK((Coordinate(1, 0) < Coordinate(1, 1)));
}
BOOST_AUTO_TEST_CASE(testLessXYZ_XYZ)
{
  BOOST_CHECK(!(Coordinate(0, 0, 0) < Coordinate(0, 0, 0)));
  BOOST_CHECK((Coordinate(0, 0, 0) < Coordinate(1, 0, 0)));
  BOOST_CHECK((Coordinate(1, 0, 0) < Coordinate(1, 1, 0)));
  BOOST_CHECK(!(Coordinate(1, 1, 0) < Coordinate(1, 1, 0)));
  BOOST_CHECK((Coordinate(1, 1, 0) < Coordinate(1, 1, 1)));
}
BOOST_AUTO_TEST_CASE(testLessXY_XYZ)
{
  BOOST_CHECK_THROW((void)(Coordinate(0, 0) < Coordinate(0, 0, 0)), Exception);
}

BOOST_AUTO_TEST_CASE(testAlmostEqual)
{
  BOOST_CHECK(Coordinate(0.0, 0.0).almostEqual(Coordinate(0.0, 0.0), 0.0));
  BOOST_CHECK(!Coordinate(0.1, 0.0).almostEqual(Coordinate(0.0, 0.0), 0.0));
  BOOST_CHECK(!Coordinate(0.0, 0.1).almostEqual(Coordinate(0.0, 0.0), 0.0));
  BOOST_CHECK(!Coordinate(0.0, 0.0).almostEqual(Coordinate(0.1, 0.0), 0.0));
  BOOST_CHECK(!Coordinate(0.0, 0.0).almostEqual(Coordinate(0.0, 0.1), 0.0));
  BOOST_CHECK_THROW(
      Coordinate(0.0, 0.0).almostEqual(Coordinate(0.0, 0.0, 0.0), 0.0),
      Exception);
  BOOST_CHECK_THROW(
      Coordinate(0.0, 0.0, 0.0).almostEqual(Coordinate(0.0, 0.0), 0.0),
      Exception);
  BOOST_CHECK(
      Coordinate(0.0, 0.0, 0.0).almostEqual(Coordinate(0.0, 0.0, 0.0), 0.0));
  BOOST_CHECK(
      !Coordinate(0.0, 0.0, 0.1).almostEqual(Coordinate(0.0, 0.0, 0.0), 0.0));
  BOOST_CHECK(
      !Coordinate(0.0, 0.0, 0.0).almostEqual(Coordinate(0.0, 0.0, 0.1), 0.0));
  BOOST_CHECK(!Coordinate(0.0, 0.0, 0.000001)
                   .almostEqual(Coordinate(0.0, 0.0, 0.000003), 0.000001));
}

BOOST_AUTO_TEST_CASE(testDropZ)
{
  Coordinate coord3D(1.0, 2.0, 3.0);
  BOOST_CHECK(!coord3D.isEmpty());
  BOOST_CHECK(coord3D.is3D());
  BOOST_CHECK(coord3D.dropZ());
  BOOST_CHECK(!coord3D.is3D());
  BOOST_CHECK_EQUAL(coord3D.x(), 1.0);
  BOOST_CHECK_EQUAL(coord3D.y(), 2.0);

  Coordinate coord2D(1.0, 4.0);
  BOOST_CHECK(!coord2D.isEmpty());
  BOOST_CHECK(!coord2D.is3D());
  BOOST_CHECK(!coord2D.dropZ());
  BOOST_CHECK_EQUAL(coord2D.x(), 1.0);
  BOOST_CHECK_EQUAL(coord2D.y(), 4.0);

  Coordinate coordEmpty;
  BOOST_CHECK(coordEmpty.isEmpty());
  BOOST_CHECK(!coordEmpty.dropZ());
}

BOOST_AUTO_TEST_CASE(testSwapXY)
{
  Coordinate coord3D(1.0, 2.0, 3.0);
  BOOST_CHECK(!coord3D.isEmpty());
  BOOST_CHECK(coord3D.is3D());
  coord3D.swapXY();
  BOOST_CHECK(coord3D.is3D());
  BOOST_CHECK_EQUAL(coord3D.x(), 2.0);
  BOOST_CHECK_EQUAL(coord3D.y(), 1.0);
  BOOST_CHECK_EQUAL(coord3D.z(), 3.0);

  Coordinate coord2D(1.0, 4.0);
  BOOST_CHECK(!coord2D.isEmpty());
  BOOST_CHECK(!coord2D.is3D());
  coord2D.swapXY();
  BOOST_CHECK(!coord2D.is3D());
  BOOST_CHECK_EQUAL(coord2D.x(), 4.0);
  BOOST_CHECK_EQUAL(coord2D.y(), 1.0);

  Coordinate coordEmpty;
  BOOST_CHECK(coordEmpty.isEmpty());
  coordEmpty.swapXY();
  BOOST_CHECK(coordEmpty.isEmpty());
}

/// bool operator == ( const Coordinate & other ) const ;
/// bool operator != ( const Coordinate & other ) const ;
/// inline Kernel::Vector_2 toVector_2() const
/// inline Kernel::Vector_3 toVector_3() const
/// Kernel::Point_2 toPoint_2() const;
/// Kernel::Point_3 toPoint_3() const;

BOOST_AUTO_TEST_SUITE_END()
