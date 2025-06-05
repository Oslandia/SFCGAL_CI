// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"

#include "SFCGAL/Envelope.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/io/wkt.h"

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_PointTest)

// Point() ;
BOOST_AUTO_TEST_CASE(defaultConstructor)
{
  Point const g;
  BOOST_CHECK(g.isEmpty());
  BOOST_CHECK(!g.is3D());
  BOOST_CHECK(!g.isMeasured());
  BOOST_CHECK_EQUAL(g.numGeometries(), 0U);

  // no more access to double
  BOOST_CHECK_THROW(g.x(), Exception);
  BOOST_CHECK_THROW(g.y(), Exception);
  BOOST_CHECK_THROW(g.z(), Exception);
  BOOST_CHECK(std::isnan(g.m()));
}

// Point( const Coordinate & coordinate ) ;

// Point( const Kernel::FT & x, const Kernel::FT & y ) ;
BOOST_AUTO_TEST_CASE(xyConstructor)
{
  Point const g(2.0, 3.0);
  BOOST_CHECK(!g.isEmpty());
  BOOST_CHECK(!g.is3D());
  BOOST_CHECK_EQUAL(g.x(), 2.0);
  BOOST_CHECK_EQUAL(g.y(), 3.0);
  BOOST_CHECK_EQUAL(g.z(), 0);
}

// Point( const Kernel::FT & x, const Kernel::FT & y, const Kernel::FT & z ) ;
// Point( const double & x, const double & y, const double & z = NaN() ) ;
BOOST_AUTO_TEST_CASE(xyzConstructor)
{
  Point const g(2.0, 3.0, 4.0);
  BOOST_CHECK(!g.isEmpty());
  BOOST_CHECK(g.is3D());
  BOOST_CHECK_EQUAL(g.x(), 2.0);
  BOOST_CHECK_EQUAL(g.y(), 3.0);
  BOOST_CHECK_EQUAL(g.z(), 4.0);
}

BOOST_AUTO_TEST_CASE(dimensionConstructor)
{
  const double x = 1.0;
  const double y = 2.0;
  const double z = 3.0;
  const double m = 4.0;

  BOOST_CHECK(Point(x, y, z, m, CoordinateType::COORDINATE_XY) == Point(x, y));
  BOOST_CHECK(Point(x, y, z, m, CoordinateType::COORDINATE_XYZ) ==
              Point(x, y, z));
  BOOST_CHECK(Point(x, y, z, m, CoordinateType::COORDINATE_XYZM) ==
              Point(x, y, z, m));
  Point xym{x, y};
  xym.setM(m);
  BOOST_CHECK(Point(x, y, z, m, CoordinateType::COORDINATE_XYM) == xym);
}

// Point( const Kernel::Point_2 & other ) ;
// Point( const Kernel::Point_3 & other ) ;
// Point( const Point & other ) ;
// Point& operator = ( const Point & other ) ;
//~Point() ;

//-- tested in Coordinate
// inline Kernel::RT x() const { return _coordinate.x() ; }
// inline Kernel::RT y() const { return _coordinate.y() ; }
// inline Kernel::RT z() const { return _coordinate.z() ; }

// inline double    m() const { return _m ; }
// inline void      setM( const double & m ) { _m = m ; }
BOOST_AUTO_TEST_CASE(testGetSetM)
{
  Point p(3.0, 4.0);
  BOOST_CHECK(!p.isMeasured());
  BOOST_CHECK(std::isnan(p.m()));
  p.setM(5.0);
  BOOST_CHECK_EQUAL(p.m(), 5.0);
}

// bool operator < ( const Point & other ) const ;
// bool operator == ( const Point & other ) const ;
// bool operator != ( const Point & other ) const ;

// inline Kernel::Vector_2 toVector_2() const
// inline Kernel::Vector_3 toVector_3() const
// inline Kernel::Point_2 toPoint_2() const
// inline Kernel::Point_3 toPoint_3() const
BOOST_AUTO_TEST_CASE(emptyToVector_2)
{
  Point const                  g;
  CGAL::Vector_2<Kernel> const p = g.toVector_2();
  BOOST_CHECK_EQUAL(CGAL::to_double(p.x()), 0.0);
  BOOST_CHECK_EQUAL(CGAL::to_double(p.y()), 0.0);
}
BOOST_AUTO_TEST_CASE(xyToVector_2)
{
  Point const                  g(3.0, 4.0);
  CGAL::Vector_2<Kernel> const p = g.toVector_2();
  BOOST_CHECK_EQUAL(CGAL::to_double(p.x()), 3.0);
  BOOST_CHECK_EQUAL(CGAL::to_double(p.y()), 4.0);
}
BOOST_AUTO_TEST_CASE(xyToVector_3)
{
  Point const                  g(3.0, 4.0);
  CGAL::Vector_3<Kernel> const p = g.toVector_3();

  BOOST_CHECK_EQUAL(CGAL::to_double(p.x()), 3.0);
  BOOST_CHECK_EQUAL(CGAL::to_double(p.y()), 4.0);
  BOOST_CHECK_EQUAL(CGAL::to_double(p.z()), 0.0);
}

// template <int D> typename TypeForDimension<D>::Point toPoint_d() const;
// inline Coordinate &       coordinate() { return _coordinate; }
// inline const Coordinate & coordinate() const { return _coordinate; }

//-- SFCGAL::Geometry

// virtual Geometry *   Geometry::clone() const = 0 ;
BOOST_AUTO_TEST_CASE(testClone)
{
  Point const               p(3.0, 4.0);
  std::unique_ptr<Geometry> copy(p.clone());
  BOOST_REQUIRE(copy->is<Point>());
  BOOST_CHECK_EQUAL(copy->as<Point>().x(), 3.0);
  BOOST_CHECK_EQUAL(copy->as<Point>().y(), 4.0);
}

// virtual Geometry*    Geometry::boundary() const ;
BOOST_AUTO_TEST_CASE(testBoundary)
{
  Point const               p(3.0, 4.0);
  std::unique_ptr<Geometry> boundary(p.boundary());
  BOOST_CHECK(boundary->isEmpty());
  BOOST_CHECK(boundary->is<GeometryCollection>());
}

// Envelope             Geometry::envelope() const ;
BOOST_AUTO_TEST_CASE(testEnvelope_empty)
{
  BOOST_CHECK(Point().envelope().isEmpty());
}
BOOST_AUTO_TEST_CASE(testEnvelope_2D)
{
  Point const    g(3.0, 4.0);
  Envelope const box = g.envelope();
  BOOST_CHECK(!box.isEmpty());
  BOOST_CHECK(!box.is3D());

  BOOST_CHECK_EQUAL(box.xMin(), 3.0);
  BOOST_CHECK_EQUAL(box.xMax(), 3.0);
  BOOST_CHECK_EQUAL(box.yMin(), 4.0);
  BOOST_CHECK_EQUAL(box.yMax(), 4.0);
}
BOOST_AUTO_TEST_CASE(testEnvelope_3D)
{
  Point const    g(3.0, 4.0, 5.0);
  Envelope const box = g.envelope();
  BOOST_CHECK(!box.isEmpty());
  BOOST_CHECK(box.is3D());

  BOOST_CHECK_EQUAL(box.xMin(), 3.0);
  BOOST_CHECK_EQUAL(box.xMax(), 3.0);
  BOOST_CHECK_EQUAL(box.yMin(), 4.0);
  BOOST_CHECK_EQUAL(box.yMax(), 4.0);
  BOOST_CHECK_EQUAL(box.zMin(), 5.0);
  BOOST_CHECK_EQUAL(box.zMax(), 5.0);
}

// std::string          Geometry::asText( const int & numDecimals = -1 ) const ;
BOOST_AUTO_TEST_CASE(asTextEmpty)
{
  Point const g;
  BOOST_CHECK_EQUAL(g.asText(1), "POINT EMPTY");
}
BOOST_AUTO_TEST_CASE(asText2d)
{
  Point const g(2.0, 3.0);
  BOOST_CHECK_EQUAL(g.asText(3), "POINT (2.000 3.000)");
}
BOOST_AUTO_TEST_CASE(asText3d)
{
  Point const g(2.0, 3.0, 4.0);
  BOOST_CHECK_EQUAL(g.asText(3), "POINT Z (2.000 3.000 4.000)");
}

// virtual std::string  Geometry::geometryType() const = 0 ;
BOOST_AUTO_TEST_CASE(testGeometryType)
{
  Point const g;
  BOOST_CHECK_EQUAL(g.geometryType(), "Point");
}
// virtual GeometryType Geometry::geometryTypeId() const = 0 ;
BOOST_AUTO_TEST_CASE(testGeometryTypeId)
{
  Point const g;
  BOOST_CHECK_EQUAL(g.geometryTypeId(), TYPE_POINT);
}

// virtual int          Geometry::dimension() const = 0 ;
BOOST_AUTO_TEST_CASE(testDimension)
{
  Point const g;
  BOOST_CHECK_EQUAL(g.dimension(), 0);
}

// virtual int          Geometry::coordinateDimension() const = 0 ;
BOOST_AUTO_TEST_CASE(testCoordinateDimension)
{
  BOOST_CHECK_EQUAL(Point().coordinateDimension(), 0);
  BOOST_CHECK_EQUAL(Point(2.0, 3.0).coordinateDimension(), 2);
  BOOST_CHECK_EQUAL(Point(2.0, 3.0, 4.0).coordinateDimension(), 3);
}
// virtual bool         Geometry::isEmpty() const = 0 ;
BOOST_AUTO_TEST_CASE(testIsEmpty)
{
  BOOST_CHECK(Point().isEmpty());
  BOOST_CHECK(!Point(2.0, 3.0).isEmpty());
}
// virtual bool         Geometry::is3D() const = 0 ;
BOOST_AUTO_TEST_CASE(testIs3D)
{
  BOOST_CHECK(!Point().is3D());
  BOOST_CHECK(!Point(2.0, 3.0).is3D());
  BOOST_CHECK(Point(2.0, 3.0, 4.0).is3D());
}
// virtual bool         Geometry::isMeasured() const = 0 ;
BOOST_AUTO_TEST_CASE(testIsMeasured)
{
  BOOST_CHECK(!Point().isMeasured());
  BOOST_CHECK(!Point(2.0, 3.0).isMeasured());
  BOOST_CHECK(!Point(2.0, 3.0, 4.0).isMeasured());
  BOOST_CHECK(Point(2.0, 3.0, 4.0, 5.0).isMeasured());
}

BOOST_AUTO_TEST_CASE(testDropZM)
{
  Point ptEmpty;
  BOOST_CHECK(ptEmpty.isEmpty());
  BOOST_CHECK(!ptEmpty.is3D());
  BOOST_CHECK(!ptEmpty.isMeasured());
  BOOST_CHECK(!ptEmpty.dropZ());
  BOOST_CHECK(!ptEmpty.dropM());

  Point pt2D(2.0, 3.0);
  BOOST_CHECK(!pt2D.is3D());
  BOOST_CHECK(!pt2D.isMeasured());
  BOOST_CHECK(!pt2D.dropZ());
  BOOST_CHECK(!ptEmpty.dropM());

  Point pt3D(2.0, 3.0, 4.0);
  BOOST_CHECK(pt3D.is3D());
  BOOST_CHECK(!pt3D.isMeasured());
  BOOST_CHECK(pt3D.dropZ());
  BOOST_CHECK(!pt3D.dropM());
  BOOST_CHECK_EQUAL(pt3D.x(), 2.0);
  BOOST_CHECK_EQUAL(pt3D.y(), 3.0);
  BOOST_CHECK_EQUAL(pt3D.z(), 0.0);
  BOOST_CHECK(!pt3D.is3D());
  BOOST_CHECK(!pt3D.isMeasured());
  BOOST_CHECK(!pt3D.dropZ());

  std::unique_ptr<Geometry> ptM(io::readWkt("POINT M (2 3 4)"));
  BOOST_REQUIRE(ptM->is<Point>());
  BOOST_CHECK(!ptM->is3D());
  BOOST_CHECK(ptM->isMeasured());
  BOOST_CHECK(!ptM->dropZ());
  BOOST_CHECK(ptM->dropM());
  BOOST_CHECK(!ptM->isMeasured());
  BOOST_CHECK_EQUAL(ptM->as<Point>().x(), 2.0);
  BOOST_CHECK_EQUAL(ptM->as<Point>().y(), 3.0);
  BOOST_CHECK_EQUAL(ptM->as<Point>().z(), 0.0);
  BOOST_CHECK(!ptM->dropM());

  Point ptZM(2.0, 3.0, 4.0, 5.0);
  BOOST_CHECK(ptZM.is3D());
  BOOST_CHECK(ptZM.isMeasured());
  BOOST_CHECK(ptZM.dropM());
  BOOST_CHECK_EQUAL(ptZM.x(), 2.0);
  BOOST_CHECK_EQUAL(ptZM.y(), 3.0);
  BOOST_CHECK_EQUAL(ptZM.z(), 4.0);
  BOOST_CHECK(ptZM.is3D());
  BOOST_CHECK(!ptZM.isMeasured());
  BOOST_CHECK(!ptZM.dropM());
  BOOST_CHECK(ptZM.dropZ());
  BOOST_CHECK_EQUAL(ptZM.x(), 2.0);
  BOOST_CHECK_EQUAL(ptZM.y(), 3.0);
  BOOST_CHECK_EQUAL(ptZM.z(), 0.0);
  BOOST_CHECK(!ptZM.dropZ());
  BOOST_CHECK(!ptZM.is3D());
  BOOST_CHECK(!ptZM.isMeasured());
}

BOOST_AUTO_TEST_CASE(testSwapXY)
{
  Point ptEmpty;
  BOOST_CHECK(ptEmpty.isEmpty());
  ptEmpty.swapXY();
  BOOST_CHECK(ptEmpty.isEmpty());

  Point pt2D(2.0, 3.0);
  pt2D.swapXY();
  BOOST_CHECK_EQUAL(pt2D.x(), 3.0);
  BOOST_CHECK_EQUAL(pt2D.y(), 2.0);

  Point pt3D(5.0, 3.0, 4.0);
  pt3D.swapXY();
  BOOST_CHECK_EQUAL(pt3D.x(), 3.0);
  BOOST_CHECK_EQUAL(pt3D.y(), 5.0);

  std::unique_ptr<Geometry> ptM(io::readWkt("POINT M (9 2 4)"));
  ptM->swapXY();
  BOOST_CHECK_EQUAL(ptM->as<Point>().x(), 2.0);
  BOOST_CHECK_EQUAL(ptM->as<Point>().y(), 9.0);

  Point ptZM(-2.0, -3.0, 4.0, 5.0);
  ptZM.swapXY();
  BOOST_CHECK_EQUAL(ptZM.x(), -3.0);
  BOOST_CHECK_EQUAL(ptZM.y(), -2.0);
}

// TODO
// virtual bool         Geometry::isSimple() const = 0 ;

// template < typename Derived > inline bool Geometry::is() const
BOOST_AUTO_TEST_CASE(isPoint)
{
  Point const g;
  BOOST_CHECK(g.is<Point>());
}
// template < typename Derived > inline const Derived &  Geometry::as() const
// template < typename Derived > inline Derived &        Geometry::as()
BOOST_AUTO_TEST_CASE(asPoint)
{
  std::unique_ptr<Geometry> g(new Point());
  BOOST_CHECK(g->as<Point>().isEmpty());
}

// virtual size_t              numGeometries() const ;
BOOST_AUTO_TEST_CASE(testAccessors)
{
  Point const g;
  BOOST_CHECK_EQUAL(g.numGeometries(), 0U);

  Point const g2D(2.0, 3.0);
  BOOST_CHECK_EQUAL(g2D.numGeometries(), 1U);

  Point const g3D(2.0, 3.0, 4.0);
  BOOST_CHECK_EQUAL(g3D.numGeometries(), 1U);
}

BOOST_AUTO_TEST_CASE(getCoordinateType)
{
  BOOST_CHECK_EQUAL(io::readWkt("POINT (0 0)")->getCoordinateType(),
                    CoordinateType::COORDINATE_XY);
  BOOST_CHECK_EQUAL(io::readWkt("POINT Z (0 0 1)")->getCoordinateType(),
                    CoordinateType::COORDINATE_XYZ);
  BOOST_CHECK_EQUAL(io::readWkt("POINT M (0 0 2)")->getCoordinateType(),
                    CoordinateType::COORDINATE_XYM);
  BOOST_CHECK_EQUAL(io::readWkt("POINT ZM (0 0 1 2)")->getCoordinateType(),
                    CoordinateType::COORDINATE_XYZM);
}

BOOST_AUTO_TEST_SUITE_END()
