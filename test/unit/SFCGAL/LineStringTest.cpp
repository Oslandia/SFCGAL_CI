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
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "SFCGAL/Envelope.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_LineStringTest)

/// LineString() ;
BOOST_AUTO_TEST_CASE(defaultConstructor)
{
  LineString const g;
  BOOST_CHECK(g.isEmpty());
  BOOST_CHECK(!g.is3D());
  BOOST_CHECK(!g.isMeasured());
  BOOST_CHECK_EQUAL(g.numPoints(), 0U);
}

/// LineString( const std::vector< Point > & points ) ;
BOOST_AUTO_TEST_CASE(constructorFromVector)
{
  std::vector<Point> points;
  points.emplace_back(1.0, 1.0);
  points.emplace_back(2.0, 2.0);

  LineString g(points);
  BOOST_REQUIRE_EQUAL(g.numPoints(), 2U);
  BOOST_CHECK_EQUAL(g.pointN(0).x(), 1.0);
  BOOST_CHECK_EQUAL(g.pointN(0).y(), 1.0);
  BOOST_CHECK_EQUAL(g.pointN(1).x(), 2.0);
  BOOST_CHECK_EQUAL(g.pointN(1).y(), 2.0);
}

/// LineString( const Point & startPoint, const Point & endPoint ) ;

BOOST_AUTO_TEST_CASE(twoPointsConstructor)
{
  LineString const g(Point(0.0, 0.0), Point(2.0, 3.0));
  BOOST_CHECK(!g.isEmpty());
  BOOST_CHECK(!g.is3D());
  BOOST_CHECK_EQUAL(g.asText(1), "LINESTRING (0.0 0.0,2.0 3.0)");
}

/// LineString( LineString const& other ) ;
/// LineString& operator = ( const LineString & other ) ;
///~LineString() ;

/// void clear() ;
BOOST_AUTO_TEST_CASE(testClear)
{
  LineString g(Point(0.0, 0.0), Point(1.0, 1.0));
  g.clear();
  BOOST_CHECK(g.isEmpty());
}

BOOST_AUTO_TEST_CASE(testAccessors)
{
  LineString g;
  BOOST_CHECK_EQUAL(g.numGeometries(), 0U);

  g.addPoint(Point(0.0, 0.0));
  BOOST_CHECK_EQUAL(g.numGeometries(), 1U);

  g.addPoint(Point(1.0, 0.0));
  BOOST_CHECK_EQUAL(g.numGeometries(), 1U);

  g.addPoint(Point(1.0, 1.0));
  BOOST_CHECK_EQUAL(g.numGeometries(), 1U);

  g.addPoint(new Point(2.0, 1.0));
  BOOST_CHECK_EQUAL(g.numGeometries(), 1U);

  BOOST_CHECK_EQUAL(g.geometryN(0).asText(0), "LINESTRING (0 0,1 0,1 1,2 1)");
}

/// void reverse() ;
BOOST_AUTO_TEST_CASE(testReverse_empty)
{
  LineString g;
  g.reverse();
  BOOST_CHECK(g.isEmpty());
}
BOOST_AUTO_TEST_CASE(testReverse_twoPoints)
{
  LineString g(Point(0.0, 0.0), Point(1.0, 1.0));
  g.reverse();
  BOOST_CHECK_EQUAL(g.pointN(0).x(), 1.0);
  BOOST_CHECK_EQUAL(g.pointN(0).y(), 1.0);
  BOOST_CHECK_EQUAL(g.pointN(1).x(), 0.0);
  BOOST_CHECK_EQUAL(g.pointN(1).y(), 0.0);
}

/// inline size_t          numPoints() const { return _points.size(); }

/// size_t                 numSegments() const ;
BOOST_AUTO_TEST_CASE(testNumSegments_empty)
{
  LineString const g;
  BOOST_CHECK_EQUAL(g.numSegments(), 0U);
}
BOOST_AUTO_TEST_CASE(testNumSegments_onlyOnePoint)
{
  LineString g;
  g.addPoint(Point(0.0, 0.0));
  BOOST_CHECK_EQUAL(g.numSegments(), 0U);
}
BOOST_AUTO_TEST_CASE(testNumSegments_twoPoints)
{
  LineString const g(Point(0.0, 0.0), Point(1.0, 1.0));
  BOOST_CHECK_EQUAL(g.numSegments(), 1U);
}

/// inline const Point  &  pointN( size_t const& n ) const
/// inline Point &         pointN( size_t const& n )
/// inline const Point &   startPoint() const { return _points.front(); }
/// inline Point &         startPoint() { return _points.front(); }
/// inline const Point &   endPoint() const { return _points.back(); }
/// inline Point &         endPoint() { return _points.back(); }
/// inline void            addPoint( const Point & p )
/// inline void            addPoint( Point* p )

/// bool isClosed() const ;
BOOST_AUTO_TEST_CASE(testIsClosed_empty)
{
  LineString const g;
  BOOST_CHECK(!g.isClosed());
}
BOOST_AUTO_TEST_CASE(testIsClosed_noClosed)
{
  LineString g;
  g.addPoint(Point(0.0, 0.0));
  g.addPoint(Point(1.0, 0.0));
  g.addPoint(Point(1.0, 1.0));
  BOOST_CHECK(!g.isClosed());
}
BOOST_AUTO_TEST_CASE(testIsClosed_triangle)
{
  LineString g;
  g.addPoint(Point(0.0, 0.0));
  g.addPoint(Point(1.0, 0.0));
  g.addPoint(Point(1.0, 1.0));
  g.addPoint(Point(0.0, 0.0));
  BOOST_CHECK(g.isClosed());
}

/// inline iterator       begin()
/// inline const_iterator begin() const
/// inline iterator       end()
/// inline const_iterator end() const
BOOST_AUTO_TEST_CASE(simpleIteratorTest)
{
  LineString g;
  g.addPoint(Point(0.0, 0.0));
  g.addPoint(Point(1.0, 1.0));
  g.addPoint(Point(2.0, 2.0));

  LineString::const_iterator       it  = g.begin();
  LineString::const_iterator const end = g.end();
  BOOST_CHECK_EQUAL(it->asText(3), "POINT (0.000 0.000)");
  ++it;
  BOOST_CHECK(it != end);
  BOOST_CHECK_EQUAL(it->asText(3), "POINT (1.000 1.000)");
  ++it;
  BOOST_CHECK(it != end);
  BOOST_CHECK_EQUAL(it->asText(3), "POINT (2.000 2.000)");
  ++it;
  BOOST_CHECK(it == end);
}

BOOST_AUTO_TEST_CASE(iteratorTests)
{
  std::vector<Point> points;
  points.emplace_back(2.0, 3.0, 1.0);
  points.emplace_back(4.0, 5.0, 1.0);
  points.emplace_back(6.0, 7.0, 1.0);
  points.emplace_back(8.0, 9.0, 1.0);
  LineString g(points);

  LineString::Point_2_const_iterator pi;
  LineString::Point_2_const_iterator pi_end;
  int                                i = 0;

  for (boost::tie(pi, pi_end) = g.points_2(); pi != pi_end; pi++, i++) {
    BOOST_CHECK_EQUAL(*pi, g.pointN(i).toPoint_2());
  }

  LineString::Point_3_const_iterator ppi;
  LineString::Point_3_const_iterator ppi_end;
  i = 0;

  for (boost::tie(ppi, ppi_end) = g.points_3(); ppi != ppi_end; ppi++, i++) {
    BOOST_CHECK_EQUAL(*ppi, g.pointN(i).toPoint_3());
  }
}

/// void reserve( const size_t & n ) ;

/// CGAL::Polygon_2< Kernel > toPolygon_2() const;
BOOST_AUTO_TEST_CASE(testToPolygon_2_Empty)
{
  LineString const              g;
  CGAL::Polygon_2<Kernel> const polygon = g.toPolygon_2();
  BOOST_CHECK(polygon.is_empty());
}
BOOST_AUTO_TEST_CASE(testToPolygon_2_checkOrientation)
{
  LineString g;
  g.addPoint(Point(0.0, 0.0));
  g.addPoint(Point(1.0, 1.0));
  g.addPoint(Point(1.0, 0.0));
  g.addPoint(g.startPoint());

  // keep orientation
  {
    CGAL::Polygon_2<Kernel> const polygon = g.toPolygon_2(false);
    BOOST_CHECK(!polygon.is_empty());
    BOOST_CHECK(polygon.is_clockwise_oriented());
  }
  // fix orientation
  {
    CGAL::Polygon_2<Kernel> const polygon = g.toPolygon_2(true);
    BOOST_CHECK(!polygon.is_empty());
    BOOST_CHECK(polygon.is_counterclockwise_oriented());
  }
}

//-- Geometry tests

// virtual Geometry *   Geometry::clone() const = 0 ;
BOOST_AUTO_TEST_CASE(testClone)
{
  LineString g;
  g.addPoint(Point(0.0, 0.0));
  g.addPoint(Point(1.0, 1.0));
  std::unique_ptr<Geometry> copy(g.clone());
  BOOST_REQUIRE(copy->is<LineString>());
  BOOST_CHECK_EQUAL(copy->as<LineString>().numPoints(), 2U);
}

// virtual Geometry*    Geometry::boundary() const ;
BOOST_AUTO_TEST_CASE(testBoundary_empty)
{
  LineString const          g;
  std::unique_ptr<Geometry> boundary(g.boundary());
  BOOST_CHECK(boundary->isEmpty());
  BOOST_CHECK(boundary->is<GeometryCollection>());
}
BOOST_AUTO_TEST_CASE(testBoundary_3points)
{
  LineString g;
  g.addPoint(Point(0.0, 0.0));
  g.addPoint(Point(1.0, 1.0));
  g.addPoint(Point(2.0, 2.0));

  std::unique_ptr<Geometry> boundary(g.boundary());
  BOOST_CHECK(!boundary->isEmpty());
  BOOST_REQUIRE(boundary->is<MultiPoint>());
  BOOST_CHECK_EQUAL(boundary->numGeometries(), 2U);
}
BOOST_AUTO_TEST_CASE(testBoundary_closed)
{
  LineString g;
  g.addPoint(Point(0.0, 0.0));
  g.addPoint(Point(1.0, 1.0));
  g.addPoint(Point(2.0, 2.0));
  g.addPoint(g.startPoint());

  std::unique_ptr<Geometry> boundary(g.boundary());
  BOOST_CHECK(boundary->isEmpty());
  BOOST_CHECK(boundary->is<GeometryCollection>());
}

// Envelope             Geometry::envelope() const ;
BOOST_AUTO_TEST_CASE(testEnvelope_empty)
{
  BOOST_CHECK(LineString().envelope().isEmpty());
}
BOOST_AUTO_TEST_CASE(testEnvelope_2D)
{
  LineString g;
  g.addPoint(Point(1.0, 5.0));
  g.addPoint(Point(2.0, 7.0));
  g.addPoint(Point(3.0, 9.0));

  Envelope const box = g.envelope();
  BOOST_CHECK(!box.isEmpty());
  BOOST_CHECK(!box.is3D());

  BOOST_CHECK_EQUAL(box.xMin(), 1.0);
  BOOST_CHECK_EQUAL(box.xMax(), 3.0);
  BOOST_CHECK_EQUAL(box.yMin(), 5.0);
  BOOST_CHECK_EQUAL(box.yMax(), 9.0);
}
BOOST_AUTO_TEST_CASE(testEnvelope_3D)
{
  LineString g;
  g.addPoint(Point(1.0, 5.0, 11.0));
  g.addPoint(Point(2.0, 7.0, 15.0));
  g.addPoint(Point(3.0, 9.0, 17.0));

  Envelope const box = g.envelope();
  BOOST_CHECK(!box.isEmpty());
  BOOST_CHECK(box.is3D());

  BOOST_CHECK_EQUAL(box.xMin(), 1.0);
  BOOST_CHECK_EQUAL(box.xMax(), 3.0);
  BOOST_CHECK_EQUAL(box.yMin(), 5.0);
  BOOST_CHECK_EQUAL(box.yMax(), 9.0);
  BOOST_CHECK_EQUAL(box.zMin(), 11.0);
  BOOST_CHECK_EQUAL(box.zMax(), 17.0);
}

// std::string          Geometry::asText( const int & numDecimals = -1 ) const ;
BOOST_AUTO_TEST_CASE(asTextEmpty)
{
  LineString const g;
  BOOST_CHECK_EQUAL(g.asText(1), "LINESTRING EMPTY");
}
BOOST_AUTO_TEST_CASE(asText2d)
{
  LineString const g(Point(2.0, 3.0), Point(4.0, 5.0));
  BOOST_CHECK_EQUAL(g.asText(3), "LINESTRING (2.000 3.000,4.000 5.000)");
}
BOOST_AUTO_TEST_CASE(asText3d)
{
  LineString const g(Point(2.0, 3.0, 7.0), Point(4.0, 5.0, 8.0));
  BOOST_CHECK_EQUAL(g.asText(3),
                    "LINESTRING Z (2.000 3.000 7.000,4.000 5.000 8.000)");
}

// virtual std::string  Geometry::geometryType() const = 0 ;
BOOST_AUTO_TEST_CASE(testGeometryType)
{
  LineString const g;
  BOOST_CHECK_EQUAL(g.geometryType(), "LineString");
}

// virtual GeometryType Geometry::geometryTypeId() const = 0 ;
BOOST_AUTO_TEST_CASE(testGeometryTypeId)
{
  LineString const g;
  BOOST_CHECK_EQUAL(g.geometryTypeId(), TYPE_LINESTRING);
}

// virtual int          Geometry::dimension() const = 0 ;
// virtual int          Geometry::coordinateDimension() const = 0 ;

// virtual bool         Geometry::isEmpty() const = 0 ;
BOOST_AUTO_TEST_CASE(testIsEmpty_true)
{
  LineString const g;
  BOOST_CHECK(g.isEmpty());
}
BOOST_AUTO_TEST_CASE(testIsEmpty_false)
{
  LineString const g(Point(0.0, 0.0), Point(1.0, 1.0));
  BOOST_CHECK(!g.isEmpty());
}

// virtual bool         Geometry::is3D() const = 0 ;
BOOST_AUTO_TEST_CASE(testIs3D_false)
{
  LineString const g(Point(0.0, 0.0), Point(1.0, 1.0));
  BOOST_CHECK(!g.is3D());
}
BOOST_AUTO_TEST_CASE(testIs3D_true)
{
  LineString const g(Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 1.0));
  BOOST_CHECK(g.is3D());
}

// virtual bool         Geometry::isMeasured() const = 0 ;
// virtual bool         Geometry::isSimple() const = 0 ;

// template < typename Derived > inline bool Geometry::is() const
BOOST_AUTO_TEST_CASE(isLineString)
{
  LineString const g;
  BOOST_CHECK(g.is<LineString>());
}

BOOST_AUTO_TEST_CASE(testDropZM)
{
  LineString lineEmpty;
  BOOST_CHECK(lineEmpty.isEmpty());
  BOOST_CHECK(!lineEmpty.is3D());
  BOOST_CHECK(!lineEmpty.isMeasured());
  BOOST_CHECK(!lineEmpty.dropZ());
  BOOST_CHECK(!lineEmpty.dropM());

  LineString line2D(Point(0.0, 0.0), Point(1.0, 1.0));
  BOOST_CHECK(!line2D.is3D());
  BOOST_CHECK(!line2D.isMeasured());
  BOOST_CHECK(!line2D.dropM());
  BOOST_CHECK(!line2D.dropZ());

  LineString line3D(Point(0.0, 0.0, 2.0), Point(1.0, 1.0, 5.0));
  BOOST_CHECK(line3D.is3D());
  BOOST_CHECK(!line3D.isMeasured());
  BOOST_CHECK(line3D.dropZ());
  BOOST_CHECK_EQUAL(line3D.asText(1), "LINESTRING (0.0 0.0,1.0 1.0)");
  BOOST_CHECK(!line3D.is3D());
  BOOST_CHECK(!line3D.dropZ());

  std::unique_ptr<Geometry> lineM(
      io::readWkt("LINESTRING M (0 0 4, 1 1 5, 2 2 6)").release());
  BOOST_CHECK(!lineM->is3D());
  BOOST_CHECK(lineM->isMeasured());
  BOOST_CHECK(!lineM->dropZ());
  BOOST_CHECK(lineM->dropM());
  BOOST_CHECK_EQUAL(lineM->asText(0), "LINESTRING (0 0,1 1,2 2)");
  BOOST_CHECK(!lineM->is3D());
  BOOST_CHECK(!lineM->isMeasured());
  BOOST_CHECK(!lineM->dropZ());
  BOOST_CHECK(!lineM->dropM());

  LineString lineZM(Point(0.0, 0.0, 2.0, 4.0), Point(1.0, 1.0, 5.0, 4.0));
  BOOST_CHECK(lineZM.is3D());
  BOOST_CHECK(lineZM.isMeasured());

  BOOST_CHECK(lineZM.dropM());
  BOOST_CHECK(lineZM.is3D());
  BOOST_CHECK(!lineZM.isMeasured());
  BOOST_CHECK_EQUAL(lineZM.asText(0), "LINESTRING Z (0 0 2,1 1 5)");
  BOOST_CHECK(!lineZM.dropM());

  BOOST_CHECK(lineZM.dropZ());
  BOOST_CHECK(!lineZM.is3D());
  BOOST_CHECK(!lineZM.isMeasured());
  BOOST_CHECK_EQUAL(lineZM.asText(0), "LINESTRING (0 0,1 1)");
  BOOST_CHECK(!lineZM.dropZ());
  BOOST_CHECK(!lineZM.dropM());
}

BOOST_AUTO_TEST_CASE(testSwapXY)
{
  LineString lineEmpty;
  BOOST_CHECK(lineEmpty.isEmpty());
  lineEmpty.swapXY();
  BOOST_CHECK(lineEmpty.isEmpty());

  LineString line2D(Point(3.0, 0.0), Point(1.0, 4.0));
  line2D.swapXY();
  BOOST_CHECK_EQUAL(line2D.asText(0), "LINESTRING (0 3,4 1)");

  LineString line3D(Point(-2.0, 0.0, 2.0), Point(5.0, 1.0, 5.0));
  line3D.swapXY();
  BOOST_CHECK_EQUAL(line3D.asText(0), "LINESTRING Z (0 -2 2,1 5 5)");

  std::unique_ptr<Geometry> lineM(
      io::readWkt("LINESTRING M (0 0 4, 1 1 5, 2 2 6)").release());
  lineM->swapXY();
  BOOST_CHECK_EQUAL(lineM->asText(0), "LINESTRING M (0 0 4,1 1 5,2 2 6)");

  LineString lineZM(Point(3.0, 0.0, 2.0, 4.0), Point(1.0, 23.0, 5.0, 4.0));
  lineZM.swapXY();
  BOOST_CHECK_EQUAL(lineZM.asText(0), "LINESTRING ZM (0 3 2 4,23 1 5 4)");
}

BOOST_AUTO_TEST_CASE(getCoordinateType)
{
  BOOST_CHECK_EQUAL(io::readWkt("LINESTRING (0 0, 1 1)")->getCoordinateType(),
                    CoordinateType::COORDINATE_XY);
  BOOST_CHECK_EQUAL(
      io::readWkt("LINESTRING Z (0 0 1, 1 1 1)")->getCoordinateType(),
      CoordinateType::COORDINATE_XYZ);
  BOOST_CHECK_EQUAL(
      io::readWkt("LINESTRING M (0 0 2, 1 1 2)")->getCoordinateType(),
      CoordinateType::COORDINATE_XYM);
  BOOST_CHECK_EQUAL(
      io::readWkt("LINESTRING ZM (0 0 1 2, 1 1 1 2)")->getCoordinateType(),
      CoordinateType::COORDINATE_XYZM);
}

// template < typename Derived > inline const Derived &  Geometry::as() const
// template < typename Derived > inline Derived &        Geometry::as()

BOOST_AUTO_TEST_SUITE_END()
