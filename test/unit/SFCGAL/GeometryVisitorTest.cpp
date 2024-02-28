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

#include <cmath>

#include "SFCGAL/GeometryVisitor.h"

using namespace SFCGAL;

// always after CGAL
using namespace boost::unit_test;

/**
 * get type from Geometry
 */
class DemoVisitorGetType : public ConstGeometryVisitor {
public:
  void
  visit(const Point & /*g*/) override
  {
    type = "Point";
  }
  void
  visit(const LineString & /*g*/) override
  {
    type = "LineString";
  }
  void
  visit(const Polygon & /*g*/) override
  {
    type = "Polygon";
  }
  void
  visit(const Triangle & /*g*/) override
  {
    type = "Triangle";
  }
  void
  visit(const Solid & /*g*/) override
  {
    type = "Solid";
  }
  void
  visit(const MultiPoint & /*g*/) override
  {
    type = "MultiPoint";
  }
  void
  visit(const MultiLineString & /*g*/) override
  {
    type = "MultiLineString";
  }
  void
  visit(const MultiPolygon & /*g*/) override
  {
    type = "MultiPolygon";
  }
  void
  visit(const MultiSolid & /*g*/) override
  {
    type = "MultiSolid";
  }

  void
  visit(const GeometryCollection & /*g*/) override
  {
    type = "GeometryCollection";
  }

  void
  visit(const PolyhedralSurface & /*g*/) override
  {
    type = "PolyhedralSurface";
  }

  void
  visit(const TriangulatedSurface & /*g*/) override
  {
    type = "TriangulatedSurface";
  }

  std::string type;
};

template <typename T>
auto
getTypeWithVisitor() -> std::string
{
  std::unique_ptr<Geometry> geometry(new T());
  DemoVisitorGetType        visitor;
  geometry->accept(visitor);
  return visitor.type;
}

/*
 * base checks (mainly for compilation issues)
 */
BOOST_AUTO_TEST_SUITE(SFCGAL_GeometryVisitorTest)

BOOST_AUTO_TEST_CASE(testVisitPoint)
{
  BOOST_CHECK_EQUAL(getTypeWithVisitor<Point>(), "Point");
}
BOOST_AUTO_TEST_CASE(testVisitLineString)
{
  BOOST_CHECK_EQUAL(getTypeWithVisitor<LineString>(), "LineString");
}
BOOST_AUTO_TEST_CASE(testVisitPolygon)
{
  BOOST_CHECK_EQUAL(getTypeWithVisitor<Polygon>(), "Polygon");
}
BOOST_AUTO_TEST_CASE(testVisitTriangle)
{
  BOOST_CHECK_EQUAL(getTypeWithVisitor<Triangle>(), "Triangle");
}

BOOST_AUTO_TEST_CASE(testVisitMultiPoint)
{
  BOOST_CHECK_EQUAL(getTypeWithVisitor<MultiPoint>(), "MultiPoint");
}
BOOST_AUTO_TEST_CASE(testVisitMultiLineString)
{
  BOOST_CHECK_EQUAL(getTypeWithVisitor<MultiLineString>(), "MultiLineString");
}
BOOST_AUTO_TEST_CASE(testVisitMultiPolygon)
{
  BOOST_CHECK_EQUAL(getTypeWithVisitor<MultiPolygon>(), "MultiPolygon");
}
BOOST_AUTO_TEST_CASE(testVisitMultiSolid)
{
  BOOST_CHECK_EQUAL(getTypeWithVisitor<MultiSolid>(), "MultiSolid");
}
BOOST_AUTO_TEST_CASE(testVisitGeometryCollection)
{
  BOOST_CHECK_EQUAL(getTypeWithVisitor<GeometryCollection>(),
                    "GeometryCollection");
}

BOOST_AUTO_TEST_CASE(testVisitTriangulatedSurface)
{
  BOOST_CHECK_EQUAL(getTypeWithVisitor<TriangulatedSurface>(),
                    "TriangulatedSurface");
}

BOOST_AUTO_TEST_CASE(testVisitPolyhedralSurface)
{
  BOOST_CHECK_EQUAL(getTypeWithVisitor<PolyhedralSurface>(),
                    "PolyhedralSurface");
}

BOOST_AUTO_TEST_CASE(testVisitSolid)
{
  BOOST_CHECK_EQUAL(getTypeWithVisitor<Solid>(), "Solid");
}

BOOST_AUTO_TEST_SUITE_END()
