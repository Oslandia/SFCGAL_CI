// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

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

  void
  visit(const NURBSCurve & /*g*/) override
  {
    type = "NURBSCurve";
  }

  /** String to store the visited geometry type name */
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

BOOST_AUTO_TEST_CASE(testVisitNURBSCurve)
{
  BOOST_CHECK_EQUAL(getTypeWithVisitor<NURBSCurve>(), "NURBSCurve");
}

BOOST_AUTO_TEST_SUITE_END()
