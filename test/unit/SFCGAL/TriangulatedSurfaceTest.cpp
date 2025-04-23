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
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/Envelope.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/io/wkt.h"
#include "SFCGAL/triangulate/triangulatePolygon.h"

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_TriangulatedSurfaceTest)

// TriangulatedSurface() ;
BOOST_AUTO_TEST_CASE(defaultConstructor)
{
  TriangulatedSurface const g;
  BOOST_CHECK(g.isEmpty());
  BOOST_CHECK_EQUAL(g.numPatchs(), 0U);
  BOOST_CHECK_EQUAL(g.numGeometries(), 0U);
}
// TriangulatedSurface( const std::vector< Triangle > & triangle ) ;
BOOST_AUTO_TEST_CASE(constructorWithTriangles)
{
  std::vector<Triangle> triangles;
  triangles.emplace_back(Point(0.0, 0.0), Point(1.0, 0.0), Point(1.0, 1.0));
  triangles.emplace_back(Point(0.0, 0.0), Point(1.0, 1.0), Point(0.0, 1.0));

  TriangulatedSurface const g(triangles);
  BOOST_CHECK(!g.isEmpty());
  BOOST_CHECK_EQUAL(g.numPatchs(), 2U);
  BOOST_CHECK_EQUAL(g.numGeometries(), 1U);
}

// TriangulatedSurface( TriangulatedSurface const& other ) ;
// TriangulatedSurface& operator = ( const TriangulatedSurface & other ) ;
//~TriangulatedSurface() ;

// inline size_t             numPatchs() const { return _triangles.size(); }
// inline const Triangle  &  patchN( size_t const& n ) const {
// inline Triangle &         patchN( size_t const& n ) {
// inline void               addTriangle( const Triangle & triangle )
// inline void               addTriangle( Triangle * triangle )
// void                      addTriangles( const TriangulatedSurface & other ) ;

// virtual size_t              numGeometries() const ;
// virtual const Triangle  &   geometryN( size_t const& n ) const ;
// virtual Triangle &          geometryN( size_t const& n ) ;

// void reserve( const size_t & n ) ;

//-- iterators

// inline iterator       begin() {
// inline const_iterator begin() const {
// inline iterator       end() {
// inline const_iterator end() const {

//-- helpers

// template < typename K, typename Polyhedron > std::unique_ptr<Polyhedron>
// toPolyhedron_3() const;
//  TODO

//-- Geometry tests

// virtual Geometry *   Geometry::clone() const = 0 ;
BOOST_AUTO_TEST_CASE(testClone)
{
  std::vector<Triangle> triangles;
  triangles.emplace_back(Point(0.0, 0.0), Point(1.0, 0.0), Point(1.0, 1.0));
  triangles.emplace_back(Point(0.0, 0.0), Point(1.0, 1.0), Point(0.0, 1.0));

  TriangulatedSurface const g(triangles);

  std::unique_ptr<Geometry> copy(g.clone());
  BOOST_REQUIRE(copy->is<TriangulatedSurface>());
  BOOST_CHECK_EQUAL(copy->as<TriangulatedSurface>().numPatchs(), 2U);
  BOOST_CHECK_EQUAL(copy->as<TriangulatedSurface>().numGeometries(), 1U);
}

// virtual Geometry*    Geometry::boundary() const ;
BOOST_AUTO_TEST_CASE(testBoundary)
{
  std::vector<Triangle> triangles;
  triangles.emplace_back(Point(0.0, 0.0), Point(1.0, 0.0), Point(1.0, 1.0));
  triangles.emplace_back(Point(0.0, 0.0), Point(1.0, 1.0), Point(0.0, 1.0));

  TriangulatedSurface const g(triangles);
  std::unique_ptr<Geometry> boundary(g.boundary());
  // TODO add algorithm::lineMerge and update
  BOOST_CHECK_EQUAL(
      boundary->asText(0),
      "MULTILINESTRING ((0 0,1 0),(1 0,1 1),(1 1,0 1),(0 1,0 0))");
}
BOOST_AUTO_TEST_CASE(testBoundaryClosed)
{
  Point const a(0.0, 0.0, 0.0);
  Point const b(1.0, 0.0, 0.0);
  Point const c(0.0, 1.0, 0.0);
  Point const d(0.0, 0.0, 1.0);

  std::vector<Triangle> triangles;
  triangles.emplace_back(a, c, b);
  triangles.emplace_back(a, b, d);
  triangles.emplace_back(b, c, d);
  triangles.emplace_back(c, a, d);

  TriangulatedSurface const g(triangles);
  std::unique_ptr<Geometry> boundary(g.boundary());
  BOOST_CHECK(boundary->isEmpty());
}

// Envelope             Geometry::envelope() const ;
BOOST_AUTO_TEST_CASE(testEnvelope)
{
  Point const a(0.0, 0.0, 0.0);
  Point const b(1.0, 0.0, 0.0);
  Point const c(0.0, 1.0, 0.0);
  Point const d(0.0, 0.0, 1.0);

  std::vector<Triangle> triangles;
  triangles.emplace_back(a, c, b);
  triangles.emplace_back(a, b, d);
  triangles.emplace_back(b, c, d);
  triangles.emplace_back(c, a, d);

  TriangulatedSurface const g(triangles);
  Envelope const            bbox = g.envelope();
  BOOST_CHECK_EQUAL(bbox.xMin(), 0.0);
  BOOST_CHECK_EQUAL(bbox.xMax(), 1.0);
  BOOST_CHECK_EQUAL(bbox.yMin(), 0.0);
  BOOST_CHECK_EQUAL(bbox.yMax(), 1.0);
  BOOST_CHECK_EQUAL(bbox.zMin(), 0.0);
  BOOST_CHECK_EQUAL(bbox.zMax(), 1.0);
}
// std::string          Geometry::asText( const int & numDecimals = -1 ) const ;
// TODO

// virtual std::string  Geometry::geometryType() const = 0 ;
BOOST_AUTO_TEST_CASE(testGeometryType)
{
  TriangulatedSurface const g;
  BOOST_CHECK_EQUAL(g.geometryType(), "TriangulatedSurface");
}
// virtual GeometryType Geometry::geometryTypeId() const = 0 ;
BOOST_AUTO_TEST_CASE(testGeometryTypeId)
{
  TriangulatedSurface const g;
  BOOST_CHECK_EQUAL(g.geometryTypeId(), TYPE_TRIANGULATEDSURFACE);
}

// virtual int          Geometry::dimension() const = 0 ;
BOOST_AUTO_TEST_CASE(testDimension)
{
  TriangulatedSurface const g;
  BOOST_CHECK_EQUAL(g.dimension(), 2); // surface
}

// virtual int          Geometry::coordinateDimension() const = 0 ;
// virtual bool         Geometry::isEmpty() const = 0 ;
// virtual bool         Geometry::is3D() const = 0 ;
// virtual bool         Geometry::isMeasured() const = 0 ;
// virtual bool         Geometry::isSimple() const = 0 ;
// template < typename Derived > inline bool Geometry::is() const
BOOST_AUTO_TEST_CASE(isTriangulatedSurface)
{
  TriangulatedSurface const g;
  BOOST_CHECK(g.is<TriangulatedSurface>());
}
// template < typename Derived > inline const Derived &  Geometry::as() const
// template < typename Derived > inline Derived &        Geometry::as()

//-- other tests

BOOST_AUTO_TEST_CASE(polyhedronConversionTest)
{
  // two unit squares sharing a common edge (1,0)-(1,1)
  std::string const gstr =
      "POLYHEDRALSURFACE (((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),"
      "((1 0 0,1 1 0,2 1 0,2 0 0,1 0 0)))";
  // the following surface would generate an exception, since the two polygons
  // have opposite orientations "POLYHEDRALSURFACE (((0 0 0,0 1 0,1 1 0,1 0 0,0
  // 0 0)),((2 0 0,2 1 0,1 1 0,1 0 0,2 0 0)))";
  std::unique_ptr<Geometry> const g(io::readWkt(gstr));

  TriangulatedSurface tri;
  triangulate::triangulatePolygon3D(*g, tri);

  std::unique_ptr<CGAL::Polyhedron_3<Kernel>> poly(
      tri.toPolyhedron_3<Kernel, CGAL::Polyhedron_3<Kernel>>());
  // we check the two squares share a common edge
  BOOST_CHECK_EQUAL(poly->size_of_facets(), 4U);
  BOOST_CHECK_EQUAL(poly->size_of_vertices(), 6U);
}

BOOST_AUTO_TEST_CASE(setTriangleNTest)
{
  std::unique_ptr<Geometry> emptyGeom(io::readWkt("TIN EMPTY"));
  BOOST_CHECK(emptyGeom->is<TriangulatedSurface>());
  BOOST_CHECK(emptyGeom->isEmpty());
  BOOST_CHECK_EQUAL(emptyGeom->numGeometries(), 0);
  BOOST_CHECK_EQUAL(emptyGeom->as<TriangulatedSurface>().numPatchs(), 0);

  std::string const triangulatedSurfaceStr = "TIN Z ("
                                             "((0 0 0, 2 0 2, 1 2 4, 0 0 0)),"
                                             "((2 0 2, 3 2 3, 1 2 4, 2 0 2)),"
                                             "((1 2 4, 3 2 3, 2 4 6, 1 2 4))"
                                             ")";

  std::unique_ptr<Geometry> geom(io::readWkt(triangulatedSurfaceStr));
  BOOST_CHECK(!geom->isEmpty());
  BOOST_CHECK_EQUAL(geom->numGeometries(), 1);
  BOOST_CHECK_EQUAL(geom->as<TriangulatedSurface>().numPatchs(), 3);
  BOOST_CHECK_EQUAL(geom->as<TriangulatedSurface>().patchN(0).asText(0),
                    "TRIANGLE Z ((0 0 0,2 0 2,1 2 4,0 0 0))");
  BOOST_CHECK_EQUAL(geom->as<TriangulatedSurface>().patchN(1).asText(0),
                    "TRIANGLE Z ((2 0 2,3 2 3,1 2 4,2 0 2))");
  BOOST_CHECK_EQUAL(geom->as<TriangulatedSurface>().patchN(2).asText(0),
                    "TRIANGLE Z ((1 2 4,3 2 3,2 4 6,1 2 4))");

  // set new Polygon at index 1 from a Geometry object
  std::string const newTriangleStr =
      "TRIANGLE Z ((0 0 0, 2 0 4, 1 2 2, 0 0 0))";
  std::unique_ptr<Geometry> newGeom(io::readWkt(newTriangleStr));
  geom->setGeometryN(newGeom->clone(), 1);

  BOOST_CHECK_EQUAL(geom->numGeometries(), 1);
  BOOST_CHECK_EQUAL(geom->as<TriangulatedSurface>().numPatchs(), 3);
  BOOST_CHECK_EQUAL(geom->as<TriangulatedSurface>().patchN(0).asText(0),
                    "TRIANGLE Z ((0 0 0,2 0 2,1 2 4,0 0 0))");
  BOOST_CHECK_EQUAL(geom->as<TriangulatedSurface>().patchN(1).asText(), newGeom->asText());
  BOOST_CHECK_EQUAL(geom->as<TriangulatedSurface>().patchN(2).asText(0),
                    "TRIANGLE Z ((1 2 4,3 2 3,2 4 6,1 2 4))");

  // set New Triangle at index 2 from a Triangle
  std::string const newTriangleStr2 =
      "TRIANGLE Z ((0 0 0, 0 15 0, 15 15 15, 0 0 0))";
  std::unique_ptr<Geometry> newGeom2(io::readWkt(newTriangleStr2));
  Triangle *newTriangle2 = dynamic_cast<Triangle *>(newGeom2.get());
  BOOST_CHECK(newTriangle2);
  geom->setGeometryN(newTriangle2->clone(), 2);

  BOOST_CHECK_EQUAL(geom->numGeometries(), 1);
  BOOST_CHECK_EQUAL(geom->as<TriangulatedSurface>().numPatchs(), 3);
  BOOST_CHECK_EQUAL(geom->as<TriangulatedSurface>().patchN(0).asText(0),
                    "TRIANGLE Z ((0 0 0,2 0 2,1 2 4,0 0 0))");
  BOOST_CHECK_EQUAL(geom->as<TriangulatedSurface>().patchN(1).asText(), newGeom->asText());
  BOOST_CHECK_EQUAL(geom->as<TriangulatedSurface>().patchN(2).asText(), newGeom2->asText());
}

BOOST_AUTO_TEST_CASE(dropZMTest)
{
  TriangulatedSurface surfaceEmpty;
  BOOST_CHECK(surfaceEmpty.isEmpty());
  BOOST_CHECK(!surfaceEmpty.is3D());
  BOOST_CHECK(!surfaceEmpty.isMeasured());
  BOOST_CHECK(!surfaceEmpty.dropZ());
  BOOST_CHECK(!surfaceEmpty.dropM());

  std::string const triangulatedSurface3DStr = "TIN Z ("
                                               "((0 0 0, 2 0 2, 1 2 4, 0 0 0)),"
                                               "((2 0 2, 3 2 3, 1 2 4, 2 0 2)),"
                                               "((1 2 4, 3 2 3, 2 4 6, 1 2 4))"
                                               ")";

  std::unique_ptr<Geometry> geom3D(io::readWkt(triangulatedSurface3DStr));
  BOOST_CHECK(!geom3D->isEmpty());
  BOOST_CHECK(geom3D->is3D());
  BOOST_CHECK(geom3D->dropZ());
  BOOST_CHECK_EQUAL(geom3D->numGeometries(), 1);
  BOOST_CHECK_EQUAL(geom3D->as<TriangulatedSurface>().numPatchs(), 3);
  BOOST_CHECK_EQUAL(geom3D->as<TriangulatedSurface>().patchN(0).asText(0),
                    "TRIANGLE ((0 0,2 0,1 2,0 0))");
  BOOST_CHECK_EQUAL(geom3D->as<TriangulatedSurface>().patchN(1).asText(0),
                    "TRIANGLE ((2 0,3 2,1 2,2 0))");
  BOOST_CHECK_EQUAL(geom3D->as<TriangulatedSurface>().patchN(2).asText(0),
                    "TRIANGLE ((1 2,3 2,2 4,1 2))");

  std::string const triangulatedSurfaceMStr =
      "TIN M(((0 0 0, 0 1 0, 1 1 0, 0 0 0)))";
  std::unique_ptr<Geometry> geomM(io::readWkt(triangulatedSurfaceMStr));
  BOOST_CHECK(!geomM->is3D());
  BOOST_CHECK(geomM->isMeasured());
  BOOST_CHECK(!geomM->dropZ());
  BOOST_CHECK(geomM->dropM());
  BOOST_CHECK_EQUAL(geomM->asText(0), "TIN (((0 0,0 1,1 1,0 0)))");
  BOOST_CHECK(!geomM->is3D());
  BOOST_CHECK(!geomM->isMeasured());
  BOOST_CHECK(!geomM->dropZ());
  BOOST_CHECK(!geomM->dropM());

  std::string const triangulatedSurfaceZMStr =
      "TIN ZM (((0 0 10 1, 10 0 15 2, 5 5 12 3, 0 0 10 1)),((5 5 12 3, 10 0 15 "
      "2, 10 10 20 4, 5 5 12 3)))";
  std::unique_ptr<Geometry> geomZM(io::readWkt(triangulatedSurfaceZMStr));
  BOOST_CHECK(geomZM->is3D());
  BOOST_CHECK(geomZM->isMeasured());
  BOOST_CHECK(geomZM->dropM());
  BOOST_CHECK(geomZM->is3D());
  BOOST_CHECK(!geomZM->isMeasured());
  BOOST_CHECK_EQUAL(geomZM->asText(0),
                    "TIN Z (((0 0 10,10 0 15,5 5 12,0 0 10)),((5 5 12,10 0 "
                    "15,10 10 20,5 5 12)))");
  BOOST_CHECK(!geomZM->dropM());
  BOOST_CHECK(geomZM->dropZ());
  BOOST_CHECK_EQUAL(geomZM->asText(0),
                    "TIN (((0 0,10 0,5 5,0 0)),((5 5,10 0,10 10,5 5)))");
  BOOST_CHECK(!geomZM->isMeasured());
  BOOST_CHECK(!geomZM->is3D());
  BOOST_CHECK(!geomZM->dropZ());
  BOOST_CHECK(!geomZM->dropM());
}

BOOST_AUTO_TEST_CASE(swapXYTest)
{
  TriangulatedSurface surfaceEmpty;
  BOOST_CHECK(surfaceEmpty.isEmpty());
  surfaceEmpty.swapXY();
  BOOST_CHECK(surfaceEmpty.isEmpty());

  std::string const         triangulatedSurface3DStr = "TIN Z ("
                                                       "((0 0 0, 2 0 2, 1 2 4, 0 0 0)),"
                                                       "((2 0 2, 3 2 3, 1 2 4, 2 0 2)),"
                                                       "((1 2 4, 3 2 3, 2 4 6, 1 2 4))"
                                                       ")";
  std::unique_ptr<Geometry> geom3D(io::readWkt(triangulatedSurface3DStr));
  geom3D->swapXY();
  BOOST_CHECK_EQUAL(geom3D->asText(0), "TIN Z "
                                       "(((0 0 0,0 2 2,2 1 4,0 0 0)),"
                                       "((0 2 2,2 3 3,2 1 4,0 2 2)),"
                                       "((2 1 4,2 3 3,4 2 6,2 1 4)))");

  std::string const triangulatedSurfaceMStr =
      "TIN M(((0 0 0, 0 1 0, 1 1 0, 0 0 0)))";
  std::unique_ptr<Geometry> geomM(io::readWkt(triangulatedSurfaceMStr));
  geomM->swapXY();
  BOOST_CHECK_EQUAL(geomM->asText(0), "TIN M (((0 0 0,1 0 0,1 1 0,0 0 0)))");

  std::string const triangulatedSurfaceZMStr =
      "TIN ZM (((0 0 10 1, 10 0 15 2, 5 5 12 3, 0 0 10 1)),((5 5 12 3, 10 0 15 "
      "2, 10 10 20 4, 5 5 12 3)))";
  std::unique_ptr<Geometry> geomZM(io::readWkt(triangulatedSurfaceZMStr));
  geomZM->swapXY();
  BOOST_CHECK_EQUAL(geomZM->asText(0),
                    "TIN ZM "
                    "(((0 0 10 1,0 10 15 2,5 5 12 3,0 0 10 1)),"
                    "((5 5 12 3,0 10 15 2,10 10 20 4,5 5 12 3)))");
}

BOOST_AUTO_TEST_CASE(getCoordinateType)
{
  BOOST_CHECK_EQUAL(
      io::readWkt("TIN(((0 0, 1 0, 0 1, 0 0)))")->getCoordinateType(),
      CoordinateType::COORDINATE_XY);
  BOOST_CHECK_EQUAL(
      io::readWkt("TIN Z(((0 0 1, 1 0 1, 0 1 1, 0 0 1)))")->getCoordinateType(),
      CoordinateType::COORDINATE_XYZ);
  BOOST_CHECK_EQUAL(
      io::readWkt("TIN M(((0 0 2, 1 0 2, 0 1 2, 0 0 2)))")->getCoordinateType(),
      CoordinateType::COORDINATE_XYM);
  BOOST_CHECK_EQUAL(
      io::readWkt("TIN ZM(((0 0 1 2, 1 0 1 2, 0 1 1 2, 0 0 1 2)))")
          ->getCoordinateType(),
      CoordinateType::COORDINATE_XYZM);
}

BOOST_AUTO_TEST_SUITE_END()
