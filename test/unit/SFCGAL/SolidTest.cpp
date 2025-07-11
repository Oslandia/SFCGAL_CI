// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/Kernel.h"

#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/algorithm/covers.h"
#include "SFCGAL/io/wkt.h"

#include <boost/test/unit_test.hpp>
#include <memory>
using namespace boost::unit_test;

using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_SolidTest)

// Solid() ;
// Solid( const PolyhedralSurface & exteriorShell ) ;
// Solid( PolyhedralSurface * exteriorShell ) ;
// Solid( const std::vector< PolyhedralSurface > & shells ) ;
// Solid( Solid const& other ) ;
// Solid& operator = ( const Solid & other ) ;
//~Solid() ;

// inline const PolyhedralSurface &    exteriorShell() const { return _shells[0]
// ; } inline PolyhedralSurface &          exteriorShell() { return _shells[0] ;
// } inline size_t                       numInteriorShells() const { return
// _shells.size() - 1 ; } inline const PolyhedralSurface  &   interiorShellN(
// size_t const& n ) const { return _shells[n+1]; } inline PolyhedralSurface &
// interiorShellN( size_t const& n ) { return _shells[n+1]; } inline void
// addInteriorShell( const PolyhedralSurface & shell ) inline void
// addInteriorShell( PolyhedralSurface * shell ) inline size_t  numShells()
// const { inline const PolyhedralSurface &  shellN( const size_t & n ) const {
// inline PolyhedralSurface &        shellN( const size_t & n ) {

//-- iterators

// inline iterator       begin() { return _shells.begin() ; }
// inline const_iterator begin() const { return _shells.begin() ; }
// inline iterator       end() { return _shells.end() ; }
// inline const_iterator end() const { return _shells.end() ; }

//-- helpers

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

//-- other tests

BOOST_AUTO_TEST_CASE(solidReadTest)
{
  // the unit cube where half of a cube has been substracted
  std::string const gstr =
      "SOLID ("
      // exterior shell
      "("
      "((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0))," // front face
      "((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0))," // right face
      "((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0))," // top face
      "((0 0 1,0 1 1,0 1 0,0 0 0,0 0 1))," // left face
      "((1 0 1,1 1 1,0 1 1,0 0 1,1 0 1))," // back face
      "((1 0 0,1 0 1,0 0 1,0 0 0,1 0 0))"  // bottom face
      "),"
      // interior shell, a.k.a. a hole
      "("
      "((0 0 0,0 0.5 0,0.5 0.5 0,0.5 0 0,0 0 0)),"             // front face
      "((0.5 0 0,0.5 0.5 0,0.5 0.5 0.5,0.5 0 0.5,0.5 0 0)),"   // right face
      "((0 0.5 0,0 0.5 0.5,0.5 0.5 0.5,0.5 0.5 0,0 0.5 0)),"   // top face
      "((0 0 0.5,0 0.5 0.5,0 0.5 0,0 0 0,0 0 0.5)),"           // left face
      "((0.5 0 0.5,0.5 0.5 0.5,0 0.5 0.5,0 0 0.5,0.5 0 0.5))," // back face
      "((0.5 0 0,0.5 0 0.5,0 0 0.5,0 0 0,0.5 0 0))"            // bottom face
      ")"
      ")";

  std::unique_ptr<Geometry> g(io::readWkt(gstr));
  BOOST_CHECK_EQUAL(g->as<Solid>().numShells(), 2U);
  BOOST_CHECK_EQUAL(g->as<Solid>().numGeometries(), 1U);
}

BOOST_AUTO_TEST_CASE(solidSetExteriorRingTest)
{
  std::unique_ptr<Solid> emptySolid = std::make_unique<Solid>();
  BOOST_CHECK(emptySolid->isEmpty());

  std::string const polyhedral1Str = "POLYHEDRALSURFACE ("
                                     "((0 0 0, 0 0 1, 0 1 1, 0 1 0, 0 0 0)),"
                                     "((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)),"
                                     "((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)),"
                                     "((1 1 0, 1 1 1, 1 0 1, 1 0 0, 1 1 0)),"
                                     "((0 1 0, 0 1 1, 1 1 1, 1 1 0, 0 1 0)),"
                                     "((0 0 1, 1 0 1, 1 1 1, 0 1 1, 0 0 1))"
                                     ")";

  std::unique_ptr<Solid> solid = std::make_unique<Solid>();
  BOOST_CHECK(solid->isEmpty());

  std::unique_ptr<Geometry> shell1(io::readWkt(polyhedral1Str));
  BOOST_CHECK(!shell1->isEmpty());
  BOOST_CHECK(solid->isEmpty());
  BOOST_CHECK_EQUAL(solid->numGeometries(), 0U);

  solid->setExteriorShell(
      dynamic_cast<PolyhedralSurface *>(shell1.get())->clone());
  BOOST_CHECK_EQUAL(solid->numShells(), 1U);
  BOOST_CHECK_EQUAL(solid->numGeometries(), 1U);
  BOOST_CHECK(!solid->isEmpty());
  BOOST_CHECK(algorithm::covers3D(solid->exteriorShell(), *shell1));
}

BOOST_AUTO_TEST_CASE(solidDropZTest)
{
  std::unique_ptr<Solid> emptySolid = std::make_unique<Solid>();
  BOOST_CHECK(emptySolid->isEmpty());
  BOOST_CHECK(!emptySolid->dropZ());

  std::string const         polyhedral1Str = "POLYHEDRALSURFACE ("
                                             "((0 0 0, 0 0 1, 0 1 1, 0 1 0, 0 0 0)),"
                                             "((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)),"
                                             "((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)),"
                                             "((1 1 0, 1 1 1, 1 0 1, 1 0 0, 1 1 0)),"
                                             "((0 1 0, 0 1 1, 1 1 1, 1 1 0, 0 1 0)),"
                                             "((0 0 1, 1 0 1, 1 1 1, 0 1 1, 0 0 1))"
                                             ")";
  std::unique_ptr<Geometry> shell1(io::readWkt(polyhedral1Str));

  std::unique_ptr<Solid> solid = std::make_unique<Solid>();
  solid->setExteriorShell(
      dynamic_cast<PolyhedralSurface *>(shell1.get())->clone());
  BOOST_CHECK(!solid->isEmpty());
  BOOST_CHECK(solid->is3D());
  BOOST_CHECK(solid->dropZ());

  BOOST_CHECK_EQUAL(solid->asText(1),
                    "SOLID ("
                    "(((0.0 0.0,0.0 0.0,0.0 1.0,0.0 1.0,0.0 0.0)),"
                    "((0.0 0.0,0.0 1.0,1.0 1.0,1.0 0.0,0.0 0.0)),"
                    "((0.0 0.0,1.0 0.0,1.0 0.0,0.0 0.0,0.0 0.0)),"
                    "((1.0 1.0,1.0 1.0,1.0 0.0,1.0 0.0,1.0 1.0)),"
                    "((0.0 1.0,0.0 1.0,1.0 1.0,1.0 1.0,0.0 1.0)),"
                    "((0.0 0.0,1.0 0.0,1.0 1.0,0.0 1.0,0.0 0.0))))");

  BOOST_CHECK(!solid->dropZ());
}

BOOST_AUTO_TEST_CASE(solidSwapXYTest)
{
  std::unique_ptr<Solid> emptySolid = std::make_unique<Solid>();
  BOOST_CHECK(emptySolid->isEmpty());
  BOOST_CHECK(!emptySolid->dropZ());

  std::string const         polyhedral1Str = "POLYHEDRALSURFACE ("
                                             "((0 0 0, 0 0 1, 0 1 1, 0 1 0, 0 0 0)),"
                                             "((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)),"
                                             "((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)),"
                                             "((1 1 0, 1 1 1, 1 0 1, 1 0 0, 1 1 0)),"
                                             "((0 1 0, 0 1 1, 1 1 1, 1 1 0, 0 1 0)),"
                                             "((0 0 1, 1 0 1, 1 1 1, 0 1 1, 0 0 1))"
                                             ")";
  std::unique_ptr<Geometry> shell1(io::readWkt(polyhedral1Str));

  std::unique_ptr<Solid> solid = std::make_unique<Solid>();
  solid->setExteriorShell(
      dynamic_cast<PolyhedralSurface *>(shell1.get())->clone());
  solid->swapXY();

  BOOST_CHECK_EQUAL(
      solid->asText(1),
      "SOLID Z (("
      "((0.0 0.0 0.0,0.0 0.0 1.0,1.0 0.0 1.0,1.0 0.0 0.0,0.0 0.0 0.0)),"
      "((0.0 0.0 0.0,1.0 0.0 0.0,1.0 1.0 0.0,0.0 1.0 0.0,0.0 0.0 0.0)),"
      "((0.0 0.0 0.0,0.0 1.0 0.0,0.0 1.0 1.0,0.0 0.0 1.0,0.0 0.0 0.0)),"
      "((1.0 1.0 0.0,1.0 1.0 1.0,0.0 1.0 1.0,0.0 1.0 0.0,1.0 1.0 0.0)),"
      "((1.0 0.0 0.0,1.0 0.0 1.0,1.0 1.0 1.0,1.0 1.0 0.0,1.0 0.0 0.0)),"
      "((0.0 0.0 1.0,0.0 1.0 1.0,1.0 1.0 1.0,1.0 0.0 1.0,0.0 0.0 1.0))))");
}

BOOST_AUTO_TEST_CASE(getCoordinateType)
{
  BOOST_CHECK_EQUAL(
      io::readWkt("SOLID((((0 0, 1 0, 0 1, 0 0))))")->getCoordinateType(),
      CoordinateType::COORDINATE_XY);
  BOOST_CHECK_EQUAL(io::readWkt("SOLID Z((((0 0 1, 1 0 1, 0 1 1, 0 0 1))))")
                        ->getCoordinateType(),
                    CoordinateType::COORDINATE_XYZ);
  BOOST_CHECK_EQUAL(io::readWkt("SOLID M((((0 0 1, 1 0 1, 0 1 1, 0 0 1))))")
                        ->getCoordinateType(),
                    CoordinateType::COORDINATE_XYM);
  BOOST_CHECK_EQUAL(
      io::readWkt("SOLID ZM((((0 0 1 2, 1 0 1 2, 0 1 1 2, 0 0 1 2))))")
          ->getCoordinateType(),
      CoordinateType::COORDINATE_XYZM);
}

BOOST_AUTO_TEST_SUITE_END()
