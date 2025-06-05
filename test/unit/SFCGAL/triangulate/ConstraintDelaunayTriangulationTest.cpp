// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/Exception.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/detail/triangulate/ConstraintDelaunayTriangulation.h"

using namespace boost::unit_test;
using namespace SFCGAL;
using namespace SFCGAL::triangulate;

BOOST_AUTO_TEST_SUITE(SFCGAL_triangulate_ConstraintDelaunayTriangulationTest)

/// Coordinate() ;
BOOST_AUTO_TEST_CASE(testDefaultConstructor)
{
  ConstraintDelaunayTriangulation const triangulation;
  BOOST_CHECK_EQUAL(triangulation.numVertices(), 0U);
  BOOST_CHECK_EQUAL(triangulation.numTriangles(), 0U);
}

BOOST_AUTO_TEST_CASE(testTriangulateSquare)
{
  ConstraintDelaunayTriangulation                        triangulation;
  typedef ConstraintDelaunayTriangulation::Vertex_handle Vertex_handle;
  typedef ConstraintDelaunayTriangulation::Face_handle   Face_handle;
  typedef ConstraintDelaunayTriangulation::All_faces_iterator
      All_faces_iterator;
  // typedef ConstraintDelaunayTriangulation::Finite_faces_iterator
  // Finite_faces_iterator ;

  Vertex_handle const a = triangulation.addVertex(Coordinate(0.0, 0.0));
  Vertex_handle const b = triangulation.addVertex(Coordinate(1.0, 0.0));
  Vertex_handle const c = triangulation.addVertex(Coordinate(1.0, 1.0));
  Vertex_handle const d = triangulation.addVertex(Coordinate(0.0, 1.0));

  BOOST_CHECK_EQUAL(triangulation.numVertices(), 4U);
  BOOST_CHECK_EQUAL(triangulation.numTriangles(), 2U);

  triangulation.addConstraint(a, b);
  triangulation.addConstraint(b, c);
  triangulation.addConstraint(c, d);
  triangulation.addConstraint(d, a);

  // constraint have no impact
  BOOST_CHECK_EQUAL(triangulation.numVertices(), 4U);
  BOOST_CHECK_EQUAL(triangulation.numTriangles(), 2U);

  /*
   * nesting level
   */
  for (All_faces_iterator it = triangulation.all_faces_begin();
       it != triangulation.all_faces_end(); ++it) {
    BOOST_CHECK_EQUAL(it->info().nestingLevel, -1);
  }

  // check mark domains
  triangulation.markDomains();

  for (All_faces_iterator it = triangulation.all_faces_begin();
       it != triangulation.all_faces_end(); ++it) {
    Face_handle const face = it;

    if (triangulation.isInfinite(face)) {
      BOOST_CHECK_EQUAL(it->info().nestingLevel, 0);
    } else {
      BOOST_CHECK_EQUAL(it->info().nestingLevel, 1);
    }
  }
}

BOOST_AUTO_TEST_CASE(testProjectionPlane)
{
  ConstraintDelaunayTriangulation triangulation;
  // typedef ConstraintDelaunayTriangulation::Vertex_handle Vertex_handle ;
  // typedef ConstraintDelaunayTriangulation::Face_handle           Face_handle
  // ;

  triangulation.setProjectionPlane(Kernel::Plane_3(
      Kernel::RT(1), Kernel::RT(0), Kernel::RT(0), Kernel::RT(0)));

  triangulation.addVertex(Coordinate(1.0, 0.0, 0.0));
  triangulation.addVertex(Coordinate(1.0, 1.0, 0.0));
  triangulation.addVertex(Coordinate(1.0, 1.0, 1.0));
  triangulation.addVertex(Coordinate(1.0, 0.0, 1.0));

  BOOST_CHECK_EQUAL(triangulation.numVertices(), 4U);
  BOOST_CHECK_EQUAL(triangulation.numTriangles(), 2U);
}

BOOST_AUTO_TEST_SUITE_END()
