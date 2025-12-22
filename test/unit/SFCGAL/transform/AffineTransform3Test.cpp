// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/Kernel.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/detail/transform/AffineTransform3.h"
#include "SFCGAL/io/wkt.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_transform_AffineTransform3Test)

// Test for empty polygon - to verify AffineTransform3::transform(Polygon&)
// isEmpty() guard
BOOST_AUTO_TEST_CASE(testEmptyPolygon)
{
  std::unique_ptr<Geometry> g(io::readWkt("POLYGON EMPTY"));
  BOOST_CHECK(g->isEmpty());

  CGAL::Aff_transformation_3<Kernel> const affine(
      CGAL::TRANSLATION, CGAL::Vector_3<Kernel>(1.0, 2.0, 3.0));

  transform::AffineTransform3 transform(affine);
  g->accept(transform); // Should not crash

  BOOST_CHECK(g->isEmpty());
  BOOST_CHECK_EQUAL(g->asText(), "POLYGON EMPTY");
}

// Test for empty solid - to verify AffineTransform3::transform(Solid&)
// isEmpty() guard
BOOST_AUTO_TEST_CASE(testEmptySolid)
{
  std::unique_ptr<Geometry> g(io::readWkt("SOLID EMPTY"));
  BOOST_CHECK(g->isEmpty());

  CGAL::Aff_transformation_3<Kernel> const affine(
      CGAL::TRANSLATION, CGAL::Vector_3<Kernel>(1.0, 2.0, 3.0));

  transform::AffineTransform3 transform(affine);
  g->accept(transform); // Should not crash

  BOOST_CHECK(g->isEmpty());
  BOOST_CHECK_EQUAL(g->asText(), "SOLID EMPTY");
}

BOOST_AUTO_TEST_SUITE_END()
