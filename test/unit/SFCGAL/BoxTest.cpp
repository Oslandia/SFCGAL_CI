// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/primitive3d/Box.h"
#include "SFCGAL/algorithm/isValid.h"
#include <boost/test/unit_test.hpp>

using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(BoxTests)

BOOST_AUTO_TEST_CASE(testDefaultConstructor)
{
  Box box;
  BOOST_CHECK_EQUAL(box.xExtent(), 1.0);
  BOOST_CHECK_EQUAL(box.yExtent(), 1.0);
  BOOST_CHECK_EQUAL(box.zExtent(), 1.0);
}

BOOST_AUTO_TEST_CASE(testCustomConstructor)
{
  Box box(14, 10, 57);
  BOOST_CHECK_EQUAL(box.xExtent(), 14);
  BOOST_CHECK_EQUAL(box.yExtent(), 10);
  BOOST_CHECK_EQUAL(box.zExtent(), 57);
}

BOOST_AUTO_TEST_CASE(testSetters)
{
  Box box;
  box.setXExtent(3.0);
  box.setYExtent(5.0);
  box.setZExtent(7.0);

  BOOST_CHECK_CLOSE(box.xExtent(), 3.0, 1e-6);
  BOOST_CHECK_CLOSE(box.yExtent(), 5.0, 1e-6);
  BOOST_CHECK_CLOSE(box.zExtent(), 7.0, 1e-6);

  // extent cannot be negative
  BOOST_CHECK_THROW(box.setXExtent(-1.2), Exception);
  BOOST_CHECK_THROW(box.setYExtent(-2.2), Exception);
  BOOST_CHECK_THROW(box.setZExtent(-5.4), Exception);
}

BOOST_AUTO_TEST_CASE(testVolume)
{
  Box    box(23, 9, 91);
  double volume = box.volume();
  BOOST_CHECK_EQUAL(volume, 18837);
}

BOOST_AUTO_TEST_CASE(testSurfaceArea)
{
  Box    box(14, 4, 96);
  double area = box.area3D();
  BOOST_CHECK_EQUAL(area, 3568);
}

BOOST_AUTO_TEST_CASE(testPolyhedralSurface)
{
  // Create PolyhedralSurface from Polyhedron and check WKT output
  Box  box(24, 10, 82);
  auto polyhedral_surface = box.generatePolyhedralSurface();

  BOOST_CHECK(algorithm::isValid(polyhedral_surface));

  BOOST_CHECK_EQUAL(
      polyhedral_surface.asText(1),
      "POLYHEDRALSURFACE Z (((0.0 0.0 0.0,0.0 10.0 0.0,24.0 10.0 0.0,24.0 0.0 "
      "0.0,0.0 0.0 0.0)),((0.0 0.0 82.0,24.0 0.0 82.0,24.0 10.0 82.0,0.0 10.0 "
      "82.0,0.0 0.0 82.0)),((0.0 0.0 0.0,24.0 0.0 0.0,24.0 0.0 82.0,0.0 0.0 "
      "82.0,0.0 0.0 0.0)),((0.0 10.0 0.0,0.0 10.0 82.0,24.0 10.0 82.0,24.0 "
      "10.0 0.0,0.0 10.0 0.0)),((24.0 0.0 0.0,24.0 10.0 0.0,24.0 10.0 "
      "82.0,24.0 0.0 82.0,24.0 0.0 0.0)),((0.0 0.0 0.0,0.0 0.0 82.0,0.0 10.0 "
      "82.0,0.0 10.0 0.0,0.0 0.0 0.0)))");
}

BOOST_AUTO_TEST_CASE(testGetSetExtents)
{
  Box box(3.5, 62.78, 1.621);
  BOOST_CHECK_EQUAL(box.xExtent(), 3.5);
  BOOST_CHECK_EQUAL(box.yExtent(), 62.78);
  BOOST_CHECK_EQUAL(box.zExtent(), 1.621);

  box.setXExtent(6.904);
  box.setYExtent(9.44);
  box.setZExtent(543.8209);

  BOOST_CHECK_EQUAL(box.xExtent(), 6.904);
  BOOST_CHECK_EQUAL(box.yExtent(), 9.44);
  BOOST_CHECK_EQUAL(box.zExtent(), 543.8209);
}

BOOST_AUTO_TEST_SUITE_END()
