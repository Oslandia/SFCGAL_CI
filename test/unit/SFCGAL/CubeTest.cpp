// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/primitive3d/Cube.h"
#include "SFCGAL/algorithm/isValid.h"
#include <boost/test/unit_test.hpp>

using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(CubeTests)

BOOST_AUTO_TEST_CASE(testDefaultConstructor)
{
  Cube cube;
  BOOST_CHECK_EQUAL(cube.size(), 1.0);
}

BOOST_AUTO_TEST_CASE(testCustomConstructor)
{
  Cube cube(57);
  BOOST_CHECK_EQUAL(cube.size(), 57);
}

BOOST_AUTO_TEST_CASE(testSetters)
{
  Cube cube;
  cube.setSize(3.0);

  BOOST_CHECK_CLOSE(cube.size(), 3.0, 1e-6);

  // size cannot be negative
  BOOST_CHECK_THROW(cube.setSize(-1.2), Exception);
}

BOOST_AUTO_TEST_CASE(testVolume)
{
  Cube   cube(23);
  double volume = cube.volume();
  BOOST_CHECK_EQUAL(volume, 12167);
}

BOOST_AUTO_TEST_CASE(testSurfaceArea)
{
  Cube   cube(14);
  double area = cube.area3D();
  BOOST_CHECK_EQUAL(area, 1176);
}

BOOST_AUTO_TEST_CASE(testPolyhedralSurface)
{
  // Create PolyhedralSurface from Polyhedron and check WKT output
  Cube cube(24);
  auto polyhedral_surface = cube.generatePolyhedralSurface();

  BOOST_CHECK(algorithm::isValid(polyhedral_surface));

  BOOST_CHECK_EQUAL(
      polyhedral_surface.asText(1),
      "POLYHEDRALSURFACE Z (((0.0 0.0 0.0,0.0 24.0 0.0,24.0 24.0 0.0,24.0 0.0 "
      "0.0,0.0 0.0 0.0)),((0.0 0.0 24.0,24.0 0.0 24.0,24.0 24.0 24.0,0.0 24.0 "
      "24.0,0.0 0.0 24.0)),((0.0 0.0 0.0,24.0 0.0 0.0,24.0 0.0 24.0,0.0 0.0 "
      "24.0,0.0 0.0 0.0)),((0.0 24.0 0.0,0.0 24.0 24.0,24.0 24.0 24.0,24.0 "
      "24.0 0.0,0.0 24.0 0.0)),((24.0 0.0 0.0,24.0 24.0 0.0,24.0 24.0 "
      "24.0,24.0 0.0 24.0,24.0 0.0 0.0)),((0.0 0.0 0.0,0.0 0.0 24.0,0.0 24.0 "
      "24.0,0.0 24.0 0.0,0.0 0.0 0.0)))");
}

BOOST_AUTO_TEST_CASE(testGetSetSize)
{
  Cube cube(3.5);
  BOOST_CHECK_EQUAL(cube.size(), 3.5);

  cube.setSize(6.904);
  BOOST_CHECK_EQUAL(cube.size(), 6.904);
}

BOOST_AUTO_TEST_SUITE_END()
