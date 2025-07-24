// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/primitive3d/Torus.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/algorithm/area.h"
#include "SFCGAL/algorithm/volume.h"
#include "SFCGAL/io/OBJ.h"
#include <boost/test/unit_test.hpp>

using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(TorusTests)

BOOST_AUTO_TEST_CASE(testDefaultConstructor)
{
  Torus torus = Torus();
  BOOST_CHECK_EQUAL(torus.mainRadius(), 10.0);
  BOOST_CHECK_EQUAL(torus.tubeRadius(), 2.0);
  BOOST_CHECK_EQUAL(torus.mainNumRadial(), 32);
  BOOST_CHECK_EQUAL(torus.tubeNumRadial(), 16);

  torus.setMainNumRadial(3);
  torus.setTubeNumRadial(3);

  BOOST_CHECK_EQUAL(torus.mainNumRadial(), 3);
  BOOST_CHECK_EQUAL(torus.tubeNumRadial(), 3);
}

BOOST_AUTO_TEST_CASE(testCustomConstructor)
{
  Torus torus(57, 13, 44, 100);
  BOOST_CHECK_EQUAL(torus.mainRadius(), 57);
  BOOST_CHECK_EQUAL(torus.tubeRadius(), 13);
  BOOST_CHECK_EQUAL(torus.mainNumRadial(), 44);
  BOOST_CHECK_EQUAL(torus.tubeNumRadial(), 100);
}

BOOST_AUTO_TEST_CASE(testVolume)
{
  Torus  torus(59, 23, 39, 59);
  double expected_volume = 616080.44;
  BOOST_CHECK_CLOSE(torus.volume(), expected_volume, 1e-6);
}

BOOST_AUTO_TEST_CASE(testSurfaceArea)
{
  Torus  torus(59, 23, 34, 123);
  double expected_area = 53572.21;
  BOOST_CHECK_CLOSE(torus.area3D(), expected_area, 1e-5);
}

BOOST_AUTO_TEST_CASE(testPolyhedralSurface)
{
  Torus torus(10, 2, 4, 4);
  auto  polyhedral_surface = torus.generatePolyhedralSurface();

  // Create PolyhedralSurface from Polyhedron and check WKT output
  BOOST_CHECK_EQUAL(
      polyhedral_surface.asText(1),
      "POLYHEDRALSURFACE Z (((12.0 0.0 0.0,0.0 12.0 0.0,0.0 10.0 2.0,10.0 0.0 "
      "2.0,12.0 0.0 0.0)),((10.0 0.0 2.0,0.0 10.0 2.0,0.0 8.0 0.0,8.0 0.0 "
      "0.0,10.0 0.0 2.0)),((8.0 0.0 0.0,0.0 8.0 0.0,0.0 10.0 -2.0,10.0 0.0 "
      "-2.0,8.0 0.0 0.0)),((10.0 0.0 -2.0,0.0 10.0 -2.0,0.0 12.0 0.0,12.0 0.0 "
      "0.0,10.0 0.0 -2.0)),((0.0 12.0 0.0,-12.0 0.0 0.0,-10.0 0.0 2.0,0.0 10.0 "
      "2.0,0.0 12.0 0.0)),((0.0 10.0 2.0,-10.0 0.0 2.0,-8.0 0.0 0.0,0.0 8.0 "
      "0.0,0.0 10.0 2.0)),((0.0 8.0 0.0,-8.0 0.0 0.0,-10.0 0.0 -2.0,0.0 10.0 "
      "-2.0,0.0 8.0 0.0)),((0.0 10.0 -2.0,-10.0 0.0 -2.0,-12.0 0.0 0.0,0.0 "
      "12.0 0.0,0.0 10.0 -2.0)),((-12.0 0.0 0.0,0.0 -12.0 0.0,0.0 -10.0 "
      "2.0,-10.0 0.0 2.0,-12.0 0.0 0.0)),((-10.0 0.0 2.0,0.0 -10.0 2.0,0.0 "
      "-8.0 0.0,-8.0 0.0 0.0,-10.0 0.0 2.0)),((-8.0 0.0 0.0,0.0 -8.0 0.0,0.0 "
      "-10.0 -2.0,-10.0 0.0 -2.0,-8.0 0.0 0.0)),((-10.0 0.0 -2.0,0.0 -10.0 "
      "-2.0,0.0 -12.0 0.0,-12.0 0.0 0.0,-10.0 0.0 -2.0)),((0.0 -12.0 0.0,12.0 "
      "0.0 0.0,10.0 0.0 2.0,0.0 -10.0 2.0,0.0 -12.0 0.0)),((0.0 -10.0 2.0,10.0 "
      "0.0 2.0,8.0 0.0 0.0,0.0 -8.0 0.0,0.0 -10.0 2.0)),((0.0 -8.0 0.0,8.0 0.0 "
      "0.0,10.0 0.0 -2.0,0.0 -10.0 -2.0,0.0 -8.0 0.0)),((0.0 -10.0 -2.0,10.0 "
      "0.0 -2.0,12.0 0.0 0.0,0.0 -12.0 0.0,0.0 -10.0 -2.0)))");
}

BOOST_AUTO_TEST_CASE(testGetSetMainRadius)
{
  Torus torus(6, 3, 8, 9);
  BOOST_CHECK_EQUAL(torus.mainRadius(), 6);

  torus.setMainRadius(5.904);
  BOOST_CHECK_EQUAL(torus.mainRadius(), 5.904);

  // main radius cannot be negative
  BOOST_CHECK_THROW(torus.setMainRadius(-1.2), Exception);

  // main radius cannot be smaller than tube radius
  BOOST_CHECK_THROW(torus.setMainRadius(2.2), Exception);
}

BOOST_AUTO_TEST_CASE(testGetSetTubeRadius)
{
  Torus torus(9, 7, 8, 9);
  BOOST_CHECK_EQUAL(torus.tubeRadius(), 7);

  torus.setTubeRadius(6.904);
  BOOST_CHECK_EQUAL(torus.tubeRadius(), 6.904);

  // tube radius cannot be negative
  BOOST_CHECK_THROW(torus.setTubeRadius(-1.2), Exception);

  // tube radius cannot be greater than main radius
  BOOST_CHECK_THROW(torus.setTubeRadius(32.2), Exception);
}

BOOST_AUTO_TEST_CASE(testGetSetMainNumRadial)
{
  Torus torus(12, 7, 8, 9);
  BOOST_CHECK_EQUAL(torus.mainNumRadial(), 8);

  torus.setMainNumRadial(54);
  BOOST_CHECK_EQUAL(torus.mainNumRadial(), 54);
}

BOOST_AUTO_TEST_CASE(testGetSetTubeNumRadial)
{
  Torus torus(10, 7, 8, 9);
  BOOST_CHECK_EQUAL(torus.tubeNumRadial(), 9);

  torus.setTubeNumRadial(54);
  BOOST_CHECK_EQUAL(torus.tubeNumRadial(), 54);
}

BOOST_AUTO_TEST_SUITE_END()
