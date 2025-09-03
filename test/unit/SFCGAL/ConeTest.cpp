// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/primitive3d/Cone.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/algorithm/area.h"
#include "SFCGAL/algorithm/volume.h"
#include <boost/test/unit_test.hpp>

using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(ConeTests)

BOOST_AUTO_TEST_CASE(testDefaultConstructor)
{
  Cone cone = Cone();
  BOOST_CHECK_EQUAL(cone.height(), 1.0);
  BOOST_CHECK_EQUAL(cone.bottomRadius(), 1.0);
  BOOST_CHECK_EQUAL(cone.topRadius(), 0.0);
  BOOST_CHECK_EQUAL(cone.numRadial(), 32);
}

BOOST_AUTO_TEST_CASE(testCustomConstructor)
{
  Cone cone(57, 13, 44, 100);
  BOOST_CHECK_EQUAL(cone.bottomRadius(), 57);
  BOOST_CHECK_EQUAL(cone.topRadius(), 13);
  BOOST_CHECK_EQUAL(cone.height(), 44);
  BOOST_CHECK_EQUAL(cone.numRadial(), 100);
}

BOOST_AUTO_TEST_CASE(testSetters)
{
  Cone cone;
  cone.setBottomRadius(15.0);
  cone.setNumRadial(12);

  BOOST_CHECK_CLOSE(cone.topRadius(), 0.0, 1e-6);
  BOOST_CHECK_CLOSE(cone.bottomRadius(), 15.0, 1e-6);
  BOOST_CHECK_EQUAL(cone.numRadial(), 12);

  // test setter and getter from Primitive
  cone.setParameter(std::string("bottom_radius"), 5.3);
  BOOST_CHECK_CLOSE(cone.bottomRadius(), 5.3, 1e-6);

  cone.setParameter(std::string("top_radius"), 2.3);
  BOOST_CHECK_CLOSE(cone.topRadius(), 2.3, 1e-6);

  // cone does not have a parameter called foo
  BOOST_CHECK_THROW(cone.setParameter(std::string("foo"), 12.3), Exception);
  BOOST_CHECK_THROW(static_cast<void>(cone.parameter("foo")), Exception);

  // bottom_radius is not a point
  BOOST_CHECK_THROW(
      cone.setParameter(std::string("bottom_radius"), Point_3(1, 1, 1)),
      Exception);

  // height, bottom and top radius cannot be negartive
  BOOST_CHECK_THROW(cone.setParameter(std::string("height"), -2.3), Exception);
  BOOST_CHECK_THROW(cone.setParameter(std::string("bottom_radius"), -20.38),
                    Exception);
  BOOST_CHECK_THROW(cone.setParameter(std::string("top_radius"), -5.2),
                    Exception);
}

BOOST_AUTO_TEST_CASE(testVolume)
{
  // small base, big top
  Cone   cone(23, 59, 68, 45);
  double expected_volume             = 382181.029;
  double expected_discretized_volume = 380940.437;
  BOOST_CHECK_CLOSE(cone.volume(), expected_volume, 1e-6);
  BOOST_CHECK_CLOSE(cone.volume(false), expected_volume, 1e-6);
  BOOST_CHECK_CLOSE(cone.volume(true), expected_discretized_volume, 1e-6);

  // big base, small top
  cone.setBottomRadius(59);
  cone.setTopRadius(23);
  BOOST_CHECK_CLOSE(cone.volume(), expected_volume, 1e-6);
  BOOST_CHECK_CLOSE(cone.volume(false), expected_volume, 1e-6);
  BOOST_CHECK_CLOSE(cone.volume(true), expected_discretized_volume, 1e-6);

  // same base and top (cylinder)
  cone.setBottomRadius(33);
  cone.setTopRadius(33);
  expected_volume             = 232641.219;
  expected_discretized_volume = 231886.045;
  BOOST_CHECK_CLOSE(cone.volume(), expected_volume, 1e-6);
  BOOST_CHECK_CLOSE(cone.volume(false), expected_volume, 1e-6);
  BOOST_CHECK_CLOSE(cone.volume(true), expected_discretized_volume, 1e-6);

  // top radius 0, real cone
  cone.setBottomRadius(36);
  cone.setTopRadius(0);
  expected_volume             = 92287.426;
  expected_discretized_volume = 91987.853;
  BOOST_CHECK_CLOSE(cone.volume(), expected_volume, 1e-6);
  BOOST_CHECK_CLOSE(cone.volume(false), expected_volume, 1e-6);
  BOOST_CHECK_CLOSE(cone.volume(true), expected_discretized_volume, 1e-6);

  // bottom radius 0, real cone upside-down
  cone.setBottomRadius(0);
  cone.setTopRadius(36);
  BOOST_CHECK_CLOSE(cone.volume(), expected_volume, 1e-6);
  BOOST_CHECK_CLOSE(cone.volume(false), expected_volume, 1e-6);
  BOOST_CHECK_CLOSE(cone.volume(true), expected_discretized_volume, 1e-6);

  // all radiuses 0
  cone.setBottomRadius(0);
  cone.setTopRadius(0);
  expected_volume             = 0;
  expected_discretized_volume = 0;
  BOOST_CHECK_EQUAL(cone.volume(), expected_volume);
  BOOST_CHECK_EQUAL(cone.volume(false), expected_volume);
  BOOST_CHECK_EQUAL(cone.volume(true), expected_discretized_volume);
}

BOOST_AUTO_TEST_CASE(testSurfaceArea)
{
  // small base, big top
  Cone   cone(23, 59, 68, 45);
  double expected_area             = 32418.741;
  double expected_discretized_area = 32351.199;
  BOOST_CHECK_CLOSE(cone.area3D(), expected_area, 1e-5);
  BOOST_CHECK_CLOSE(cone.area3D(false), expected_area, 1e-5);
  BOOST_CHECK_CLOSE(cone.area3D(true), expected_discretized_area, 1e-5);

  // big base, small top
  cone.setBottomRadius(59);
  cone.setTopRadius(23);
  BOOST_CHECK_CLOSE(cone.area3D(), expected_area, 1e-5);
  BOOST_CHECK_CLOSE(cone.area3D(false), expected_area, 1e-5);
  BOOST_CHECK_CLOSE(cone.area3D(true), expected_discretized_area, 1e-5);

  // same base and top (cylinder)
  cone.setBottomRadius(32);
  expected_area             = 16730.913;
  expected_discretized_area = 16704.954;
  BOOST_CHECK_CLOSE(cone.area3D(), expected_area, 1e-5);
  BOOST_CHECK_CLOSE(cone.area3D(false), expected_area, 1e-5);
  BOOST_CHECK_CLOSE(cone.area3D(true), expected_discretized_area, 1e-5);

  // top radius 0, real cone
  cone.setBottomRadius(87);
  cone.setTopRadius(0);
  cone.setHeight(13);
  cone.setNumRadial(30);
  expected_area             = 47821.428;
  expected_discretized_area = 47475.459;
  BOOST_CHECK_CLOSE(cone.area3D(), expected_area, 1e-5);
  BOOST_CHECK_CLOSE(cone.area3D(false), expected_area, 1e-5);
  BOOST_CHECK_CLOSE(cone.area3D(true), expected_discretized_area, 1e-5);

  // bottom radius 0, real cone upside-down
  cone.setTopRadius(32);
  cone.setBottomRadius(0);
  expected_area             = 6689.313;
  expected_discretized_area = 6643.212;
  BOOST_CHECK_CLOSE(cone.area3D(), expected_area, 1e-5);
  BOOST_CHECK_CLOSE(cone.area3D(false), expected_area, 1e-5);
  BOOST_CHECK_CLOSE(cone.area3D(true), expected_discretized_area, 1e-5);

  // all radiuses 0
  cone.setBottomRadius(0);
  cone.setTopRadius(0);
  expected_area             = 0;
  expected_discretized_area = 0;
  BOOST_CHECK_EQUAL(cone.area3D(), expected_area);
  BOOST_CHECK_EQUAL(cone.area3D(false), expected_area);
  BOOST_CHECK_EQUAL(cone.area3D(true), expected_discretized_area);
}

BOOST_AUTO_TEST_CASE(testPolyhedralSurface)
{
  Cone cone(54, 21, 87, 7);
  auto polyhedral_surface = cone.generatePolyhedralSurface();

  // Create PolyhedralSurface from Polyhedron and check WKT output
  BOOST_CHECK_EQUAL(
      polyhedral_surface.asText(1),
      "POLYHEDRALSURFACE Z (((54.0 0.0 0.0,33.7 -42.2 0.0,-12.0 -52.6 "
      "0.0,-48.7 -23.4 0.0,-48.7 23.4 0.0,-12.0 52.6 0.0,33.7 42.2 0.0,54.0 "
      "0.0 0.0)),((21.0 0.0 87.0,13.1 16.4 87.0,-4.7 20.5 87.0,-18.9 9.1 "
      "87.0,-18.9 -9.1 87.0,-4.7 -20.5 87.0,13.1 -16.4 87.0,21.0 0.0 "
      "87.0)),((54.0 0.0 0.0,21.0 0.0 87.0,13.1 -16.4 87.0,33.7 -42.2 0.0,54.0 "
      "0.0 0.0)),((33.7 -42.2 0.0,13.1 -16.4 87.0,-4.7 -20.5 87.0,-12.0 -52.6 "
      "0.0,33.7 -42.2 0.0)),((-12.0 -52.6 0.0,-4.7 -20.5 87.0,-18.9 -9.1 "
      "87.0,-48.7 -23.4 0.0,-12.0 -52.6 0.0)),((-48.7 -23.4 0.0,-18.9 -9.1 "
      "87.0,-18.9 9.1 87.0,-48.7 23.4 0.0,-48.7 -23.4 0.0)),((-48.7 23.4 "
      "0.0,-18.9 9.1 87.0,-4.7 20.5 87.0,-12.0 52.6 0.0,-48.7 23.4 "
      "0.0)),((-12.0 52.6 0.0,-4.7 20.5 87.0,13.1 16.4 87.0,33.7 42.2 "
      "0.0,-12.0 52.6 0.0)),((33.7 42.2 0.0,13.1 16.4 87.0,21.0 0.0 87.0,54.0 "
      "0.0 0.0,33.7 42.2 0.0)))");

  BOOST_CHECK_CLOSE(CGAL::to_double(algorithm::volume(
                        Solid(cone.generatePolyhedralSurface()))),
                    cone.volume(true), 1e-6);

  BOOST_CHECK_CLOSE(
      CGAL::to_double(algorithm::area3D(cone.generatePolyhedralSurface())),
      cone.area3D(true), 1e-6);
}

BOOST_AUTO_TEST_CASE(testGetSetHeight)
{
  Cone cone(6, 7, 8, 9);
  BOOST_CHECK_EQUAL(cone.height(), 8);

  cone.setHeight(6.904);
  BOOST_CHECK_EQUAL(cone.height(), 6.904);
}

BOOST_AUTO_TEST_CASE(testGetSetBottomRadius)
{
  Cone cone(6, 7, 8, 9);
  BOOST_CHECK_EQUAL(cone.bottomRadius(), 6);

  cone.setBottomRadius(5.904);
  BOOST_CHECK_EQUAL(cone.bottomRadius(), 5.904);
}

BOOST_AUTO_TEST_CASE(testGetSetTopRadius)
{
  Cone cone(6, 7, 8, 9);
  BOOST_CHECK_EQUAL(cone.topRadius(), 7);

  cone.setTopRadius(6.904);
  BOOST_CHECK_EQUAL(cone.topRadius(), 6.904);
}

BOOST_AUTO_TEST_CASE(testGetSetNumRadial)
{
  Cone cone(6, 7, 8, 9);
  BOOST_CHECK_EQUAL(cone.numRadial(), 9);

  cone.setNumRadial(54);
  BOOST_CHECK_EQUAL(cone.numRadial(), 54);
}

BOOST_AUTO_TEST_SUITE_END()
