#include "SFCGAL/Sphere.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/io/wkt.h"
#include <boost/test/unit_test.hpp>
#include <cmath>

using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SphereTests)

BOOST_AUTO_TEST_CASE(testDefaultConstructor)
{
  Sphere sphere;
  BOOST_CHECK_CLOSE(sphere.radius(), 1.0, 1e-6);
  BOOST_CHECK_EQUAL(sphere.center(), Point_3(0, 0, 0));
  BOOST_CHECK_EQUAL(sphere.numVertical(), 16);
  BOOST_CHECK_EQUAL(sphere.numHorizontal(), 32);
}

BOOST_AUTO_TEST_CASE(testCustomConstructor)
{
  Point_3 center(1, 2, 3);
  Sphere  sphere(2.0, center, 20, 40);
  BOOST_CHECK_CLOSE(sphere.radius(), 2.0, 1e-6);
  BOOST_CHECK_EQUAL(sphere.center(), center);
  BOOST_CHECK_EQUAL(sphere.numVertical(), 20);
  BOOST_CHECK_EQUAL(sphere.numHorizontal(), 40);
}

BOOST_AUTO_TEST_CASE(testSetters)
{
  Sphere sphere;
  sphere.setRadius(3.0);
  sphere.setCenter(Point_3(1, 1, 1));
  sphere.setNumVertical(24);
  sphere.setNumHorizontal(48);

  BOOST_CHECK_CLOSE(sphere.radius(), 3.0, 1e-6);
  BOOST_CHECK_EQUAL(sphere.center(), Point_3(1, 1, 1));
  BOOST_CHECK_EQUAL(sphere.numVertical(), 24);
  BOOST_CHECK_EQUAL(sphere.numHorizontal(), 48);
}

BOOST_AUTO_TEST_CASE(testGeneratePolyhedron)
{
  Sphere sphere(1.0, Point_3(0, 0, 0), 8, 16);
  auto   polyhedron = sphere.generatePolyhedron();

  // The number of vertices should be (numVertical - 1) * numHorizontal + 2
  BOOST_CHECK_EQUAL(polyhedron.size_of_vertices(), 7 * 16 + 2);

  // The number of faces should be (numVertical -1) * numHorizontal * 2
  BOOST_CHECK_EQUAL(polyhedron.size_of_facets(), 7 * 16 * 2);
}

BOOST_AUTO_TEST_CASE(testGeneratePoints)
{
  Sphere sphere(1.0, Point_3(0, 0, 0), 8, 16);
  auto   points = sphere.generatePoints();

  // The number of points should be numVertical * numHorizontal
  BOOST_CHECK_EQUAL(points.size(), 8 * 16);

  // Check that all points are at a distance of 1 from the center
  for (const auto &point : points) {
    BOOST_CHECK_CLOSE(CGAL::sqrt(CGAL::to_double(
                          CGAL::squared_distance(Point_3(0, 0, 0), point))),
                      1.0, 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(testVolume)
{
  Sphere sphere(2.0, Point_3(0, 0, 0), 32, 64);
  double volume          = sphere.volume();
  double expected_volume = 4.0 / 3.0 * M_PI * 2.0 * 2.0 * 2.0;
  BOOST_CHECK_CLOSE(volume, expected_volume, 0.1); // 0.1% tolerance
}

BOOST_AUTO_TEST_CASE(testSurfaceArea)
{
  Sphere sphere(2.0, Point_3(0, 0, 0), 32, 64);
  double area          = sphere.area();
  double expected_area = 4.0 * M_PI * 2.0 * 2.0;
  BOOST_CHECK_CLOSE(area, expected_area, 0.1); // 0.1% tolerance
}

BOOST_AUTO_TEST_CASE(testWKT)
{
  Sphere            sphere(1.0, Point_3(0, 0, 0), 4, 8);
  auto              polyhedron = sphere.generatePolyhedron();
  PolyhedralSurface surface(polyhedron);
  std::string       wkt = surface.asText(1);

  BOOST_TEST_MESSAGE("Sphere WKT: " << wkt);
  BOOST_CHECK(wkt.find("POLYHEDRALSURFACE Z") == 0);
}

BOOST_AUTO_TEST_SUITE_END()
