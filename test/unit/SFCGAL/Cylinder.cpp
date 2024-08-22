#include "SFCGAL/Cylinder.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Sphere.h"
#include "SFCGAL/algorithm/area.h"
#include "SFCGAL/algorithm/volume.h"
#include "SFCGAL/io/OBJ.h"
#include "SFCGAL/io/wkt.h"
#include <boost/test/unit_test.hpp>
#include <cmath>

using namespace SFCGAL;
using Point_3  = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;

BOOST_AUTO_TEST_SUITE(CylinderTests)

BOOST_AUTO_TEST_CASE(testDefaultConstructor)
{
  Cylinder cyl;
  BOOST_CHECK_CLOSE(cyl.radius(), 1.0, 1e-6);
  BOOST_CHECK_CLOSE(cyl.height(), 1.0, 1e-6);
  BOOST_CHECK_EQUAL(cyl.numRadial(), 32);
  BOOST_CHECK_EQUAL(cyl.baseCenter(), Point_3(0, 0, 0));
  BOOST_CHECK_EQUAL(cyl.axis(), Vector_3(0, 0, 1));
}

BOOST_AUTO_TEST_CASE(testCustomConstructor)
{
  Point_3  base(1, 2, 3);
  Vector_3 axis(0, 1, 0);
  Cylinder cyl(base, axis, 2.0, 5.0, 16);
  BOOST_CHECK_CLOSE(cyl.radius(), 2.0, 1e-6);
  BOOST_CHECK_CLOSE(cyl.height(), 5.0, 1e-6);
  BOOST_CHECK_EQUAL(cyl.numRadial(), 16);
  BOOST_CHECK_EQUAL(cyl.baseCenter(), base);
  BOOST_CHECK_EQUAL(cyl.axis(), axis);
}

BOOST_AUTO_TEST_CASE(testSetters)
{
  Cylinder cyl;
  cyl.setRadius(3.0);
  cyl.setHeight(4.0);
  cyl.setNumRadial(24);
  cyl.setBaseCenter(Point_3(1, 1, 1));
  cyl.setAxis(Vector_3(1, 1, 1));

  BOOST_CHECK_CLOSE(cyl.radius(), 3.0, 1e-6);
  BOOST_CHECK_CLOSE(cyl.height(), 4.0, 1e-6);
  BOOST_CHECK_EQUAL(cyl.numRadial(), 24);
  BOOST_CHECK_EQUAL(cyl.baseCenter(), Point_3(1, 1, 1));
  BOOST_CHECK_EQUAL(cyl.axis(), Vector_3(1, 1, 1));
}

BOOST_AUTO_TEST_CASE(testGenerateSurfaceMesh)
{
  Cylinder cyl(Point_3(0, 0, 0), Vector_3(0, 0, 1), 1.0, 2.0, 4);
  auto     mesh = cyl.generateSurfaceMesh();

  BOOST_CHECK_EQUAL(mesh.number_of_vertices(), cyl.numRadial() * 2 + 2);
  BOOST_CHECK_EQUAL(mesh.number_of_edges(), cyl.numRadial() * 6);
  BOOST_CHECK_EQUAL(mesh.number_of_faces(), cyl.numRadial() * 4);
}

BOOST_AUTO_TEST_CASE(testVolume)
{
  Cylinder cyl(Point_3(0, 0, 0), Vector_3(0, 0, 1), 2.0, 5.0, 32);
  double   volume          = cyl.volume();
  double   expected_volume = M_PI * 2.0 * 2.0 * 5.0;
  BOOST_CHECK_CLOSE(volume, expected_volume, 0.01);
}

BOOST_AUTO_TEST_CASE(testSurfaceArea)
{
  Cylinder cyl(Point_3(0, 0, 0), Vector_3(0, 0, 1), 2.0, 5.0, 32);
  double   area          = cyl.area();
  double   expected_area = 2 * M_PI * 2.0 * 2.0 + 2 * M_PI * 2.0 * 5.0;
  BOOST_CHECK_CLOSE(area, expected_area, 0.01);
}

BOOST_AUTO_TEST_CASE(testTiltedCylinder)
{
  Point_3  base(1, 1, 1);
  Vector_3 axis(1, 1, 1);
  Cylinder cyl(base, axis, 1.0, std::sqrt(3.0), 16);
  auto     mesh = cyl.generateSurfaceMesh();

  // Check that the top center is where we expect it to be
  Point_3 expected_top(2, 2, 2);
  bool    found_top = false;
  for (auto v : mesh.vertices()) {
    if (CGAL::squared_distance(mesh.point(v), expected_top) < 1e-10) {
      found_top = true;
      break;
    }
  }
  BOOST_CHECK(found_top);

  // Create PolyhedralSurface and output WKT for visual inspection
  PolyhedralSurface poly_surface(mesh);
  BOOST_TEST_MESSAGE("Tilted Cylinder WKT: " << poly_surface.asText(1));
  std::cout << "Tilted Cylinder WKT: " << poly_surface.asText(1) << std::endl;
  SFCGAL::io::OBJ::save(poly_surface, "/tmp/surface.obj");
}

BOOST_AUTO_TEST_CASE(testPolyhedron)
{
  Cylinder cyl(Point_3(0, 0, 0), Vector_3(0, 0, 1), 1.0, 2.0, 8);
  auto     polyhedron = cyl.generatePolyhedron();

  BOOST_CHECK_EQUAL(polyhedron.size_of_vertices(), cyl.numRadial() * 2 + 2);
  BOOST_CHECK_EQUAL(polyhedron.size_of_facets(), cyl.numRadial() * 4);

  // Create PolyhedralSurface from Polyhedron and check WKT output
  PolyhedralSurface poly_surface(polyhedron);
  std::cout << "Generated Polyhedron WKT: " << poly_surface.asText(1)
            << std::endl;

  SFCGAL::io::OBJ::save(poly_surface, "/tmp/polyhedron.obj");
}

BOOST_AUTO_TEST_CASE(testBuffer3D)
{

  std::vector<Point_3> points = {Point_3(-100, 0, 0), Point_3(40, -70, 0),
                                 Point_3(40, 50, 40), Point_3(-90, -60, 60),
                                 Point_3(0, 0, -100), Point_3(30, 0, 150)};
  Kernel::FT           radius(10.0);
  int                  segments = 16;

  std::vector<CGAL::Polyhedron_3<Kernel>> meshes;

  for (size_t i = 0; i < points.size() - 1; ++i) {
    // Create and add cylinder
    Vector_3   axis   = points[i + 1] - points[i];
    Kernel::FT height = CGAL::sqrt(CGAL::to_double(axis.squared_length()));
    Cylinder   cyl(points[i], axis, radius, height, segments);
    meshes.push_back(cyl.generatePolyhedron());

    // Add sphere at joints
    if (i > 0) {
      Sphere sphere(radius, points[i], segments / 2, segments);
      meshes.push_back(sphere.generatePolyhedron());
    }
  }

  // Combine all meshes into a single PolyhedralSurface
  PolyhedralSurface poly_surface;
  for (const auto &mesh : meshes) {
    poly_surface.addPolygons(mesh);
  }

  // Output WKT
  BOOST_TEST_MESSAGE("Buffer 3D WKT: " << poly_surface.asText(1));
  std::cout << "Buffer 3D WKT: " << poly_surface.asText(1) << std::endl;

  // Save as OBJ
  std::string obj_filename = "/tmp/buffer_3d.obj";
  SFCGAL::io::OBJ::save(poly_surface, obj_filename);
  BOOST_TEST_MESSAGE("Saved Buffer 3D as OBJ: " << obj_filename);

  // Add some assertions
  BOOST_CHECK(poly_surface.numGeometries() > 0);
  BOOST_CHECK(poly_surface.is3D());

  // Check the number of geometries
  size_t expected_geometries =
      (points.size() - 1) * 2 -
      1; // cylinders + spheres (except for the first and last points)
  BOOST_CHECK_EQUAL(poly_surface.numGeometries(), expected_geometries);
}

BOOST_AUTO_TEST_SUITE_END()
