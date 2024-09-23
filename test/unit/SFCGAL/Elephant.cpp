#include <SFCGAL/Cylinder.h>
#include <SFCGAL/Kernel.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Sphere.h>
#include <SFCGAL/algorithm/boolean3D.h>
#include <SFCGAL/algorithm/buffer3D.h>
#include <SFCGAL/algorithm/rotate.h>
#include <SFCGAL/algorithm/scale.h>
#include <SFCGAL/algorithm/translate.h>
#include <SFCGAL/io/OBJ.h>
#include <SFCGAL/io/wkb.h>
#include <boost/test/unit_test.hpp>
#include <cmath>

using namespace SFCGAL;
using namespace SFCGAL::algorithm;
using Point_3          = Kernel::Point_3;
using Vector_3         = Kernel::Vector_3;
using Polyhedron_3     = CGAL::Polyhedron_3<Kernel>;
using Nef_Polyhedron_3 = CGAL::Nef_polyhedron_3<Kernel>;

struct ElephantConstants {
  static constexpr double SCALE_FACTOR = 1.0;

  // Body
  static constexpr double BODY_RADIUS = 22 * SCALE_FACTOR;

  // Head
  static constexpr double HEAD_RADIUS = 18 * SCALE_FACTOR;

  // Ear
  static constexpr double EAR_RADIUS = 10 * SCALE_FACTOR;

  // Eye
  static constexpr double EYE_RADIUS = 2 * SCALE_FACTOR;

  // Trunk
  static constexpr double TRUNK_RADIUS = 3 * SCALE_FACTOR;
  static constexpr int    TRUNK_POINTS = 30;
  static constexpr double GLOBE_RADIUS = 5 * SCALE_FACTOR;

  // Legs
  static constexpr double LEG_RADIUS = 5 * SCALE_FACTOR;

  // Tail
  static constexpr double TAIL_RADIUS = 1.5 * SCALE_FACTOR;
};

BOOST_AUTO_TEST_SUITE(ElephantTest)

std::unique_ptr<Geometry>
createSphere(double radius, double x, double y, double z, double scaleX = 1.0,
             double scaleY = 1.0, double scaleZ = 1.0)
{
  Sphere sphere(radius, Point_3(x, y, z));
  auto   generated = sphere.generatePolyhedron();
  auto   ph        = std::make_unique<PolyhedralSurface>(generated);
  algorithm::scale(*ph, scaleX, scaleY, scaleZ);
  return ph;
}

std::unique_ptr<Geometry>
createCylinder(double x1, double y1, double z1, double x2, double y2, double z2,
               double radius)
{
  Point_3  base_center(x1, y1, z1);
  Point_3  top_center(x2, y2, z2);
  Vector_3 axis   = top_center - base_center;
  double   height = std::sqrt(CGAL::to_double(axis.squared_length()));

  Cylinder cylinder(base_center, axis, radius, height);
  auto     generated = cylinder.generatePolyhedron();
  return std::make_unique<PolyhedralSurface>(generated);
}

std::unique_ptr<Geometry>
createEar(double rotateX, double rotateY, double rotateZ, double posX,
          double posY, double posZ)
{
  Sphere sphear(ElephantConstants::EAR_RADIUS, Point_3(0, 0, 0));
  auto   generated = sphear.generatePolyhedron();
  auto   ear       = std::make_unique<PolyhedralSurface>(generated);
  algorithm::scale(*ear, 1.0, 1.0, 0.2);
  algorithm::rotate(*ear, rotateX * M_PI / 180.0, Vector_3(1.0, 0.0, 0.0));
  algorithm::rotate(*ear, rotateY * M_PI / 180.0, Vector_3(0.0, 1.0, 0.0));
  algorithm::rotate(*ear, rotateZ * M_PI / 180.0, Vector_3(0.0, 0.0, 1.0));
  algorithm::translate(*ear, posX, posY, posZ);
  return ear;
}

std::unique_ptr<Geometry>
createTrunk()
{
  const double       x0     = 27;
  const double       y0     = 0;
  const double       z0     = 20;
  const double       length = 20;
  std::vector<Point> trunk_points;
  for (int i = 0; i <= ElephantConstants::TRUNK_POINTS; i++) {
    double t = static_cast<double>(i) / ElephantConstants::TRUNK_POINTS;
    double x = x0 + length * t;
    double z = z0 - 15 * std::sin(M_PI * t);
    trunk_points.emplace_back(Point(x, 0, z));
  }
  LineString trunk_line(trunk_points);
  auto       buf = SFCGAL::algorithm::Buffer3D(trunk_line,
                                               ElephantConstants::TRUNK_RADIUS, 16);
  return buf.compute(SFCGAL::algorithm::Buffer3D::BufferType::ROUND);
}

BOOST_AUTO_TEST_CASE(testElephant)
{
  // Body
  std::unique_ptr<Geometry> elephant =
      createSphere(ElephantConstants::BODY_RADIUS, 0, 0, 0, 1.1, 1.0, 0.9);

  // Head
  auto head =
      createSphere(ElephantConstants::HEAD_RADIUS, 20, 0, 15, 0.9, 0.9, 0.8);
  elephant = Boolean3D::union3D(*elephant, *head);

  // Ears
  auto ear_right = createEar(20, -10, -20, 15, 20, 20);
  elephant       = Boolean3D::union3D(*elephant, *ear_right);

  auto ear_left = createEar(-20, -10, 20, 15, -20, 20);
  elephant      = Boolean3D::union3D(*elephant, *ear_left);

  // Eyes
  auto eye_right = createSphere(ElephantConstants::EYE_RADIUS, 25, 8, 22);
  elephant       = Boolean3D::union3D(*elephant, *eye_right);

  auto eye_left = createSphere(ElephantConstants::EYE_RADIUS, 25, -8, 22);
  elephant      = Boolean3D::union3D(*elephant, *eye_left);

  // Trunk
  auto trunk = createTrunk();
  elephant   = Boolean3D::union3D(*elephant, *trunk);

  // Globe
  auto globe = createSphere(ElephantConstants::GLOBE_RADIUS, 47, 0, 25);
  elephant   = Boolean3D::union3D(*elephant, *globe);

  // Legs
  auto leg_fr =
      createCylinder(10, 12, -20, 10, 12, 0, ElephantConstants::LEG_RADIUS);
  elephant = Boolean3D::union3D(*elephant, *leg_fr);

  auto leg_fl =
      createCylinder(10, -12, -20, 10, -12, 0, ElephantConstants::LEG_RADIUS);
  elephant = Boolean3D::union3D(*elephant, *leg_fl);

  auto leg_br =
      createCylinder(-15, 10, -20, -15, 10, 0, ElephantConstants::LEG_RADIUS);
  elephant = Boolean3D::union3D(*elephant, *leg_br);

  auto leg_bl =
      createCylinder(-15, -10, -20, -15, -10, 0, ElephantConstants::LEG_RADIUS);
  elephant = Boolean3D::union3D(*elephant, *leg_bl);

  // Tail
  auto tail =
      createCylinder(-22, 0, 5, -37, 0, -10, ElephantConstants::TAIL_RADIUS);
  elephant = Boolean3D::union3D(*elephant, *tail);

  // Save to OBJ file
  io::OBJ::save(*elephant, "/tmp/elephant.obj");

  // Add a simple test to ensure the geometry is not empty
  BOOST_CHECK(!elephant->isEmpty());

  std::cout << elephant->asWkb(boost::endian::order::native, true) << std::endl;
}
BOOST_AUTO_TEST_SUITE_END()
