#include <SFCGAL/Cylinder.h>
#include <SFCGAL/Kernel.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Sphere.h>
#include <SFCGAL/algorithm/boolean3D.h>
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

BOOST_AUTO_TEST_CASE(testCigale)
{
  // Slender body
  Sphere body_sphere(20);
  auto body = std::make_unique<PolyhedralSurface>(body_sphere.generatePolyhedron());
  algorithm::scale(*body, 1.2, 0.5, 0.4);

  // Head
  Sphere head_sphere(17);
  auto head = std::make_unique<PolyhedralSurface>(head_sphere.generatePolyhedron());
  algorithm::scale(*head, 0.4, 0.6, 0.4);
  algorithm::translate(*head, 18.0, 0.0, 0.0);
  auto cigale = Boolean3D::union3D(*body, *head);

  // Torso
  Sphere torso_sphere(18);
  auto torso = std::make_unique<PolyhedralSurface>(torso_sphere.generatePolyhedron());
  algorithm::scale(*torso, 0.6, 0.6, 0.4);
  algorithm::translate(*torso, 11.0, 0.0, 0.0);
  cigale = Boolean3D::union3D(*cigale, *torso);

  // Eyes
  Sphere eye_sphere(3.5);
  auto eye_right = std::make_unique<PolyhedralSurface>(eye_sphere.generatePolyhedron());
  algorithm::translate(*eye_right, 18.75, 6.5, 1.5);
  cigale = Boolean3D::union3D(*cigale, *eye_right);

  auto eye_left = std::make_unique<PolyhedralSurface>(eye_sphere.generatePolyhedron());
  algorithm::translate(*eye_left, 18.75, -6.5, 1.5);
  cigale = Boolean3D::union3D(*cigale, *eye_left);

  // Wings
  // Complex part: Creating wings by subtracting an inner sphere from an outer sphere
  auto createWing = []() {
    Sphere outer_sphere(12, Point_3(0, 0, 0));
    auto outer = std::make_unique<PolyhedralSurface>(outer_sphere.generatePolyhedron());
    algorithm::scale(*outer, 3, 0.9, 0.7);

    Sphere inner_sphere(12, Point_3(0, 0, -3));
    auto inner = std::make_unique<PolyhedralSurface>(inner_sphere.generatePolyhedron());
    algorithm::scale(*inner, 3.5, 1, 0.9);

    return Boolean3D::difference3D(*outer, *inner);
  };

  auto wing_right = createWing();
  algorithm::translate(*wing_right, -20, 1.5, 4);
  algorithm::rotateX(*wing_right, -55 * M_PI / 180);
  algorithm::rotateY(*wing_right, -5 * M_PI / 180);
  algorithm::rotateZ(*wing_right, 2 * M_PI / 180);
  cigale = Boolean3D::union3D(*cigale, *wing_right);

  auto wing_left = createWing();
  algorithm::translate(*wing_left, -20, -1.5, 4);
  algorithm::rotateX(*wing_left, 55 * M_PI / 180);
  algorithm::rotateY(*wing_left, -5 * M_PI / 180);
  algorithm::rotateZ(*wing_left, -2 * M_PI / 180);
  cigale = Boolean3D::union3D(*cigale, *wing_left);

  // Complex part: Creating legs with two segments and proper bending
  auto createLeg = [](double x, double y, double z, double angle_horizontal, double angle_vertical, bool reverse_bend) {
    // Helper function to convert degrees to radians
    auto toRadians = [](double degrees) { return degrees * M_PI / 180.0; };

    // Helper function to create a cylinder between two points
    auto createCylinder = [](const Point_3& start, const Point_3& end, double radius) {
      Vector_3 axis = end - start;
      double height = std::sqrt(CGAL::to_double(axis.squared_length()));
      Cylinder cylinder(start, axis, radius, height);
      return std::make_unique<PolyhedralSurface>(cylinder.generatePolyhedron());
    };

    // First segment (thigh)
    double thigh_length = reverse_bend ? 8.0 : 10.0;
    double thigh_radius = 1.0;

    double cx = std::cos(toRadians(angle_horizontal)) * std::cos(toRadians(angle_vertical));
    double cy = std::sin(toRadians(angle_horizontal));
    double cz = -std::sin(toRadians(angle_vertical));

    Point_3 thigh_start(x, y, z);
    Point_3 thigh_end(x + thigh_length * cx, y + thigh_length * cy, z + thigh_length * cz);
    auto upper = createCylinder(thigh_start, thigh_end, thigh_radius);

    // Second segment (tibia) - bent relative to the first
    double tibia_length = reverse_bend ? 6.0 : 8.0;
    double tibia_radius = 0.8;
    double knee_angle = reverse_bend ? -30.0 : 30.0;  // Knee bending angle

    // Calculate the new direction for the tibia
    double new_angle_vertical = angle_vertical + knee_angle;
    cx = std::cos(toRadians(angle_horizontal)) * std::cos(toRadians(new_angle_vertical));
    cy = std::sin(toRadians(angle_horizontal));
    cz = -std::sin(toRadians(new_angle_vertical));

    Point_3 tibia_end(thigh_end.x() + tibia_length * cx, 
                      thigh_end.y() + tibia_length * cy, 
                      thigh_end.z() + tibia_length * cz);
    auto lower = createCylinder(thigh_end, tibia_end, tibia_radius);

    return Boolean3D::union3D(*upper, *lower);
  };

  // Legs
  // Front legs (reversed bending)
  cigale = Boolean3D::union3D(*cigale, *createLeg(15, 5, -4, 30, 45, true));
  cigale = Boolean3D::union3D(*cigale, *createLeg(15, -5, -4, -30, 45, true));

  // Middle legs (unchanged)
  cigale = Boolean3D::union3D(*cigale, *createLeg(5, 6, -4, 20, 30, false));
  cigale = Boolean3D::union3D(*cigale, *createLeg(5, -6, -4, -20, 30, false));

  // Rear legs (moved back)
  cigale = Boolean3D::union3D(*cigale, *createLeg(-15, 7, -2, 30, 15, false));
  cigale = Boolean3D::union3D(*cigale, *createLeg(-15, -7, -2, -30, 15, false));

  // Antennae
  Point_3 antenna_base_right(22, 2, 5);
  Vector_3 antenna_dir_right(std::cos(65 * M_PI / 180), std::sin(10 * M_PI / 180), std::sin(65 * M_PI / 180));
  Cylinder antenna_right(antenna_base_right, antenna_dir_right, 0.4, 15);
  auto antenna_r = std::make_unique<PolyhedralSurface>(antenna_right.generatePolyhedron());
  cigale = Boolean3D::union3D(*cigale, *antenna_r);

  Point_3 antenna_base_left(22, -2, 5);
  Vector_3 antenna_dir_left(std::cos(65 * M_PI / 180), -std::sin(10 * M_PI / 180), std::sin(65 * M_PI / 180));
  Cylinder antenna_left(antenna_base_left, antenna_dir_left, 0.4, 15);
  auto antenna_l = std::make_unique<PolyhedralSurface>(antenna_left.generatePolyhedron());
  cigale = Boolean3D::union3D(*cigale, *antenna_l);

  // Save as OBJ file
  io::OBJ::save(*cigale, "/tmp/cigale.obj");

  BOOST_CHECK(!cigale->isEmpty());
}
