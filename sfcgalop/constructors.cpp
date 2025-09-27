#include "constructors.hpp"

#include <SFCGAL/Kernel.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/primitive3d/Box.h>
#include <SFCGAL/primitive3d/Cone.h>
#include <SFCGAL/primitive3d/Cylinder.h>
#include <SFCGAL/primitive3d/Sphere.h>
#include <SFCGAL/primitive3d/Torus.h>

namespace Constructors {

auto
make_sphere(double x, double y, double z, double radius,
            unsigned int num_vertical, unsigned int num_horizontal)
    -> std::unique_ptr<SFCGAL::Geometry>
{
  SFCGAL::Kernel::Point_3  center(x, y, z);
  SFCGAL::Kernel::Vector_3 direction(0, 0, 1); // Default up direction

  SFCGAL::Sphere sphere(radius, center, num_vertical, num_horizontal,
                        direction);

  auto polyhedral_surface = sphere.generatePolyhedralSurface();
  return std::make_unique<SFCGAL::PolyhedralSurface>(
      std::move(polyhedral_surface));
}

auto
make_box(double x_extent, double y_extent, double z_extent)
    -> std::unique_ptr<SFCGAL::Geometry>
{
  SFCGAL::Box box(x_extent, y_extent, z_extent);

  auto polyhedral_surface = box.generatePolyhedralSurface();
  return std::make_unique<SFCGAL::PolyhedralSurface>(
      std::move(polyhedral_surface));
}

auto
make_cube(double size) -> std::unique_ptr<SFCGAL::Geometry>
{
  return make_box(size, size, size);
}

auto
make_cylinder(double base_x, double base_y, double base_z, double axis_x,
              double axis_y, double axis_z, double radius, double height,
              unsigned int num_radial) -> std::unique_ptr<SFCGAL::Geometry>
{
  SFCGAL::Kernel::Point_3  base_center(base_x, base_y, base_z);
  SFCGAL::Kernel::Vector_3 axis(axis_x, axis_y, axis_z);

  SFCGAL::Cylinder cylinder(base_center, axis, radius, height, num_radial);

  auto polyhedral_surface = cylinder.generatePolyhedralSurface();
  return std::make_unique<SFCGAL::PolyhedralSurface>(
      std::move(polyhedral_surface));
}

auto
make_cone(double base_x, double base_y, double base_z, double axis_x,
          double axis_y, double axis_z, double bottom_radius, double top_radius,
          double height, unsigned int num_radial)
    -> std::unique_ptr<SFCGAL::Geometry>
{
  // Note: SFCGAL::Cone constructor doesn't take center/axis, only dimensions
  // Parameters base_x, base_y, base_z, axis_x, axis_y, axis_z are for future
  // use
  (void)base_x;
  (void)base_y;
  (void)base_z;
  (void)axis_x;
  (void)axis_y;
  (void)axis_z;

  SFCGAL::Cone cone(bottom_radius, top_radius, height, num_radial);

  auto polyhedral_surface = cone.generatePolyhedralSurface();
  return std::make_unique<SFCGAL::PolyhedralSurface>(
      std::move(polyhedral_surface));
}

auto
make_torus(double center_x, double center_y, double center_z, double axis_x,
           double axis_y, double axis_z, double major_radius,
           double minor_radius, unsigned int num_major, unsigned int num_minor)
    -> std::unique_ptr<SFCGAL::Geometry>
{
  // Note: SFCGAL::Torus constructor doesn't take center/axis, only dimensions
  // Parameters center_x, center_y, center_z, axis_x, axis_y, axis_z are for
  // future use
  (void)center_x;
  (void)center_y;
  (void)center_z;
  (void)axis_x;
  (void)axis_y;
  (void)axis_z;

  SFCGAL::Torus torus(major_radius, minor_radius, num_major, num_minor);

  auto polyhedral_surface = torus.generatePolyhedralSurface();
  return std::make_unique<SFCGAL::PolyhedralSurface>(
      std::move(polyhedral_surface));
}

auto
make_solid(std::unique_ptr<SFCGAL::Geometry> polyhedralsurface)
    -> std::unique_ptr<SFCGAL::Geometry>
{
  // Check if input is a PolyhedralSurface
  auto *polyhedral_ptr =
      dynamic_cast<SFCGAL::PolyhedralSurface *>(polyhedralsurface.get());
  if (polyhedral_ptr == nullptr) {
    throw std::invalid_argument(
        "make_solid: input geometry must be a PolyhedralSurface");
  }

  // Create a copy of the PolyhedralSurface for the Solid
  auto polyhedral_copy =
      std::make_unique<SFCGAL::PolyhedralSurface>(*polyhedral_ptr);

  // Create and return the Solid with the PolyhedralSurface as exterior shell
  return std::make_unique<SFCGAL::Solid>(polyhedral_copy.release());
}

} // namespace Constructors