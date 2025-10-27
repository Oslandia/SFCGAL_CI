// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "constructors.hpp"

#include <SFCGAL/Kernel.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/algorithm/rotate.h>
#include <SFCGAL/algorithm/translate.h>
#include <SFCGAL/primitive3d/Box.h>
#include <SFCGAL/primitive3d/Cone.h>
#include <SFCGAL/primitive3d/Cylinder.h>
#include <SFCGAL/primitive3d/Sphere.h>
#include <SFCGAL/primitive3d/Torus.h>

#include <algorithm>
#include <cmath>

namespace Constructors {

auto
make_sphere(double x, double y, double z, double radius,
            unsigned int num_subdivisions) -> std::unique_ptr<SFCGAL::Geometry>
{
  SFCGAL::Kernel::Point_3  center(x, y, z);
  SFCGAL::Kernel::Vector_3 direction(0, 0, 1); // Default up direction

  SFCGAL::Sphere sphere(radius, center, num_subdivisions, direction);

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
  // SFCGAL::Cone constructor creates a cone with base centered at origin,
  // oriented along +Z axis. We need to apply transformations for custom
  // position/orientation.

  SFCGAL::Cone cone(bottom_radius, top_radius, height, num_radial);
  auto         polyhedral_surface = cone.generatePolyhedralSurface();

  // Normalize the axis vector
  double axis_length =
      std::sqrt((axis_x * axis_x) + (axis_y * axis_y) + (axis_z * axis_z));

  // If axis is non-zero and not aligned with +Z, apply rotation
  constexpr double epsilon = 1e-10;
  if (axis_length > epsilon) {
    double nx = axis_x / axis_length;
    double ny = axis_y / axis_length;
    double nz = axis_z / axis_length;

    // Check if axis is not already aligned with +Z (0, 0, 1)
    if (std::abs(nx) > epsilon || std::abs(ny) > epsilon ||
        std::abs(nz - 1.0) > epsilon) {
      // Compute rotation axis: cross product of +Z with desired axis
      // +Z = (0, 0, 1), so cross product is (ny, -nx, 0) (normalized later)
      SFCGAL::Kernel::Vector_3 rotation_axis(ny, -nx, 0);

      // Compute rotation angle using dot product: +Z · normalized_axis = nz
      double angle = std::acos(std::clamp(nz, -1.0, 1.0));

      if (std::abs(angle) > epsilon) {
        SFCGAL::algorithm::rotate(polyhedral_surface, angle, rotation_axis);
      }
    }
  }

  // Apply translation to position the base at (base_x, base_y, base_z)
  if (std::abs(base_x) > epsilon || std::abs(base_y) > epsilon ||
      std::abs(base_z) > epsilon) {
    SFCGAL::algorithm::translate(polyhedral_surface, base_x, base_y, base_z);
  }

  return std::make_unique<SFCGAL::PolyhedralSurface>(
      std::move(polyhedral_surface));
}

auto
make_torus(double center_x, double center_y, double center_z, double axis_x,
           double axis_y, double axis_z, double major_radius,
           double minor_radius, unsigned int num_major, unsigned int num_minor)
    -> std::unique_ptr<SFCGAL::Geometry>
{
  // SFCGAL::Torus constructor creates a torus centered at origin,
  // oriented in the XY plane (normal along +Z axis).
  // We need to apply transformations for custom position/orientation.

  SFCGAL::Torus torus(major_radius, minor_radius, num_major, num_minor);
  auto          polyhedral_surface = torus.generatePolyhedralSurface();

  // Normalize the axis vector
  double axis_length =
      std::sqrt((axis_x * axis_x) + (axis_y * axis_y) + (axis_z * axis_z));

  // If axis is non-zero and not aligned with +Z, apply rotation
  constexpr double epsilon = 1e-10;
  if (axis_length > epsilon) {
    double nx = axis_x / axis_length;
    double ny = axis_y / axis_length;
    double nz = axis_z / axis_length;

    // Check if axis is not already aligned with +Z (0, 0, 1)
    if (std::abs(nx) > epsilon || std::abs(ny) > epsilon ||
        std::abs(nz - 1.0) > epsilon) {
      // Compute rotation axis: cross product of +Z with desired axis
      // +Z = (0, 0, 1), so cross product is (ny, -nx, 0)
      SFCGAL::Kernel::Vector_3 rotation_axis(ny, -nx, 0);

      // Compute rotation angle using dot product: +Z · normalized_axis = nz
      double angle = std::acos(std::clamp(nz, -1.0, 1.0));

      if (std::abs(angle) > epsilon) {
        SFCGAL::algorithm::rotate(polyhedral_surface, angle, rotation_axis);
      }
    }
  }

  // Apply translation to position the center at (center_x, center_y, center_z)
  if (std::abs(center_x) > epsilon || std::abs(center_y) > epsilon ||
      std::abs(center_z) > epsilon) {
    SFCGAL::algorithm::translate(polyhedral_surface, center_x, center_y,
                                 center_z);
  }

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
