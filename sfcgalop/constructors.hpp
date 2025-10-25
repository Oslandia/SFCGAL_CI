// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGALOP_CONSTRUCTORS_HPP
#define SFCGALOP_CONSTRUCTORS_HPP

#include <SFCGAL/Geometry.h>
#include <memory>
#include <string>

namespace Constructors {

/**
 * @brief Create a sphere primitive using icosahedron subdivision
 * @param x X coordinate of center (default: 0.0)
 * @param y Y coordinate of center (default: 0.0)
 * @param z Z coordinate of center (default: 0.0)
 * @param radius Sphere radius (default: 1.0)
 * @param num_subdivisions Number of icosahedron subdivisions (default: 2)
 * @return Unique pointer to PolyhedralSurface representing the sphere
 */
auto
make_sphere(double x = 0.0, double y = 0.0, double z = 0.0, double radius = 1.0,
            unsigned int num_subdivisions = 2)
    -> std::unique_ptr<SFCGAL::Geometry>;

/**
 * @brief Create a box/cube primitive
 * @param x_extent Length in X direction (default: 1.0)
 * @param y_extent Length in Y direction (default: 1.0)
 * @param z_extent Length in Z direction (default: 1.0)
 * @return Unique pointer to PolyhedralSurface representing the box
 */
auto
make_box(double x_extent = 1.0, double y_extent = 1.0, double z_extent = 1.0)
    -> std::unique_ptr<SFCGAL::Geometry>;

/**
 * @brief Create a cube primitive (equal extents in all directions)
 * @param size Edge length of the cube (default: 1.0)
 * @return Unique pointer to PolyhedralSurface representing the cube
 */
auto
make_cube(double size = 1.0) -> std::unique_ptr<SFCGAL::Geometry>;

/**
 * @brief Create a cylinder primitive
 * @param base_x X coordinate of base center (default: 0.0)
 * @param base_y Y coordinate of base center (default: 0.0)
 * @param base_z Z coordinate of base center (default: 0.0)
 * @param axis_x X component of cylinder axis (default: 0.0)
 * @param axis_y Y component of cylinder axis (default: 0.0)
 * @param axis_z Z component of cylinder axis (default: 1.0)
 * @param radius Cylinder radius (default: 1.0)
 * @param height Cylinder height (default: 1.0)
 * @param num_radial Number of radial divisions (default: 32)
 * @return Unique pointer to PolyhedralSurface representing the cylinder
 */
auto
make_cylinder(double base_x = 0.0, double base_y = 0.0, double base_z = 0.0,
              double axis_x = 0.0, double axis_y = 0.0, double axis_z = 1.0,
              double radius = 1.0, double height = 1.0,
              unsigned int num_radial = 32)
    -> std::unique_ptr<SFCGAL::Geometry>;

/**
 * @brief Create a cone primitive (supports truncated cones/frustum)
 * @param base_x X coordinate of base center (default: 0.0)
 * @param base_y Y coordinate of base center (default: 0.0)
 * @param base_z Z coordinate of base center (default: 0.0)
 * @param axis_x X component of cone axis (default: 0.0)
 * @param axis_y Y component of cone axis (default: 0.0)
 * @param axis_z Z component of cone axis (default: 1.0)
 * @param bottom_radius Cone bottom radius (default: 1.0)
 * @param top_radius Cone top radius - 0.0 for regular cone (default: 0.0)
 * @param height Cone height (default: 1.0)
 * @param num_radial Number of radial divisions (default: 32)
 * @return Unique pointer to PolyhedralSurface representing the cone
 */
auto
make_cone(double base_x = 0.0, double base_y = 0.0, double base_z = 0.0,
          double axis_x = 0.0, double axis_y = 0.0, double axis_z = 1.0,
          double bottom_radius = 1.0, double top_radius = 0.0,
          double height = 1.0, unsigned int num_radial = 32)
    -> std::unique_ptr<SFCGAL::Geometry>;

/**
 * @brief Create a torus primitive
 * @param center_x X coordinate of center (default: 0.0)
 * @param center_y Y coordinate of center (default: 0.0)
 * @param center_z Z coordinate of center (default: 0.0)
 * @param axis_x X component of torus axis (default: 0.0)
 * @param axis_y Y component of torus axis (default: 0.0)
 * @param axis_z Z component of torus axis (default: 1.0)
 * @param major_radius Major radius of torus (default: 2.0)
 * @param minor_radius Minor radius of torus (default: 0.5)
 * @param num_major Number of major divisions (default: 32)
 * @param num_minor Number of minor divisions (default: 16)
 * @return Unique pointer to PolyhedralSurface representing the torus
 */
auto
make_torus(double center_x = 0.0, double center_y = 0.0, double center_z = 0.0,
           double axis_x = 0.0, double axis_y = 0.0, double axis_z = 1.0,
           double major_radius = 2.0, double minor_radius = 0.5,
           unsigned int num_major = 32, unsigned int num_minor = 16)
    -> std::unique_ptr<SFCGAL::Geometry>;

/**
 * @brief Convert a PolyhedralSurface to a Solid
 * @param polyhedralsurface Input PolyhedralSurface geometry
 * @return Unique pointer to Solid representing the input as a solid volume
 */
auto
make_solid(std::unique_ptr<SFCGAL::Geometry> polyhedralsurface)
    -> std::unique_ptr<SFCGAL::Geometry>;

} // namespace Constructors

#endif // SFCGALOP_CONSTRUCTORS_HPP
