// Copyright (c) 2024-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_CYLINDER_H_
#define SFCGAL_CYLINDER_H_

#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/export.h"
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <optional>

namespace SFCGAL {

/**
 * @class Cylinder
 * @brief Represents a cylinder in 3D space
 *
 * This class provides methods to generate a polyhedron and a surface mesh
 * representation of a cylinder. It uses SFCGAL's Kernel for exact computations.
 */
class SFCGAL_API Cylinder {
public:
  /**
   * @brief Constructs a Cylinder object
   * @param base_center The center point of the base of the cylinder
   * @param axis The axis of the cylinder
   * @param radius The radius of the cylinder
   * @param height The height of the cylinder
   * @param num_radial The number of radial divisions
   */
  Cylinder(const Point_3    &base_center = Point_3(0, 0, 0),
           const Vector_3   &axis        = Vector_3(0, 0, 1),
           const Kernel::FT &radius = 1.0, const Kernel::FT &height = 1.0,
           unsigned int num_radial = 32);

  /**
   * @brief Copy constructor
   */
  Cylinder(const Cylinder &other) = default;

  /**
   * @brief Assignment operator
   */
  auto
  operator=(Cylinder other) -> Cylinder &;

  /**
   * @brief Destructor
   */
  ~Cylinder() = default;

  /**
   * @brief Sets the base center of the cylinder
   * @param base_center The new base center point
   */
  void
  setBaseCenter(const Point_3 &base_center);

  /**
   * @brief Sets the axis of the cylinder
   * @param axis The new axis vector
   */
  void
  setAxis(const Vector_3 &axis);

  /**
   * @brief Sets the radius of the cylinder
   * @param radius The new radius
   */
  void
  setRadius(const Kernel::FT &radius);

  /**
   * @brief Sets the height of the cylinder
   * @param height The new height
   */
  void
  setHeight(const Kernel::FT &height);

  /**
   * @brief Sets the number of radial divisions
   * @param num The new number of radial divisions
   */
  void
  setNumRadial(unsigned int num);

  /**
   * @brief Gets the base center of the cylinder
   * @return The base center point
   */
  [[nodiscard]] auto
  baseCenter() const -> const Point_3 &
  {
    return m_base_center;
  }

  /**
   * @brief Gets the axis of the cylinder
   * @return The axis vector
   */
  [[nodiscard]] auto
  axis() const -> const Vector_3 &
  {
    return m_axis;
  }

  /**
   * @brief Gets the radius of the cylinder
   * @return The radius
   */
  [[nodiscard]] auto
  radius() const -> const Kernel::FT &
  {
    return m_radius;
  }

  /**
   * @brief Gets the height of the cylinder
   * @return The height
   */
  [[nodiscard]] auto
  height() const -> const Kernel::FT &
  {
    return m_height;
  }

  /**
   * @brief Gets the number of radial divisions
   * @return The number of radial divisions
   */
  [[nodiscard]] auto
  numRadial() const -> unsigned int
  {
    return m_num_radial;
  }

  /**
   * @brief Generates a polyhedron representation of the cylinder
   * @return A CGAL::Polyhedron_3 object representing the cylinder
   */
  auto
  generatePolyhedron() const -> Polyhedron_3;

  /**
   * @brief Generates a surface mesh representation of the cylinder
   * @return A CGAL::Surface_mesh object representing the cylinder
   */
  auto
  generateSurfaceMesh() const -> Surface_mesh_3;

  /**
   * @brief Generates a surface mesh representation of the cylinder
   * @return A SFCGAL::PolyhedralSurface object representing the cylinder
   */
  auto
  generatePolyhedralSurface() const -> PolyhedralSurface;

  [[nodiscard]] auto
  volume() const -> double
  {
    return CGAL::to_double(m_radius * m_radius * m_height * CGAL_PI);
  }

  [[nodiscard]] auto
  area3D() const -> double
  {
    return CGAL::to_double(2 * m_radius * m_radius * CGAL_PI +
                           2 * m_radius * m_height * CGAL_PI);
  }

private:
  Point_3                                  m_base_center;
  Vector_3                                 m_axis;
  Kernel::FT                               m_radius;
  Kernel::FT                               m_height;
  unsigned int                             m_num_radial;
  mutable std::optional<Polyhedron_3>      m_polyhedron;
  mutable std::optional<Surface_mesh_3>    m_surface_mesh;
  mutable std::optional<PolyhedralSurface> m_polyhedral_surface;

  /**
   * @brief Invalidates the cached polyhedron and surface mesh
   */
  void
  invalidateCache();

  /**
   * @brief Normalizes a vector
   * @param vector The vector to normalize
   * @return The normalized vector
   */
  static auto
  normalize(const Vector_3 &vector) -> Vector_3;
};

} // namespace SFCGAL

#endif // SFCGAL_CYLINDER_H_
