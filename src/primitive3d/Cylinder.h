// Copyright (c) 2024-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_CYLINDER_H_
#define SFCGAL_CYLINDER_H_

#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/export.h"
#include "SFCGAL/primitive3d/Primitive.h"

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
class SFCGAL_API Cylinder : public PrimitiveImpl<Cylinder, Primitive> {
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
   * @param other copy from
   */
  Cylinder(const Cylinder &other) = default;

  /**
   * @brief returns the primitive type
   * @warning use CamelCase (Cylinder, not CYLINDER)
   * @return the primitive type as string
   */
  [[nodiscard]] auto
  primitiveType() const -> std::string override;

  /**
   * @brief returns a code corresponding to the type
   * @return a code corresponding to the type
   */
  [[nodiscard]] auto
  primitiveTypeId() const -> PrimitiveType override;

  /**
   * @brief Destructor
   */
  ~Cylinder() override = default;

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
    return std::get<Point_3>(m_parameters.at("base_center"));
  }

  /**
   * @brief Gets the axis of the cylinder
   * @return The axis vector
   */
  [[nodiscard]] auto
  axis() const -> const Vector_3 &
  {
    return std::get<Vector_3>(m_parameters.at("axis"));
  }

  /**
   * @brief Gets the radius of the cylinder
   * @return The radius
   */
  [[nodiscard]] auto
  radius() const -> const Kernel::FT &
  {
    return std::get<Kernel::FT>(m_parameters.at("radius"));
  }

  /**
   * @brief Gets the height of the cylinder
   * @return The height
   */
  [[nodiscard]] auto
  height() const -> const Kernel::FT &
  {
    return std::get<Kernel::FT>(m_parameters.at("height"));
  }

  /**
   * @brief Gets the number of radial divisions
   * @return The number of radial divisions
   */
  [[nodiscard]] auto
  numRadial() const -> unsigned int
  {
    return std::get<unsigned int>(m_parameters.at("num_radial"));
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
  generatePolyhedralSurface() const -> PolyhedralSurface override;

  /**
   * @brief Returns the primitive volume
   * @param withDiscretization Computes volume with discretization (true) or as
   * perfect primitive (false). Defaults to false.
   * @return The cylinder volume
   * @note only perfect primitive version is available
   */
  [[nodiscard]] auto
  volume(bool withDiscretization = false) const -> double override;

  /**
   * @brief Returns the primitive area
   * @param withDiscretization Computes area with discretization (true) or as
   * perfect primitive (false). Defaults to false.
   * @return The cylinder area
   * @note only perfect primitive version is available
   */
  [[nodiscard]] auto
  area3D(bool withDiscretization = false) const -> double override;

protected:
  /**
   * @brief Verifies that all parameters are valid. For instance, it raises an
   * error if a radius is negative.
   * @param tempParameters a temp map of parameter with new values
   * @throws SFCGAL::Exception if one of the parameters if not valid
   * provided variant type is not compatible with the parameter.
   */
  void
  validateParameters(std::unordered_map<std::string, PrimitiveParameter> const
                         &tempParameters) const override;

  /**
   * @brief Invalidates the cached polyhedron and surface mesh
   */
  void
  invalidateCache() override;

private:
  mutable std::optional<Polyhedron_3>   m_polyhedron;
  mutable std::optional<Surface_mesh_3> m_surface_mesh;

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
