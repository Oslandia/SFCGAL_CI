// Copyright (c) 2024-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_SPHERE_H_
#define SFCGAL_SPHERE_H_

#include <optional>
#include <vector>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/squared_distance_3.h>

#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/export.h"
#include "SFCGAL/numeric.h"
#include "SFCGAL/primitive3d/Primitive.h"

namespace SFCGAL {

class PolyhedralSurface;

/**
 * @brief Represents a sphere in 3D space
 *
 * This class provides methods to generate a polyhedron and a point cloud
 * representation of a sphere. It uses SFCGAL's Kernel for exact computations.
 */
class SFCGAL_API Sphere : public PrimitiveImpl<Sphere, Primitive> {
public:
  /**
   * @brief Constructs a Sphere object
   * @param radius The radius of the sphere
   * @param center The center point of the sphere
   * @param num_subdivisions The number of icosahedron subdivisions (0=12
   * vertices, 1=42, 2=162, etc.)
   * @param direction The direction vector for sphere orientation
   */
  Sphere(const Kernel::FT       &radius           = 1.0,
         const Kernel::Point_3  &center           = Kernel::Point_3(0, 0, 0),
         unsigned int            num_subdivisions = 2,
         const Kernel::Vector_3 &direction        = Kernel::Vector_3(0, 0, 1));

  /**
   * @brief Copy constructor
   * @param other copy from
   */
  Sphere(const Sphere &other) = default;

  /**
   * @brief Destructor
   */
  ~Sphere() override = default;

  /**
   * @brief returns the primitive type
   * @warning use CamelCase (Sphere, not SPHERE)
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
   * @brief Sets the radius of the sphere
   * @param radius The new radius
   */
  void
  setRadius(const Kernel::FT &radius)
  {
    validateAndSetParameter("radius", radius);
  }

  /**
   * @brief Sets the center of the sphere
   * @param center The new center point
   */
  void
  setCenter(const Kernel::Point_3 &center)
  {
    validateAndSetParameter("center", center);
  }

  /**
   * @brief Sets the number of subdivisions
   * @param num The new number of subdivisions
   */
  void
  setNumSubdivisions(unsigned int num)
  {
    validateAndSetParameter("num_subdivisions", num);
  }

  /**
   * @brief Sets the direction of the sphere
   * @param direction The new direction vector
   */
  void
  setDirection(const Kernel::Vector_3 &direction)
  {
    validateAndSetParameter("direction", direction);
  }

  /**
   * @brief Gets the radius of the sphere
   * @return The radius
   */
  [[nodiscard]] auto
  radius() const -> const Kernel::FT &
  {
    return std::get<Kernel::FT>(m_parameters.at("radius"));
  }

  /**
   * @brief Gets the center of the sphere
   * @return The center point
   */
  [[nodiscard]] auto
  center() const -> const Kernel::Point_3 &
  {
    return std::get<Kernel::Point_3>(m_parameters.at("center"));
  }

  /**
   * @brief Gets the number of subdivisions
   * @return The number of subdivisions
   */
  [[nodiscard]] auto
  numSubdivisions() const -> unsigned int
  {
    return std::get<unsigned int>(m_parameters.at("num_subdivisions"));
  }

  /**
   * @brief Gets the direction of the sphere
   * @return The direction vector
   */
  [[nodiscard]] auto
  direction() const -> const Kernel::Vector_3 &
  {
    return std::get<Kernel::Vector_3>(m_parameters.at("direction"));
  }

  /**
   * @brief Generates a polyhedron representation of the sphere
   * @return A CGAL::Polyhedron_3 object representing the sphere
   */
  auto
  generatePolyhedron() const -> CGAL::Polyhedron_3<Kernel>;

  /**
   * @brief Generates a point cloud representation of the sphere
   * @return A vector of Point_3 objects representing points on the sphere's
   * surface
   */
  auto
  generatePoints() const -> std::vector<Kernel::Point_3>;

  /**
   * @brief Generates a surface mesh representation of the sphere
   * @return A SFCGAL::PolyhedralSurface object representing the sphere
   */
  auto
  generatePolyhedralSurface() const -> PolyhedralSurface override;

  /**
   * @brief Returns the sphere volume
   * @param withDiscretization Computes volume with discretization (true) or as
   * perfect primitive (false). Defaults to false.
   * @return The sphere volume
   * @note only perfect primitive version is available
   */
  auto
  volume(bool withDiscretization = false) const -> double override;

  /**
   * @brief Returns the sphere area
   * @param withDiscretization Computes area with discretization (true) or as
   * perfect primitive (false). Defaults to false.
   * @return The sphere area
   * @note only perfect primitive version is available
   */
  auto
  area3D(bool withDiscretization = false) const -> double override;

protected:
  /**
   * @brief Invalidates the cached geometries
   */
  void
  invalidateCache() override;

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

private:
  mutable std::optional<CGAL::Polyhedron_3<Kernel>>   m_polyhedron;
  mutable std::optional<std::vector<Kernel::Point_3>> m_points;
  mutable std::optional<PolyhedralSurface>            m_polyhedral_surface;

  auto
  generateSpherePolyhedron() const -> CGAL::Polyhedron_3<Kernel>;

  auto
  generateSpherePoints() const -> std::vector<Kernel::Point_3>;

  // Function to get an orthogonal vector in the XY plane
  static auto
  get_orthogonal_vector(const Kernel::Vector_3 &vec) -> Kernel::Vector_3
  {
    if (vec.x() != 0 || vec.y() != 0) {
      return {-vec.y(), vec.x(), 0};
    }
    return {0, -vec.z(), vec.y()};
  }
};

} // namespace SFCGAL

#endif // SFCGAL_SPHERE_H_
