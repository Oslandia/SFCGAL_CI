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
class SFCGAL_API Sphere : public Primitive {
public:
  /**
   * @brief Constructs a Sphere object
   * @param radius The radius of the sphere
   * @param center The center point of the sphere
   * @param num_vertical The number of vertical divisions
   * @param num_horizontal The number of horizontal divisions
   * @param direction The direction vector for sphere orientation
   */
  Sphere(const Kernel::FT      &radius = 1.0,
         const Kernel::Point_3 &center = Kernel::Point_3(0, 0, 0),
         unsigned int num_vertical = 16, unsigned int num_horizontal = 32,
         const Kernel::Vector_3 &direction = Kernel::Vector_3(0, 0, 1));

  /**
   * @brief Copy constructor
   * @param other copy from
   */
  Sphere(const Sphere &other) = default;

  /**
   * @brief Assignment operator
   * @param other copy from
   * @return ref on this
   */
  auto
  operator=(Sphere other) -> Sphere &;

  /**
   * @brief Destructor
   */
  ~Sphere() override = default;

  [[nodiscard]] auto
  primitiveType() const -> std::string override;

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
   * @brief Sets the number of vertical divisions
   * @param num The new number of vertical divisions
   */
  void
  setNumVertical(unsigned int num)
  {
    validateAndSetParameter("num_vertical", num);
  }

  /**
   * @brief Sets the number of horizontal divisions
   * @param num The new number of horizontal divisions
   */
  void
  setNumHorizontal(unsigned int num)
  {
    validateAndSetParameter("num_horizontal", num);
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
   * @brief Gets the number of vertical divisions
   * @return The number of vertical divisions
   */
  [[nodiscard]] auto
  numVertical() const -> unsigned int
  {
    return std::get<unsigned int>(m_parameters.at("num_vertical"));
  }

  /**
   * @brief Gets the number of horizontal divisions
   * @return The number of horizontal divisions
   */
  [[nodiscard]] auto
  numHorizontal() const -> unsigned int
  {
    return std::get<unsigned int>(m_parameters.at("num_horizontal"));
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
   * @copydoc SFCGAL::Primitive::volume
   * @note only perfect primitive version is available
   */
  auto
  volume(bool withDiscretization = false) const -> double override;

  /**
   * @copydoc SFCGAL::Primitive::area3D
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
