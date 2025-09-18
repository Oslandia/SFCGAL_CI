// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TORUS_H_
#define SFCGAL_TORUS_H_

#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/export.h"
#include "SFCGAL/primitive3d/Primitive.h"

#include <CGAL/Surface_mesh.h>
#include <optional>

namespace SFCGAL {

/**
 * @class Torus
 * @brief Represents a torus in 3D space
 *
 * This class provides methods to generate a polyhedral surface
 * representation of a torus. It uses SFCGAL's Kernel for
 * exact computations.
 */
class SFCGAL_API Torus : public Primitive {
public:
  /**
   * @brief Constructs a Torus object centered on the origin and lying in
   * the XY plane
   *
   * @param main_radius the main radius of the torus
   * @param tube_radius the tube radius of the torus
   * @param main_num_radial The number of main radial divisions
   * @param tube_num_radial The number of tube radial divisions
   */
  Torus(const Kernel::FT &main_radius = 10.0,
        const Kernel::FT &tube_radius = 2.0, unsigned int main_num_radial = 32,
        unsigned int tube_num_radial = 16);

  /**
   * @brief Copy constructor
   */
  Torus(const Torus &other) = default;

  [[nodiscard]] auto
  primitiveType() const -> std::string override;

  [[nodiscard]] auto
  primitiveTypeId() const -> PrimitiveType override;

  /**
   * @brief Assignment operator
   */
  auto
  operator=(Torus &other) -> Torus &;

  /**
   * @brief Destructor
   */
  ~Torus() = default;

  /**
   * @brief Sets the main radius of the torus
   * @param main_radius The new main radius
   */
  void
  setMainRadius(const Kernel::FT &main_radius);

  /**
   * @brief Sets the tube radius of the torus
   * @param tube_radius The new tube radius
   */
  void
  setTubeRadius(const Kernel::FT &tube_radius);

  /**
   * @brief Sets the main number of radial divisions
   * @param main_num_radial The new number of main radial
   */
  void
  setMainNumRadial(const unsigned int &main_num_radial);

  /**
   * @brief Sets the number of radial divisions for the tube
   * @param tube_num_radial The new number of tube radial divisions
   */
  void
  setTubeNumRadial(const unsigned int &tube_num_radial);

  /**
   * @brief Gets the main radius of the torus
   * @return The main radius
   */
  [[nodiscard]] auto
  mainRadius() const -> const Kernel::FT &;

  /**
   * @brief Gets the tube radius of the torus
   * @return The tube radius
   */
  [[nodiscard]] auto
  tubeRadius() const -> const Kernel::FT &;

  /**
   * @brief Gets the number of main radial divisions
   * @return The number of main radial divisions
   */
  [[nodiscard]] auto
  mainNumRadial() const -> unsigned int;

  /**
   * @brief Gets the number of radial divisions for the tube
   * @return The number of tube radial divisions
   */
  [[nodiscard]] auto
  tubeNumRadial() const -> unsigned int;

  /**
   * @brief Generates a surface mesh representation of the torus
   * @return A CGAL::Surface_mesh object representing the torus
   */
  auto
  generatePolyhedralSurface() const -> PolyhedralSurface override;

  /**
   * @brief Returns the perfect torus volume (without discretization)
   * @return The perfect torus volume (without discretization)
   */
  [[nodiscard]] auto
  volume(bool withDiscretization = false) const -> double override;

  /**
   * @brief Returns the perfect torus area (without discretization)
   * @return The perfect torus area (without discretization)
   */
  [[nodiscard]] auto
  area3D(bool withDiscretization = false) const -> double override;

protected:
  /**
   * @brief Verifies that all parameters are valid. For instance, it raises an
   * error if a radius is negative.
   * @throws SFCGAL::Exception if one of the parameters if not valid
   * provided variant type is not compatible with the parameter.
   */
  void
  validateParameters(std::unordered_map<std::string, PrimitiveParameter> const
                         &tempParameters) const override;
};

} // namespace SFCGAL

#endif // SFCGAL_TORUS_H_
