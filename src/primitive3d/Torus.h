// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TORUS_H_
#define SFCGAL_TORUS_H_

#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/export.h"
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
class SFCGAL_API Torus {
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
  Torus(Kernel::FT main_radius = 10.0, Kernel::FT tube_radius = 2.0,
        int main_num_radial = 32, int tube_num_radial = 16);

  /**
   * @brief Copy constructor
   */
  Torus(const Torus &other) = default;

  /**
   * @brief Assignment operator
   */
  auto
  operator=(Torus other) -> Torus &;

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
  setMainNumRadial(const int &main_num_radial);

  /**
   * @brief Sets the number of radial divisions for the tube
   * @param tube_num_radial The new number of tube radial divisions
   */
  void
  setTubeNumRadial(const int &tube_num_radial);

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
  mainNumRadial() const -> int;

  /**
   * @brief Gets the number of radial divisions for the tube
   * @return The number of tube radial divisions
   */
  [[nodiscard]] auto
  tubeNumRadial() const -> int;

  /**
   * @brief Generates a surface mesh representation of the torus
   * @return A CGAL::Surface_mesh object representing the torus
   */
  auto
  generatePolyhedralSurface() -> PolyhedralSurface;

  /**
   * @brief Returns the perfect torus volume (without discretization)
   * @return The perfect torus volume (without discretization)
   */
  [[nodiscard]] auto
  volume() const -> double;

  /**
   * @brief Returns the perfect torus area (without discretization)
   * @return The perfect torus area (without discretization)
   */
  [[nodiscard]] auto
  area() const -> double;

private:
  Kernel::FT                       m_main_radius;
  Kernel::FT                       m_tube_radius;
  int                              m_main_num_radial;
  int                              m_tube_num_radial;
  std::optional<PolyhedralSurface> m_polyhedral_surface;

  /**
   * @brief Invalidates the cached polyhedral surface
   */
  void
  invalidateCache();
};

} // namespace SFCGAL

#endif // SFCGAL_TORUS_H_
