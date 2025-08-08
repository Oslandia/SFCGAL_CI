// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_CONE_H_
#define SFCGAL_CONE_H_

#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/export.h"
#include <CGAL/Surface_mesh.h>
#include <optional>

namespace SFCGAL {

/**
 * @class Cone
 * @brief Represents a cone in 3D space
 *
 * This class provides methods to generate a polyhedral surface
 * representation of a cone. It uses SFCGAL's Kernel for
 * exact computations.
 */
class SFCGAL_API Cone {
public:
  /**
   * @brief Constructs a Cone object with its base centered on the origin and in
   * the XY plane
   *
   * @param bottom_radius the bottom face radius of the cone
   * @param top_radius the top face radius of the cone
   * @param height the height of the cone
   * @param num_radial The number of radial divisions
   */
  Cone(Kernel::FT bottom_radius = 1.0, Kernel::FT top_radius = 0.0,
       Kernel::FT height = 1.0, int num_radial = 32);

  /**
   * @brief Copy constructor
   */
  Cone(const Cone &other) = default;

  /**
   * @brief Assignment operator
   */
  auto
  operator=(Cone other) -> Cone &;

  /**
   * @brief Destructor
   */
  ~Cone() = default;

  /**
   * @brief Sets the cone height
   * @param height The new cone height
   */
  void
  setHeight(const Kernel::FT &height);

  /**
   * @brief Sets the bottom face radius of the cone
   * @param bottom_radius The new bottom face radius
   */
  void
  setBottomRadius(const Kernel::FT &bottom_radius);

  /**
   * @brief Sets the top face radius of the cone
   * @param top_radius The new top face radius
   */
  void
  setTopRadius(const Kernel::FT &top_radius);

  /**
   * @brief Sets the number of radial divisions
   * @param num The new number of radial divisions
   */
  void
  setNumRadial(const int &num);

  /**
   * @brief Gets the cone height
   * @return The cone height
   */
  [[nodiscard]] auto
  height() const -> const Kernel::FT &;

  /**
   * @brief Gets the bottom face radius of the cone
   * @return The bottom face radius
   */
  [[nodiscard]] auto
  bottomRadius() const -> const Kernel::FT &;

  /**
   * @brief Gets the top face radius of the cone
   * @return The top face radius
   */
  [[nodiscard]] auto
  topRadius() const -> const Kernel::FT &;

  /**
   * @brief Gets the number of radial divisions
   * @return The number of radial divisions
   */
  [[nodiscard]] auto
  numRadial() const -> int;

  /**
   * @brief Generates a surface mesh representation of the cone
   * @return A CGAL::Surface_mesh object representing the cone
   */
  auto
  generatePolyhedralSurface() -> PolyhedralSurface;

  /**
   * @brief Returns the cone volume
   * @param withDiscretization If true, the volume is computed
   * using the real discretization with radial segments. If false, the volume is
   * computed for a perfect cone. Defaults to false.
   * @return The cone volume
   */
  [[nodiscard]] auto
  volume(bool withDiscretization = false) const -> double;

  /**
   * @brief Returns the cone area
   * @param withDiscretization If true, the area is computed
   * using the real discretization with radial segments. If false, the area is
   * computed for a perfect cone. Defaults to false.
   * @return The cone volume
   */
  [[nodiscard]] auto
  area(bool withDiscretization = false) const -> double;

private:
  Kernel::FT                       m_height;
  Kernel::FT                       m_bottom_radius;
  Kernel::FT                       m_top_radius;
  int                              m_num_radial;
  std::optional<PolyhedralSurface> m_polyhedral_surface;

  /**
   * @brief Invalidates the cached polyhedral surface
   */
  void
  invalidateCache();
};

} // namespace SFCGAL

#endif // SFCGAL_CONE_H_
