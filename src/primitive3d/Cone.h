// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_CONE_H_
#define SFCGAL_CONE_H_

#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/export.h"
#include "SFCGAL/primitive3d/Primitive.h"

#include <CGAL/Surface_mesh.h>

namespace SFCGAL {

/**
 * @class Cone
 * @brief Represents a cone in 3D space
 *
 * This class provides methods to generate a polyhedral surface
 * representation of a cone. It uses SFCGAL's Kernel for
 * exact computations.
 */
class SFCGAL_API Cone : public Primitive {
public:
  /**
   * @brief Constructs a Cone object with its base centered on the origin and in
   * the XY plane
   * A cone frustum is a cone with a top_radius greater than zero.
   *
   * @param bottom_radius the bottom face radius of the cone
   * @param top_radius the top face radius of the cone
   * @param height the height of the cone
   * @param num_radial The number of radial divisions
   */
  Cone(const Kernel::FT &bottom_radius = 1.0,
       const Kernel::FT &top_radius = 0.0, const Kernel::FT &height = 1.0,
       unsigned int num_radial = 32);

  /**
   * @brief Copy constructor
   * @param other copy from
   */
  Cone(const Cone &other) = default;

  [[nodiscard]] auto
  primitiveType() const -> std::string override;

  [[nodiscard]] auto
  primitiveTypeId() const -> PrimitiveType override;

  /**
   * @brief Assignment operator
   * @param other copy from
   * @return ref on this
   */
  auto
  operator=(Cone &other) -> Cone &;

  /**
   * @brief Destructor
   */
  ~Cone() override = default;

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
   * @param num_radial The new number of radial divisions
   */
  void
  setNumRadial(const unsigned int &num_radial);

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
  numRadial() const -> unsigned int;

  /**
   * @brief Generates a surface mesh representation of the cone
   * @return A CGAL::Surface_mesh object representing the cone
   */
  auto
  generatePolyhedralSurface() const -> PolyhedralSurface override;

  /**
   * @brief Returns the cone volume
   * @param withDiscretization If true, the volume is computed
   * using the real discretization with radial segments. If false, the volume is
   * computed for a perfect cone. Defaults to false.
   * @return The cone volume
   */
  [[nodiscard]] auto
  volume(bool withDiscretization = false) const -> double override;

  /**
   * @brief Returns the cone area
   * @param withDiscretization If true, the area is computed
   * using the real discretization with radial segments. If false, the area is
   * computed for a perfect cone. Defaults to false.
   * @return The cone volume
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
};

} // namespace SFCGAL

#endif // SFCGAL_CONE_H_
