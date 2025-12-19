// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_BOX_H_
#define SFCGAL_BOX_H_

#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/export.h"
#include "SFCGAL/primitive3d/Primitive.h"

#include <CGAL/Surface_mesh.h>

namespace SFCGAL {

/**
 * @class Box
 * @brief Represents a box in 3D space
 *
 * This class provides methods to generate a polyhedral surface
 * representation of a box (rectangular prism). It uses SFCGAL's Kernel for
 * exact computations.
 */
class SFCGAL_API Box : public Primitive {
public:
  /**
   * @brief Constructs a Box object in the first octant (+, +, +)
   * with one corner at (0, 0, 0)
   *
   * @param x_extent The box length in x direction
   * @param y_extent The box length in y direction
   * @param z_extent The box length in z direction
   */
  Box(const Kernel::FT &x_extent = 1.0, const Kernel::FT &y_extent = 1.0,
      const Kernel::FT &z_extent = 1.0);

  /**
   * @brief Copy constructor
   * @param other copy from
   */
  Box(const Box &other) = default;

  [[nodiscard]] auto
  primitiveType() const -> std::string override;

  [[nodiscard]] auto
  primitiveTypeId() const -> PrimitiveType override;

  /**
   * @brief Destructor
   */
  ~Box() override = default;

  /**
   * @brief Sets the x_extent of the box
   * @param x_extent The new x_extent
   */
  void
  setXExtent(const Kernel::FT &x_extent);

  /**
   * @brief Sets the y_extent of the box
   * @param y_extent The new y_extent
   */
  void
  setYExtent(const Kernel::FT &y_extent);

  /**
   * @brief Sets the z_extent of the box
   * @param z_extent The new z_extent
   */
  void
  setZExtent(const Kernel::FT &z_extent);

  /**
   * @brief Gets the x_extent of the box
   * @return The x_extent
   */
  [[nodiscard]] auto
  xExtent() const -> const Kernel::FT &;

  /**
   * @brief Gets the y_extent of the box
   * @return The y_extent
   */
  [[nodiscard]] auto
  yExtent() const -> const Kernel::FT &;

  /**
   * @brief Gets the z_extent of the box
   * @return The z_extent
   */
  [[nodiscard]] auto
  zExtent() const -> const Kernel::FT &;

  /**
   * @brief Generates a surface mesh representation of the box
   * @return A CGAL::Surface_mesh object representing the box
   */
  auto
  generatePolyhedralSurface() const -> PolyhedralSurface override;

  /**
   * @copydoc SFCGAL::Primitive::volume
   * @note only perfect primitive version is available
   */
  [[nodiscard]] auto
  volume(bool withDiscretization = false) const -> double override;

  /**
   * @copydoc SFCGAL::Primitive::area3D
   * @note only perfect primitive version is available
   */
  [[nodiscard]] auto
  area3D(bool withDiscretization = false) const -> double override;

protected:
  /**
   * @brief Verifies that all parameters are valid. For instance, it raises an
   * error if an extent is negative.
   * @param tempParameters a temp map of parameter with new values
   * @throws SFCGAL::Exception if one of the parameters if not valid
   * provided variant type is not compatible with the parameter.
   */
  void
  validateParameters(std::unordered_map<std::string, PrimitiveParameter> const
                         &tempParameters) const override;
};

} // namespace SFCGAL

#endif // SFCGAL_BOX_H_
