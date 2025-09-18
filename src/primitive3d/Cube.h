// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_CUBE_H_
#define SFCGAL_CUBE_H_

#include "Box.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/export.h"

namespace SFCGAL {

/**
 * @class Cube
 * @brief Represents a cube in 3D space, which is a box where all 3 sides are of
 * equal size
 *
 * This class provides methods to generate a polyhedral surface
 * representation of a cube. It uses SFCGAL's Kernel for
 * exact computations.
 */
class SFCGAL_API Cube : public Primitive {
public:
  /**
   * @brief Constructs a Cube object in the first octant (+, +, +)
   * with one corner at (0, 0, 0)
   *
   * @param size The cube size
   */
  Cube(const Kernel::FT &size = 1.0);

  /**
   * @brief Copy constructor
   */
  Cube(const Cube &other) = default;

  [[nodiscard]] auto
  primitiveType() const -> std::string override;

  [[nodiscard]] auto
  primitiveTypeId() const -> PrimitiveType override;

  /**
   * @brief Assignment operator
   */
  auto
  operator=(Cube &other) -> Cube &;

  /**
   * @brief Destructor
   */
  ~Cube() override = default;

  /**
   * @brief Sets the size of the cube
   * @param size The new size
   */
  void
  setSize(const Kernel::FT &size);

  /**
   * @brief Gets the size of the cube
   * @return The size
   */
  [[nodiscard]] auto
  size() const -> const Kernel::FT &;

  /**
   * @brief Generates a surface mesh representation of the cube
   * @return A CGAL::Surface_mesh object representing the cube
   */
  auto
  generatePolyhedralSurface() const -> PolyhedralSurface override;

  /**
   * @brief Returns the cube volume
   * @return The cube volume
   */
  [[nodiscard]] auto
  volume(bool withDiscretization = false) const -> double override;

  /**
   * @brief Returns the cube area
   * @return The cube area
   */
  [[nodiscard]] auto
  area3D(bool withDiscretization = false) const -> double override;

  /**
   * Returns string representation of this object.
   *
   * \return string representation of this object
   */
  [[nodiscard]] virtual auto
  toString() const -> std::string override;

protected:
  /**
   * @brief Verifies that all parameters are valid. For instance, it raises an
   * error if an extent is negative.
   * @throws SFCGAL::Exception if one of the parameters if not valid
   * provided variant type is not compatible with the parameter.
   */
  void
  validateParameters(std::unordered_map<std::string, PrimitiveParameter> const
                         &tempParameters) const override;

  void
  onValidatedAndSetParameter(const std::string        &name,
                             const PrimitiveParameter &parameter) override;

private:
  // prefer composition over inheritance to avoid inheriting x/y/zExtents
  // methods that could allow to modify the cube into NOT a cube.
  Box m_box;
};

} // namespace SFCGAL

#endif // SFCGAL_CUBE_H_
