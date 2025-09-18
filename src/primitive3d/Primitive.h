// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_PRIMITIVE_H_
#define SFCGAL_PRIMITIVE_H_

#include <string>
#include <unordered_map>
#include <variant>

#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/export.h"

namespace SFCGAL {

/**
 * @brief Enumeration of available primitive types.
 */
enum class PrimitiveType : std::int8_t {
  TYPE_CYLINDER = 0,
  TYPE_SPHERE   = 1,
  TYPE_TORUS    = 2,
  TYPE_BOX      = 3,
  TYPE_CUBE     = 4,
  TYPE_CONE     = 5
};

/**
 * @brief Variant type representing the possible values of a primitive
 * parameter.
 *
 * A primitive parameter can be one of the following types:
 * - SFCGAL::Kernel::FT        : a scalar floating-point value
 * - SFCGAL::Kernel::Point_3   : a 3D point
 * - SFCGAL::Kernel::Vector_3  : a 3D vector
 * - unsigned int              : an integer value
 */
using PrimitiveParameter =
    std::variant<Kernel::FT, Kernel::Point_3, Kernel::Vector_3, unsigned int>;

/**
 * @brief Abstract base class for primitive object (Cube, Cylinder, Cone, etc.)
 *
 * This class provides generic parameter functions to set, get, list member
 * fields. All parameters needed by inherited classes must be declared in the
 * constructor:
 *
 * \code{.cpp}
 * SubClass::SubClass(const Kernel::FT &radius, const Kernel::FT &height,
 *                    unsigned int num)
 * {
 *   m_parameters["radius"] = radius;
 *   m_parameters["height"] = height;
 *   m_parameters["num"]    = num;
 *
 *   SubClass::validateParameters(m_parameters);
 * }
 * \endcode
 */
class SFCGAL_API Primitive {

public:
  /**
   * @brief Default constructor.
   */
  Primitive() = default;

  /**
   * @brief Assignment operator
   */
  auto
  operator=(const Primitive &other) -> Primitive &;

  /**
   * @brief Destructor.
   */
  virtual ~Primitive() = default;

  /**
   * @brief returns the primitive type
   * @warning use CamelCase (Cylinder, not CYLINDER)
   */
  [[nodiscard]] virtual auto
  primitiveType() const -> std::string = 0;

  /**
   * @brief Returns a code corresponding to the type
   */
  [[nodiscard]] virtual auto
  primitiveTypeId() const -> PrimitiveType = 0;

  /**
   * @brief Generates a surface mesh representation of the primitive
   * @return A PolyhedralSurface object representing the torus
   */
  virtual auto
  generatePolyhedralSurface() const -> PolyhedralSurface = 0;

  /**
   * @brief Sets the value of a primitive parameter.
   * Assigns a new value to the parameter identified by @p name.
   * @param name The name of the parameter to set.
   * @param parameter The new parameter value as a ::PrimitiveParameter.
   * @pre The parameter identified by @p name must exist and accept the provided
   * variant type.
   * @throws SFCGAL::Exception if the parameter does not exist or if the
   * provided variant type is not compatible with the parameter.
   */
  void
  setParameter(const std::string &name, const PrimitiveParameter &parameter);

  /**
   * @brief Retrieves the value of a primitive parameter based ont its @p name.
   * @param name The name of the parameter to retrieve.
   * @return The parameter value as a ::PrimitiveParameter variant.
   * @pre The parameter identified by @p name must exist.
   * @throws SFCGAL::Exception if the parameter does not exist or if the
   * provided variant type is not compatible with the parameter.
   */
  [[nodiscard]] auto
  parameter(const std::string &name) const -> PrimitiveParameter;

  /**
   * @brief Retrieves the list of parameters
   */
  [[nodiscard]] auto
  parameters() const -> std::unordered_map<std::string, PrimitiveParameter>;

  /**
   * Equality operator
   * @param other other geometry to compare to
   * @param epsilon allowed tolerance
   */
  [[nodiscard]] auto
  almostEqual(const Primitive &other, double epsilon) const -> bool;

  /**
   * Returns string representation of this object. Can be override by
   * sub-classes.
   *
   * \return string representation of this object
   */
  [[nodiscard]] virtual auto
  toString() const -> std::string;

  /**
   * overloaded operator<<
   * \param out stream to fill
   * \param prim primitive to print via Primitive::toString()
   * \return filled stream
   */
  friend std::ostream &
  operator<<(std::ostream &out, const Primitive &prim)
  {
    // Call virtual function toString() to get the string to be printed
    out << prim.toString();
    return out;
  }

protected:
  /**
   * @brief Invalidates the cached geometries
   */
  virtual void
  invalidateCache();

  /**
   * @brief Verifies that all parameters are valid. For instance, it raises an
   * error if a radius is negative.
   * @param tempParameters a temp map of parameter with new values
   * @throws SFCGAL::Exception if one of the parameters if not valid
   * provided variant type is not compatible with the parameter.
   */
  virtual void
  validateParameters(std::unordered_map<std::string, PrimitiveParameter> const
                         &tempParameters) const = 0;

  /**
   * @brief Checks that the new value is valid and sets it.
   * This also invalidates the cache.
   * Unlike @ref setParameter, this does not check that name exists and that it
   * accepts the provided variant type. This assumes that name exists and the
   * variant type is compatible with the parameter
   * @param name parameter name to check
   * @param parameter value
   * @throws SFCGAL::Exception if the parameter is not valid
   */
  virtual void
  validateAndSetParameter(const std::string        &name,
                          const PrimitiveParameter &parameter);

  /**
   * @brief Callback to propagate paramater change after a parameter has been
   * validated and set.
   * @param name parameter name to propagate
   * @param parameter value
   * @throws SFCGAL::Exception if propagation failed
   */
  virtual void
  onValidatedAndSetParameter(const std::string        &name,
                             const PrimitiveParameter &parameter);

  std::unordered_map<std::string, PrimitiveParameter> m_parameters;
  mutable std::optional<PolyhedralSurface>            m_polyhedral_surface;
};

/**
 * Equality operator
 * @todo only compare coordinate points
 * @pre the two geometries must be valid
 */
SFCGAL_API auto
operator==(const Primitive &, const Primitive &) -> bool;

} // namespace SFCGAL

#endif // SFCGAL_PRIMITIVE_H_
