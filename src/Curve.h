// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_CURVE_H
#define SFCGAL_CURVE_H

#include "SFCGAL/Exception.h"
#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

namespace SFCGAL {

/**
 * @brief Abstract base class for parametric curves
 *
 * This class provides a common interface for parametric curves.
 * All curves are represented as NURBS (Non-Uniform Rational B-Splines)
 * which can degenerate to B-splines (uniform weights) or Bezier curves.
 *
 * @since 2.3.0
 */
class SFCGAL_API Curve : public Geometry {
public:
  using FT        = Kernel::FT; ///< Floating-point type for curve calculations
  using Parameter = FT;         ///< Parameter type for curve parameterization

  /**
   * @brief Curve continuity types
   */
  enum class Continuity : std::uint8_t {
    C0, ///< Position continuity (sharp corners allowed)
    C1, ///< Tangent continuity
    C2, ///< Curvature continuity
    G1, ///< Geometric tangent continuity
    G2  ///< Geometric curvature continuity
  };

  /**
   * @brief Default constructor
   */
  Curve();

  /**
   * @brief Copy constructor
   * @param other Curve to copy from
   */
  Curve(const Curve &other);

  /**
   * @brief Assignment operator
   * @param other Curve to assign from
   * @return Reference to this curve
   */
  auto
  operator=(const Curve &other) -> Curve &;

  /**
   * @brief Virtual destructor
   */
  ~Curve() override;

  // Geometry interface implementation

  /**
   * Return the topological dimension of this geometry.
   *
   * For Curve objects this is always 1.
   *
   * @return The integer dimension (always 1 for curves).
   */
  [[nodiscard]] auto
  dimension() const -> int override
  {
    return 1;
  }

  /**
   * Return the coordinate dimension of the curve.
   *
   * Computes 2 for non-3D curves or 3 for 3D curves, and adds 1 if the curve is
   * measured.
   *
   * @return Coordinate dimension (2 or 3, plus 1 when measured).
   */
  [[nodiscard]] auto
  coordinateDimension() const -> int override
  {
    return (is3D() ? 3 : 2) + (isMeasured() ? 1 : 0);
  }

  // Pure virtual methods for curves

  /**
   * @brief Evaluate the curve at parameter value
   * @param parameter Parameter value (usually in [tMin, tMax])
   * @return Point on the curve at parameter
   * @throw Exception if parameter is outside valid range
   */
  [[nodiscard]] virtual auto
  evaluate(Parameter parameter) const -> Point = 0;

  /**
   * @brief Compute the derivative of the curve at parameter value
   * @param parameter Parameter value
   * @param order Derivative order (1 for first derivative, 2 for second, etc.)
   * @return Vector representing the derivative at parameter
   */
  [[nodiscard]] virtual auto
  derivative(Parameter parameter, unsigned int order = 1) const -> Point = 0;

  /**
   * @brief Compute the tangent unit vector at parameter value
   * @param parameter Parameter value
   * @return Normalized tangent vector
   */
  [[nodiscard]] virtual auto
  tangent(Parameter parameter) const -> Point = 0;

  /**
   * @brief Compute the normal vector at parameter value (for planar curves)
   * @param parameter Parameter value
   * @return Normal vector perpendicular to tangent
   */
  [[nodiscard]] virtual auto
  normal(Parameter parameter) const -> Point = 0;

  /**
   * @brief Compute the binormal vector at parameter value (for 3D curves)
   * @param parameter Parameter value
   * @return Binormal vector (tangent × normal)
   */
  [[nodiscard]] virtual auto
  binormal(Parameter parameter) const -> Point = 0;

  /**
   * @brief Compute the curvature at parameter value
   * @param parameter Parameter value
   * @return Curvature value (1/radius of curvature)
   */
  [[nodiscard]] virtual auto
  curvature(Parameter parameter) const -> FT = 0;

  /**
   * @brief Compute the torsion at parameter value (for 3D curves)
   * @param parameter Parameter value
   * @return Torsion value (measure of curve twisting)
   */
  [[nodiscard]] virtual auto
  torsion(Parameter parameter) const -> FT = 0;

  /**
   * @brief Get the Frenet-Serret frame at parameter value
   * @param parameter Parameter value
   * @return Tuple of (tangent, normal, binormal) vectors
   */
  [[nodiscard]] virtual auto
  frenetFrame(Parameter parameter) const -> std::tuple<Point, Point, Point> = 0;

  /**
   * @brief Convert curve to a LineString approximation
   * @param numSegments Number of line segments to use
   * @return LineString approximating the curve
   */
  [[nodiscard]] virtual auto
  toLineString(unsigned int numSegments = 32) const
      -> std::unique_ptr<LineString> = 0;

  /**
   * @brief Convert curve to LineString with adaptive sampling
   * @param tolerance Maximum deviation from actual curve
   * @param minSegments Minimum number of segments
   * @param maxSegments Maximum number of segments
   * @return Adaptively sampled LineString
   */
  [[nodiscard]] virtual auto
  toLineStringAdaptive(FT tolerance = FT(1e-3), unsigned int minSegments = 8,
                       unsigned int maxSegments = 256) const
      -> std::unique_ptr<LineString> = 0;

  /**
   * @brief Get the parameter bounds of the curve
   * @return Pair containing (minParameter, maxParameter)
   */
  [[nodiscard]] virtual auto
  parameterBounds() const -> std::pair<Parameter, Parameter> = 0;

  /**
   * @brief Check if the curve is periodic
   * @return true if curve has periodic parameterization
   */
  [[nodiscard]] virtual auto
  isPeriodic() const -> bool = 0;

  /**
   * @brief Check if the curve is planar
   * @param plane Output parameter for the plane equation if planar
   * @return true if curve lies in a plane
   */
  [[nodiscard]] virtual auto
  isPlanar(std::vector<FT> *plane = nullptr) const -> bool = 0;

  /**
   * @brief Check if the curve is linear (straight line)
   * @return true if curve is a straight line segment
   */
  [[nodiscard]] virtual auto
  isLinear() const -> bool = 0;

  /**
   * @brief Get the degree of the curve
   * @return Polynomial degree
   */
  [[nodiscard]] virtual auto
  degree() const -> unsigned int = 0;

  /**
   * @brief Get the curve length
   * @param from Start parameter (default: curve start)
   * @param to End parameter (default: curve end)
   * @param tolerance Integration tolerance
   * @return Arc length between parameters
   */
  [[nodiscard]] virtual auto
  length(Parameter from = Parameter(-1), Parameter to = Parameter(-1),
         FT tolerance = FT(1e-6)) const -> FT = 0;

  /**
   * @brief Find parameter value at given arc length from start
   * @param arcLength Target arc length from curve start
   * @param tolerance Tolerance for root finding
   * @return Parameter value at specified arc length
   */
  [[nodiscard]] virtual auto
  parameterAtLength(FT arcLength, FT tolerance = FT(1e-6)) const
      -> Parameter = 0;

  /**
   * @brief Reparameterize curve by arc length
   * @return New curve parameterized by arc length [0, totalLength]
   */
  [[nodiscard]] virtual auto
  reparameterizeByArcLength() const -> std::unique_ptr<Curve> = 0;

  /**
   * @brief Split curve at parameter value
   * @param parameter Split location
   * @return Pair of curves (before split, after split)
   */
  [[nodiscard]] virtual auto
  split(Parameter parameter) const
      -> std::pair<std::unique_ptr<Curve>, std::unique_ptr<Curve>> = 0;

  /**
   * @brief Extract subcurve between two parameters
   * @param from Start parameter
   * @param to End parameter
   * @return Subcurve between parameters
   */
  [[nodiscard]] virtual auto
  subcurve(Parameter from, Parameter to) const -> std::unique_ptr<Curve> = 0;

  /**
   * @brief Reverse the curve direction
   * @return Reversed curve
   */
  [[nodiscard]] virtual auto
  reverse() const -> std::unique_ptr<Curve> = 0;

  /**
   * @brief Join with another curve
   * @param other Curve to join with
   * @param continuity Desired continuity at join
   * @param tolerance Tolerance for join
   * @return Joined curve or nullptr if join impossible
   */
  [[nodiscard]] virtual auto
  join(const Curve &other, Continuity continuity = Continuity::C0,
       FT tolerance = FT(1e-6)) const -> std::unique_ptr<Curve> = 0;

  /**
   * @brief Project point onto curve
   * @param point Point to project
   * @param parameter Output parameter at closest point
   * @return Closest point on curve
   */
  [[nodiscard]] virtual auto
  closestPoint(const Point &point, Parameter *parameter = nullptr) const
      -> Point = 0;

  // NURBS-specific interface methods (virtual with default implementations)

  /**
   * @brief Get number of control points
   * @return Control point count
   * @throws Exception if not implemented for curve type
   */
  [[nodiscard]] virtual auto
  numControlPoints() const -> size_t
  {
    BOOST_THROW_EXCEPTION(
        Exception("numControlPoints() not supported for this curve type"));
  }

  /**
   * @brief Check if curve uses non-uniform rational weights
   * @return false for non-NURBS curves, true if NURBS with varying weights
   * @throws Exception if not implemented for curve type
   */
  [[nodiscard]] virtual auto
  isRational() const -> bool
  {
    BOOST_THROW_EXCEPTION(
        Exception("isRational() not supported for this curve type"));
  }

  /**
   * Return the rational weight of the curve's control point at the given index.
   *
   * Default implementation throws an Exception indicating that weight access
   * is not supported for this curve type; concrete NURBS-derived classes must
   * override to provide valid weights.
   *
   * @param index Index of the control point weight.
   * @return Weight value associated with the control point.
   * @throws Exception if this curve type does not support weights or if the
   *         index is out of range.
   */
  [[nodiscard]] virtual auto
  weight(size_t index) const -> FT
  {
    static_cast<void>(index); // Suppress unused parameter warning
    BOOST_THROW_EXCEPTION(
        Exception("weight() not supported for this curve type"));
  }

  /**
   * Return a read-only reference to the control point at the given index.
   *
   * This method is part of the NURBS-specific interface. The base
   * implementation always throws an Exception to indicate the operation is not
   * supported for curve types that do not expose control points.
   *
   * @param index Index of the control point.
   * @return Const reference to the control point at `index`.
   * @throws Exception if the curve type does not support control points or if
   *         `index` is out of range for implementations that override this
   * method.
   */
  [[nodiscard]] virtual auto
  controlPointN(size_t index) const -> const Point &
  {
    static_cast<void>(index); // Suppress unused parameter warning
    BOOST_THROW_EXCEPTION(
        Exception("controlPointN() not supported for this curve type"));
  }

  /**
   * Return a modifiable reference to the control point at the given index.
   *
   * Concrete NURBS-style curve implementations should override this to provide
   * writable access to their control points. The base implementation throws an
   * Exception indicating the operation is not supported for this curve type.
   *
   * @param index Zero-based control point index.
   * @return Reference to the control point that can be modified.
   * @throws Exception If the curve type does not support control point access
   *         or if index is out of range.
   */
  [[nodiscard]] virtual auto
  controlPointN(size_t index) -> Point &
  {
    static_cast<void>(index); // Suppress unused parameter warning
    BOOST_THROW_EXCEPTION(
        Exception("controlPointN() not supported for this curve type"));
  }

  // Convenience methods (implemented in base class)

  /**
   * @brief Get the start point of the curve
   * @return First point of the curve
   */
  [[nodiscard]] auto
  startPoint() const -> Point;

  /**
   * @brief Get the end point of the curve
   * @return Last point of the curve
   */
  [[nodiscard]] auto
  endPoint() const -> Point;

  /**
   * @brief Get the midpoint of the curve
   * @return Point at middle parameter value
   */
  [[nodiscard]] auto
  midPoint() const -> Point;

  /**
   * @brief Check if a parameter value is within curve bounds
   * @param parameter Parameter to check
   * @return true if parameter is valid for this curve
   */
  [[nodiscard]] auto
  isValidParameter(const Parameter &parameter) const -> bool;

  /**
   * @brief Normalize parameter to [0,1] range
   * @param parameter Raw parameter value
   * @return Normalized parameter in [0,1]
   */
  [[nodiscard]] auto
  normalizeParameter(const Parameter &parameter) const -> FT;

  /**
   * @brief Denormalize parameter from [0,1] to actual range
   * @param normalizedParameter Parameter in [0,1]
   * @return Actual parameter value
   */
  [[nodiscard]] auto
  denormalizeParameter(const FT &normalizedParameter) const -> Parameter;
};

} // namespace SFCGAL

#endif // SFCGAL_CURVE_H
