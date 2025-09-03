#ifndef SFCGAL_NURBSCURVE_H
#define SFCGAL_NURBSCURVE_H

#include "SFCGAL/BSplineCurve.h"
#include <vector>

namespace SFCGAL {

/**
 * @brief Non-Uniform Rational B-spline (NURBS) curve implementation
 *
 * NURBS curves are a generalization of B-spline curves that include rational
 * basis functions through the use of weights. They can exactly represent
 * conic sections (circles, ellipses, parabolas, hyperbolas) and provide
 * additional shape control through weight manipulation.
 *
 * A NURBS curve is defined by:
 * - A set of control points
 * - A weight for each control point
 * - A degree (order - 1)
 * - A knot vector
 *
 * Key properties:
 * - If all weights are equal, reduces to standard B-spline
 * - Can exactly represent conic sections
 * - Weights provide additional shape control
 * - Rational evaluation using homogeneous coordinates
 * - Projective invariance
 *
 * @ingroup public_api
 * @since 2.3.0
 */
class SFCGAL_API NURBSCurve : public BSplineCurve {
public:
  /**
   * @brief Default constructor creating empty curve
   */
  NURBSCurve();

  /**
   * @brief Constructor with uniform weights
   * @param controlPoints Vector of control points
   * @param degree Curve degree
   */
  NURBSCurve(const std::vector<Point> &controlPoints, unsigned int degree);

  /**
   * @brief Constructor with explicit weights
   * @param controlPoints Vector of control points
   * @param degree Curve degree
   * @param weights Vector of weights (one per control point)
   * @throws Exception if weights size doesn't match control points or contains
   * non-positive values
   */
  NURBSCurve(const std::vector<Point> &controlPoints, unsigned int degree,
             const std::vector<double> &weights);

  /**
   * @brief Constructor with weights and knot vector
   * @param controlPoints Vector of control points
   * @param degree Curve degree
   * @param weights Vector of weights
   * @param knotVector Knot vector
   * @throws Exception if parameters are incompatible
   */
  NURBSCurve(const std::vector<Point> &controlPoints, unsigned int degree,
             const std::vector<double> &weights,
             const std::vector<double> &knotVector);

  /**
   * @brief Copy constructor
   * @param other Source curve to copy
   */
  NURBSCurve(const NURBSCurve &other) = default;

  /**
   * @brief Assignment operator
   * @param other Source curve to copy
   * @return Reference to this curve
   */
  auto
  operator=(const NURBSCurve &other) -> NURBSCurve & = default;

  /**
   * @brief Virtual destructor
   */
  ~NURBSCurve() override = default;

  // Visitor pattern implementation
  void
  accept(GeometryVisitor &visitor) override;

  void
  accept(ConstGeometryVisitor &visitor) const override;

  // Geometry interface

  /**
   * @brief Get geometry type ID
   * @return TYPE_NURBSCURVE
   */
  [[nodiscard]] auto
  geometryTypeId() const -> GeometryType override;

  /**
   * @brief Get geometry type as string
   * @return "NURBSCurve"
   */
  [[nodiscard]] auto
  geometryType() const -> std::string override;

  /**
   * @brief Clone the curve
   * @return New NURBSCurve instance
   */
  [[nodiscard]] auto
  clone() const -> NURBSCurve * override;

  // Curve interface

  /**
   * @brief Evaluate curve at parameter t
   *
   * Uses rational evaluation with homogeneous coordinates when weights
   * are non-uniform, otherwise falls back to B-spline evaluation.
   *
   * @param parameter Parameter value within parameter bounds
   * @return Point on curve at parameter t
   * @throws Exception if t outside parameter bounds or curve is empty
   */
  [[nodiscard]] auto
  evaluate(double parameter) const -> Point override;

  /**
   * @brief Compute derivative at parameter t
   *
   * For rational NURBS, uses the quotient rule for derivatives.
   * For non-rational cases, uses B-spline derivative algorithm.
   *
   * @param parameter Parameter value within parameter bounds
   * @param order Derivative order (1, 2, 3, ...)
   * @return Derivative vector at parameter t
   * @throws Exception if t outside parameter bounds
   */
  [[nodiscard]] auto
  derivative(double parameter, unsigned int order = 1) const -> Point override;

  /**
   * @brief Validate curve
   * @return true if curve and weights are valid
   */
  [[nodiscard]] auto
  isValid() const -> bool override;

  // NURBS specific methods

  /**
   * @brief Get all weights
   * @return Const reference to weights vector
   */
  [[nodiscard]] auto
  weights() const -> const std::vector<double> &;

  /**
   * @brief Set all weights
   * @param weights New weights vector
   * @throws Exception if weights size doesn't match control points or contains
   * non-positive values
   */
  void
  setWeights(const std::vector<double> &weights);

  /**
   * @brief Get weight at index
   * @param index Weight index
   * @return Weight value
   * @throws Exception if index out of bounds
   */
  [[nodiscard]] auto
  weight(size_t index) const -> double;

  /**
   * @brief Set weight at index
   * @param index Weight index
   * @param weight New weight value (must be positive)
   * @throws Exception if index out of bounds or weight <= 0
   */
  void
  setWeight(size_t index, double weight);

  /**
   * @brief Check if curve is rational
   *
   * A NURBS curve is rational if not all weights are equal.
   * Non-rational NURBS curves are equivalent to B-splines.
   *
   * @return true if weights are non-uniform
   */
  [[nodiscard]] auto
  isRational() const -> bool;

  /**
   * @brief Get control points in homogeneous coordinates
   *
   * Returns control points multiplied by their weights, with weight
   * as an additional coordinate. Useful for certain algorithms.
   *
   * @return Vector of homogeneous coordinate arrays [wx, wy, wz, w]
   */
  [[nodiscard]] auto
  homogeneousControlPoints() const -> std::vector<std::vector<double>>;

  // Knot manipulation

  /**
   * @brief Insert a knot into the curve
   *
   * Adds a knot to the knot vector and computes new control points
   * and weights such that the curve shape remains unchanged.
   *
   * @param parameter Knot value to insert (must be within parameter bounds)
   * @throws Exception if t outside parameter bounds
   */
  void
  insertKnot(double parameter);

  /**
   * @brief Insert multiple knots
   * @param newKnots Vector of knot values to insert
   */
  void
  refineKnots(const std::vector<double> &newKnots);

  // Factory methods for common NURBS curves

  /**
   * @brief Create a circular arc using NURBS representation
   *
   * Creates a NURBS curve that exactly represents a circular arc.
   * This demonstrates the power of NURBS for representing conic sections.
   *
   * @param center Center point of circle
   * @param radius Radius of circle
   * @param startAngle Start angle in radians
   * @param endAngle End angle in radians
   * @param is3D Whether to create 2D or 3D arc
   * @return Unique pointer to NURBS curve representing the arc
   * @throws Exception if radius <= 0 or angles are invalid
   */
  static auto
  createCircularArc(const Point &center, double radius, double startAngle,
                    double endAngle, bool is3D = false)
      -> std::unique_ptr<NURBSCurve>;

private:
  std::vector<double> _weights; ///< Weight for each control point

  /**
   * @brief Evaluate curve using rational B-spline algorithm
   * @param parameter Parameter value
   * @return Point on curve
   */
  [[nodiscard]] auto
  evaluateRational(double parameter) const -> Point;

  /**
   * @brief Compute first derivative for rational NURBS
   * @param parameter Parameter value
   * @return First derivative vector
   */
  [[nodiscard]] auto
  computeRationalDerivative(double parameter) const -> Point;

  /**
   * @brief Compute basis functions for given span and parameter
   * @param span Knot span index
   * @param parameter Parameter value
   * @param basisFunctions Output vector for basis function values
   */
  void
  computeBasisFunctions(size_t span, double parameter,
                        std::vector<double> &basisFunctions) const;

  /**
   * @brief Compute basis function derivatives
   * @param span Knot span index
   * @param parameter Parameter value
   * @param order Derivative order
   * @param derivativeBasisFunctions Output vector for basis function
   * derivatives
   */
  void
  computeBasisDerivatives(size_t span, double parameter, unsigned int order,
                          std::vector<double> &derivativeBasisFunctions) const;

  /**
   * @brief Validate weights vector
   * @return true if weights are valid (correct size, all positive)
   */
  [[nodiscard]] auto
  validateWeights() const -> bool;

  /**
   * @brief Linear interpolation between two points using CoordinateType
   * @param point1 First point
   * @param point2 Second point
   * @param parameter Interpolation parameter [0,1]
   * @return Interpolated point
   */
  [[nodiscard]] static auto
  lerp(const Point &point1, const Point &point2, double parameter) -> Point;
};

} // namespace SFCGAL

#endif // SFCGAL_NURBSCURVE_H
