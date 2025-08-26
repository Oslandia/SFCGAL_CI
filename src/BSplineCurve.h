#ifndef SFCGAL_BSPLINECURVE_H
#define SFCGAL_BSPLINECURVE_H

#include "SFCGAL/Curve.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include <memory>
#include <vector>

namespace SFCGAL {

/**
 * @brief B-spline curve implementation
 *
 * A B-spline (Basis spline) curve is a generalization of Bezier curves defined
 * by:
 * - A set of control points
 * - A degree (order - 1)
 * - A knot vector that determines parameter intervals
 *
 * B-spline curves offer several advantages over Bezier curves:
 * - Local control: moving a control point affects only nearby portions of the
 * curve
 * - Efficient evaluation using De Boor's algorithm
 * - Support for non-uniform parameterization through knot vectors
 * - Can represent both curves that pass through control points (interpolating)
 *   and curves that approximate control points
 *
 * Properties:
 * - Degree can be less than (number of control points - 1)
 * - Generally do not pass through control points
 * - Lie within convex hull of relevant control points
 * - C^(degree-1) continuity at knots (except where knots have multiplicity)
 *
 * @ingroup public_api
 * @since 2.3.0
 */
class SFCGAL_API BSplineCurve : public Curve {
public:
  /**
   * @brief Default constructor creating empty curve
   */
  BSplineCurve();

  /**
   * @brief Constructor with control points and degree
   *
   * Creates a B-spline with uniform knot vector. The degree must be less than
   * the number of control points.
   *
   * @param controlPoints Vector of control points
   * @param degree Curve degree (must be < controlPoints.size())
   * @throws Exception if degree >= number of control points
   */
  BSplineCurve(const std::vector<Point> &controlPoints, unsigned int degree);

  /**
   * @brief Constructor with control points, degree, and knot vector
   *
   * @param controlPoints Vector of control points
   * @param degree Curve degree
   * @param knotVector Knot vector (size must equal controlPoints.size() +
   * degree + 1)
   * @throws Exception if parameters are incompatible
   */
  BSplineCurve(const std::vector<Point> &controlPoints, unsigned int degree,
               const std::vector<double> &knotVector);

  /**
   * @brief Copy constructor
   * @param other Source curve to copy
   */
  BSplineCurve(const BSplineCurve &other);

  /**
   * @brief Assignment operator
   * @param other Source curve to copy
   * @return Reference to this curve
   */
  auto
  operator=(const BSplineCurve &other) -> BSplineCurve &;

  /**
   * @brief Virtual destructor
   */
  ~BSplineCurve() override = default;

  // Geometry interface

  /**
   * @brief Get geometry type ID
   * @return TYPE_BSPLINECURVE
   */
  [[nodiscard]] auto
  geometryTypeId() const -> GeometryType override;

  /**
   * @brief Get geometry type as string
   * @return "BSplineCurve"
   */
  [[nodiscard]] auto
  geometryType() const -> std::string override;

  /**
   * @brief Clone the curve
   * @return New BSplineCurve instance
   */
  [[nodiscard]] auto
  clone() const -> BSplineCurve * override;

  /**
   * @brief Check if curve is empty
   * @return true if no control points defined
   */
  [[nodiscard]] auto
  isEmpty() const -> bool override;

  /**
   * @brief Check if curve is 3D
   * @return true if control points have Z coordinates
   */
  [[nodiscard]] auto
  is3D() const -> bool override;

  /**
   * @brief Check if curve has measure values
   * @return true if control points have M coordinates
   */
  [[nodiscard]] auto
  isMeasured() const -> bool override;

  // Curve interface

  /**
   * @brief Evaluate curve at parameter t using De Boor's algorithm
   * @param parameter Parameter value within parameter bounds
   * @return Point on curve at parameter t
   * @throws Exception if t outside parameter bounds or curve is empty
   */
  [[nodiscard]] auto
  evaluate(double parameter) const -> Point override;

  /**
   * @brief Compute derivative at parameter t
   * @param parameter Parameter value within parameter bounds
   * @param order Derivative order (1, 2, 3, ...)
   * @return Derivative vector at parameter t
   * @throws Exception if t outside parameter bounds
   */
  [[nodiscard]] auto
  derivative(double parameter, unsigned int order = 1) const -> Point override;

  /**
   * @brief Convert to LineString approximation
   * @param numSegments Number of line segments for approximation
   * @return LineString approximating the curve
   * @throws Exception if numSegments < 2
   */
  [[nodiscard]] auto
  toLineString(unsigned int numSegments = 32) const
      -> std::unique_ptr<LineString> override;

  /**
   * @brief Get parameter bounds from knot vector
   * @return Pair containing (first_internal_knot, last_internal_knot)
   */
  [[nodiscard]] auto
  parameterBounds() const -> std::pair<double, double> override;

  /**
   * @brief Validate curve
   * @return true if curve is valid (consistent dimensions, valid knot vector,
   * etc.)
   */
  [[nodiscard]] auto
  isValid() const -> bool override;

  // BSplineCurve specific methods

  /**
   * @brief Get curve degree
   * @return Degree of the curve
   */
  [[nodiscard]] auto
  degree() const -> unsigned int;

  /**
   * @brief Get number of control points
   * @return Number of control points
   */
  [[nodiscard]] auto
  numControlPoints() const -> size_t override;

  /**
   * @brief Get control point at index (const)
   * @param index Control point index
   * @return Const reference to control point
   * @throws Exception if index out of bounds
   */
  [[nodiscard]] auto
  controlPointAt(size_t index) const -> const Point & override;

  /**
   * @brief Get control point at index (mutable)
   * @param index Control point index
   * @return Reference to control point
   * @throws Exception if index out of bounds
   */
  auto
  controlPointAt(size_t index) -> Point & override;

  /**
   * @brief Set control point at index
   * @param index Control point index
   * @param point New control point value
   * @throws Exception if index out of bounds
   */
  void
  setControlPoint(size_t index, const Point &point);

  /**
   * @brief Get all control points
   * @return Vector of control points
   */
  [[nodiscard]] auto
  controlPoints() const -> std::vector<Point> override;

  /**
   * @brief Get the knot vector
   * @return Const reference to knot vector
   */
  [[nodiscard]] auto
  knotVector() const -> const std::vector<double> &;

  /**
   * @brief Set the knot vector
   * @param knots New knot vector
   * @throws Exception if knot vector is invalid for current control points and
   * degree
   */
  void
  setKnotVector(const std::vector<double> &knots);

  /**
   * @brief Generate uniform knot vector
   *
   * Creates a uniform knot vector with appropriate multiplicities at the ends.
   * This is the default knot vector type for most applications.
   */
  void
  generateUniformKnotVector();

protected:
  std::vector<Point>  _controlPoints; ///< Control points defining the curve
  unsigned int        _degree;        ///< Degree of the curve
  std::vector<double> _knotVector; ///< Knot vector determining parameterization

  /**
   * @brief Find the knot span for parameter t
   *
   * Uses binary search to find the knot span containing parameter t.
   * This is a key step in De Boor's algorithm.
   *
   * @param parameter Parameter value
   * @return Index of knot span containing t
   */
  [[nodiscard]] auto
  findSpan(double parameter) const -> size_t;

  /**
   * @brief Validate the knot vector
   *
   * Checks that the knot vector:
   * - Has correct size (n + p + 1 where n = number of control points, p =
   * degree)
   * - Values are non-decreasing
   * - Compatible with current control points and degree
   *
   * @return true if knot vector is valid
   */
  [[nodiscard]] auto
  validateKnotVector() const -> bool;
};

} // namespace SFCGAL

#endif // SFCGAL_BSPLINECURVE_H
