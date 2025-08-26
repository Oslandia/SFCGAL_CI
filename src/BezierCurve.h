#ifndef SFCGAL_BEZIERCURVE_H
#define SFCGAL_BEZIERCURVE_H

#include "SFCGAL/Curve.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include <memory>
#include <vector>

namespace SFCGAL {

/**
 * @brief Bezier curve implementation
 *
 * A Bezier curve is a parametric curve defined by a set of control points.
 * The curve passes through the first and last control points and is influenced
 * by the intermediate control points. The curve is evaluated using De
 * Casteljau's algorithm, which provides numerical stability.
 *
 * Properties:
 * - Passes through first and last control points
 * - Lies within convex hull of control points
 * - Degree = number of control points - 1
 * - C∞ continuity within each curve segment
 *
 * @ingroup public_api
 * @since 2.3.0
 */
class SFCGAL_API BezierCurve : public Curve {
public:
  /**
   * @brief Default constructor creating empty curve
   */
  BezierCurve() = default;

  /**
   * @brief Constructor with control points
   * @param controlPoints Vector of control points defining the curve
   */
  explicit BezierCurve(const std::vector<Point> &controlPoints);

  /**
   * @brief Copy constructor
   * @param other Source curve to copy
   */
  BezierCurve(const BezierCurve &other) = default;

  /**
   * @brief Assignment operator
   * @param other Source curve to copy
   * @return Reference to this curve
   */
  auto
  operator=(const BezierCurve &other) -> BezierCurve & = default;

  /**
   * @brief Virtual destructor
   */
  ~BezierCurve() override = default;

  void
  accept(GeometryVisitor &visitor) override;

  void
  accept(ConstGeometryVisitor &visitor) const override;

  // Geometry interface

  /**
   * @brief Get geometry type ID
   * @return TYPE_BEZIERCURVE
   */
  [[nodiscard]] auto
  geometryTypeId() const -> GeometryType override;

  /**
   * @brief Get geometry type as string
   * @return "BezierCurve"
   */
  [[nodiscard]] auto
  geometryType() const -> std::string override;

  /**
   * @brief Clone the curve
   * @return New BezierCurve instance
   */
  [[nodiscard]] auto
  clone() const -> BezierCurve * override;

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
   * @brief Evaluate curve at parameter t using De Casteljau's algorithm
   * @param parameter Parameter value in [0,1]
   * @return Point on curve at parameter t
   * @throws Exception if t outside [0,1] or curve is empty
   */
  [[nodiscard]] auto
  evaluate(double parameter) const -> Point override;

  /**
   * @brief Compute derivative at parameter t
   * @param parameter Parameter value in [0,1]
   * @param order Derivative order (1, 2, 3, ...)
   * @return Derivative vector at parameter t
   * @throws Exception if t outside [0,1]
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
   * @brief Get parameter bounds
   * @return Always returns (0.0, 1.0) for Bezier curves
   */
  [[nodiscard]] auto
  parameterBounds() const -> std::pair<double, double> override;

  /**
   * @brief Validate curve
   * @return true if all control points have consistent dimensions
   */
  [[nodiscard]] auto
  isValid() const -> bool override;

  // BezierCurve specific methods

  /**
   * @brief Get curve degree
   * @return Degree of curve (number of control points - 1)
   */
  [[nodiscard]] auto
  degree() const -> unsigned int;

  /**
   * @brief Get number of control points
   * @return Number of control points
   */
  [[nodiscard]] auto
  numControlPoints() const -> size_t override
  {
    return _controlPoints.size();
  }

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
   * @brief Add control point to end
   * @param point Control point to add
   * @throws Exception if point dimensions inconsistent with existing points
   */
  void
  addControlPoint(const Point &point);

  /**
   * @brief Remove control point at index
   * @param index Index of control point to remove
   * @throws Exception if index out of bounds
   */
  void
  removeControlPoint(size_t index);

  /**
   * @brief Clear all control points
   */
  void
  clear();

  /**
   * @brief Evaluate curve showing intermediate steps of De Casteljau's
   * algorithm
   *
   * This method returns all intermediate point arrays computed during
   * De Casteljau's algorithm evaluation. Useful for visualization and
   * debugging.
   *
   * @param parameter Parameter value in [0,1]
   * @return Vector of point arrays, one for each level of the algorithm
   */
  [[nodiscard]] auto
  evaluateWithSteps(double parameter) const -> std::vector<std::vector<Point>>;

private:
  std::vector<Point> _controlPoints; ///< Control points defining the curve

  /**
   * @brief Linear interpolation between two points
   * @param point1 First point
   * @param point2 Second point
   * @param parameter Interpolation parameter [0,1]
   * @return Interpolated point
   */
  [[nodiscard]] static auto
  lerp(const Point &point1, const Point &point2, double parameter) -> Point;
};

} // namespace SFCGAL

#endif // SFCGAL_BEZIERCURVE_H
