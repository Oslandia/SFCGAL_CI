#ifndef SFCGAL_CURVE_H
#define SFCGAL_CURVE_H

#include "SFCGAL/Geometry.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include <memory>
#include <vector>

namespace SFCGAL {

// Forward declarations
class BezierCurve;
class BSplineCurve;
class NURBSCurve;
class HermiteCurve;
class CatmullRomCurve;

/**
 * @brief Abstract base class for parametric curves
 *
 * This class provides a common interface for all parametric curve types
 * including Bezier, B-spline, NURBS, Hermite, and Catmull-Rom curves.
 *
 * A parametric curve is defined by a parameter t that varies over a specific
 * range (usually [0,1]) and maps to points in 2D or 3D space.
 *
 * @ingroup public_api
 * @since 2.3.0
 */
class SFCGAL_API Curve : public Geometry {
public:
  /**
   * @brief Default constructor
   */
  Curve() = default;

  /**
   * @brief Virtual destructor
   */
  ~Curve() override = default;

  // Abstract Geometry interface implementation in derived classes

  /**
   * @brief Returns the dimension of the curve (always 1)
   * @return 1 for curves
   */
  [[nodiscard]] auto
  dimension() const -> int override
  {
    return 1;
  }

  /**
   * @brief Returns the coordinate dimension (2D or 3D)
   * @return 2 for 2D curves, 3 for 3D curves
   */
  [[nodiscard]] auto
  coordinateDimension() const -> int override
  {
    return is3D() ? 3 : 2;
  }

  /**
   * @brief Drop Z coordinate (not supported for curves yet)
   * @return false (not implemented)
   */
  auto
  dropZ() -> bool override
  {
    return false;
  }

  /**
   * @brief Drop M coordinate (not supported for curves yet)
   * @return false (not implemented)
   */
  auto
  dropM() -> bool override
  {
    return false;
  }

  /**
   * @brief Swap X and Y coordinates (not supported for curves yet)
   */
  void
  swapXY() override
  {
  }

  // Visitor pattern - implement in derived classes

  /**
   * @brief Accept a geometry visitor
   * @param visitor The visitor to accept
   */
  void
  accept(GeometryVisitor &visitor) override
  {
  }

  /**
   * @brief Accept a const geometry visitor
   * @param visitor The const visitor to accept
   */
  void
  accept(ConstGeometryVisitor &visitor) const override
  {
  }

  // Pure virtual methods for curves

  /**
   * @brief Evaluate the curve at parameter t
   * @param parameter Parameter value (usually in [0,1])
   * @return Point on the curve at parameter t
   * @throws Exception if t is outside valid parameter range
   */
  [[nodiscard]] virtual auto
  evaluate(double parameter) const -> Point = 0;

  /**
   * @brief Compute the derivative of the curve at parameter t
   * @param parameter Parameter value (usually in [0,1])
   * @param order Derivative order (1 for first derivative, 2 for second, etc.)
   * @return Vector representing the derivative at parameter t
   * @throws Exception if t is outside valid parameter range
   */
  [[nodiscard]] virtual auto
  derivative(double parameter, unsigned int order = 1) const -> Point = 0;

  /**
   * @brief Convert curve to a LineString approximation
   * @param numSegments Number of line segments to use for approximation
   * @return LineString approximating the curve
   */
  [[nodiscard]] virtual auto
  toLineString(unsigned int numSegments = 32) const
      -> std::unique_ptr<LineString> = 0;

  /**
   * @brief Get the parameter bounds of the curve
   * @return Pair containing (min_parameter, max_parameter)
   */
  [[nodiscard]] virtual auto
  parameterBounds() const -> std::pair<double, double> = 0;

  /**
   * @brief Check if the curve is valid
   * @return true if curve is valid, false otherwise
   */
  [[nodiscard]] virtual auto
  isValid() const -> bool = 0;

  // Control points interface - for curves that use control points

  /**
   * @brief Get the number of control points
   * @return Number of control points (0 for curves without control points)
   */
  [[nodiscard]] virtual auto
  numControlPoints() const -> size_t
  {
    return 0;
  }

  /**
   * @brief Get control point at index (const version)
   * @param index Index of control point
   * @return Const reference to control point
   * @throws Exception if index is out of bounds
   */
  [[nodiscard]] virtual auto
  controlPointAt(size_t index) const -> const Point &
  {
    static Point dummy;
    return dummy;
  }

  /**
   * @brief Get control point at index (mutable version)
   * @param index Index of control point
   * @return Reference to control point
   * @throws Exception if index is out of bounds
   */
  virtual auto
  controlPointAt(size_t index) -> Point &
  {
    static Point dummy;
    return dummy;
  }

  /**
   * @brief Get all control points
   * @return Vector of control points
   */
  [[nodiscard]] virtual auto
  controlPoints() const -> std::vector<Point>
  {
    return {};
  }

  // Convenience methods

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

  // Static factory methods for curve conversion

  /**
   * @brief Convert any curve to Bezier representation
   * @param curve Source curve to convert
   * @param degree Target degree for Bezier curve
   * @return Unique pointer to new Bezier curve
   */
  [[nodiscard]] static auto
  toBezier(const Curve &curve, unsigned int degree = 3)
      -> std::unique_ptr<BezierCurve>;

  /**
   * @brief Convert any curve to B-spline representation
   * @param curve Source curve to convert
   * @param degree Target degree for B-spline curve
   * @return Unique pointer to new B-spline curve
   */
  [[nodiscard]] static auto
  toBSpline(const Curve &curve, unsigned int degree = 3)
      -> std::unique_ptr<BSplineCurve>;

  /**
   * @brief Convert any curve to NURBS representation
   * @param curve Source curve to convert
   * @param degree Target degree for NURBS curve
   * @return Unique pointer to new NURBS curve
   */
  [[nodiscard]] static auto
  toNURBS(const Curve &curve, unsigned int degree = 3)
      -> std::unique_ptr<NURBSCurve>;

  /**
   * @brief Convert any curve to Hermite representation
   * @param curve Source curve to convert
   * @return Unique pointer to new Hermite curve
   */
  [[nodiscard]] static auto
  toHermite(const Curve &curve) -> std::unique_ptr<HermiteCurve>;

  /**
   * @brief Convert any curve to Catmull-Rom representation
   * @param curve Source curve to convert
   * @return Unique pointer to new Catmull-Rom curve
   */
  [[nodiscard]] static auto
  toCatmullRom(const Curve &curve) -> std::unique_ptr<CatmullRomCurve>;

  /**
   * @brief Fit a Bezier curve to a set of points
   * @param points Points to fit curve to
   * @param degree Degree of the Bezier curve
   * @return Unique pointer to fitted Bezier curve
   */
  [[nodiscard]] static auto
  fitBezier(const std::vector<Point> &points, unsigned int degree)
      -> std::unique_ptr<BezierCurve>;

  /**
   * @brief Fit a B-spline curve to a set of points
   * @param points Points to fit curve to
   * @param degree Degree of the B-spline curve
   * @return Unique pointer to fitted B-spline curve
   */
  [[nodiscard]] static auto
  fitBSpline(const std::vector<Point> &points, unsigned int degree)
      -> std::unique_ptr<BSplineCurve>;

  /**
   * @brief Fit a NURBS curve to a set of points
   * @param points Points to fit curve to
   * @param degree Degree of the NURBS curve
   * @return Unique pointer to fitted NURBS curve
   */
  [[nodiscard]] static auto
  fitNURBS(const std::vector<Point> &points, unsigned int degree)
      -> std::unique_ptr<NURBSCurve>;
};

} // namespace SFCGAL

#endif // SFCGAL_CURVE_H
