#ifndef SFCGAL_HERMITECURVE_H
#define SFCGAL_HERMITECURVE_H

#include "SFCGAL/Curve.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include <memory>
#include <utility>
#include <vector>

namespace SFCGAL {

/**
 * Hermite spline curve with explicit control over tangent handles
 * Provides smooth interpolation through control points with user-defined
 * tangents
 */
class SFCGAL_API HermiteCurve : public Curve {
public:
  struct ControlPoint {
    Point position;
    Point inHandle;
    Point outHandle;

    ControlPoint() = default;
    ControlPoint(const Point &pos)
        : position(pos), inHandle(pos), outHandle(pos)
    {
    }
    ControlPoint(const Point &pos, const Point &inHandle,
                 const Point &outHandle)
        : position(pos), inHandle(inHandle), outHandle(outHandle)
    {
    }
  };

  HermiteCurve();
  HermiteCurve(const std::vector<ControlPoint> &controlPoints);
  HermiteCurve(const HermiteCurve &other);
  auto
  operator=(const HermiteCurve &other) -> HermiteCurve &;
  ~HermiteCurve() override = default;

  // Geometry interface
  [[nodiscard]] auto
  geometryTypeId() const -> GeometryType override;
  [[nodiscard]] auto
  geometryType() const -> std::string override;
  [[nodiscard]] auto
  clone() const -> HermiteCurve * override;
  [[nodiscard]] auto
  isEmpty() const -> bool override;
  [[nodiscard]] auto
  is3D() const -> bool override;
  [[nodiscard]] auto
  isMeasured() const -> bool override;

  // Curve interface
  [[nodiscard]] auto
  evaluate(double parameter) const -> Point override;
  [[nodiscard]] auto
  derivative(double parameter, unsigned int order = 1) const -> Point override;
  [[nodiscard]] auto
  toLineString(unsigned int numSegments = 32) const
      -> std::unique_ptr<LineString> override;
  [[nodiscard]] auto
  parameterBounds() const -> std::pair<double, double> override;
  [[nodiscard]] auto
  isValid() const -> bool override;

  // HermiteCurve specific
  [[nodiscard]] auto
  numControlPoints() const -> size_t override
  {
    return _controlPoints.size();
  }
  [[nodiscard]] auto
  getControlPoint(size_t index) const -> const ControlPoint &;
  auto
  getControlPoint(size_t index) -> ControlPoint &;
  void
  setControlPoint(size_t index, const ControlPoint &controlPoint);
  void
  addControlPoint(const ControlPoint &controlPoint);

  // Handle manipulation
  void
  setInHandle(size_t index, const Point &handle);
  void
  setOutHandle(size_t index, const Point &handle);
  [[nodiscard]] auto
  getInHandle(size_t index) const -> const Point &;
  [[nodiscard]] auto
  getOutHandle(size_t index) const -> const Point &;

  // Auto-smooth handles
  void
  autoSmoothHandles(size_t index);
  void
  autoSmoothAllHandles();

  // Convert to/from Bezier
  [[nodiscard]] auto
  toBezierCurve() const -> std::unique_ptr<BezierCurve>;
  static auto
  fromBezierCurve(const BezierCurve &bezier) -> std::unique_ptr<HermiteCurve>;

private:
  std::vector<ControlPoint> _controlPoints;

  /**
   * @brief Performs Hermite interpolation between two control points
   * @param point1 First control point with position and tangent handles
   * @param point2 Second control point with position and tangent handles
   * @param parameter Local parameter value within the segment [0,1]
   * @return Interpolated point at the given parameter
   *
   * Uses cubic Hermite basis functions to smoothly interpolate between
   * control points while respecting their tangent handle constraints.
   */
  [[nodiscard]] auto
  hermiteInterpolate(const ControlPoint &point1, const ControlPoint &point2,
                     double parameter) const -> Point;

  /**
   * @brief Computes the tangent vector along a Hermite curve segment
   * @param point1 First control point with position and tangent handles
   * @param point2 Second control point with position and tangent handles
   * @param parameter Local parameter value within the segment [0,1]
   * @return Tangent vector at the given parameter
   *
   * Calculates the derivative of the Hermite curve using the derivatives
   * of the cubic basis functions.
   */
  [[nodiscard]] auto
  hermiteTangent(const ControlPoint &point1, const ControlPoint &point2,
                 double parameter) const -> Point;

  /**
   * @brief Calculates smooth handles for the first control point in the curve
   * @param index Index of the first control point (must be 0)
   *
   * Computes appropriate in and out handles for the curve's starting point
   * based on the direction to the next control point.
   */
  void
  calculateFirstPointHandles(size_t index);

  /**
   * @brief Calculates smooth handles for the last control point in the curve
   * @param index Index of the last control point
   *
   * Computes appropriate in and out handles for the curve's ending point
   * based on the direction from the previous control point.
   */
  void
  calculateLastPointHandles(size_t index);

  /**
   * @brief Calculates distances to neighboring control points for handle
   * computation
   * @param index Index of the middle control point
   * @return Pair containing distances to previous and next control points
   *
   * Helper function that computes the geometric distances to adjacent control
   * points, used for proportional handle length calculations.
   */
  [[nodiscard]] auto
  calculateNeighborDistances(size_t index) const -> std::pair<double, double>;

  /**
   * @brief Calculates smooth handles for middle control points in the curve
   * @param index Index of the middle control point
   *
   * Computes handles for interior control points by averaging directions
   * to neighboring points and scaling by relative distances for smooth
   * transitions.
   */
  void
  calculateMiddlePointHandles(size_t index);
};

} // namespace SFCGAL

#endif // SFCGAL_HERMITECURVE_H
