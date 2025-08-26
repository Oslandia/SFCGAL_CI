#ifndef SFCGAL_CATMULLROMCURVE_H
#define SFCGAL_CATMULLROMCURVE_H

#include "SFCGAL/Curve.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include <memory>
#include <vector>

namespace SFCGAL {

/**
 * Catmull-Rom spline - interpolating curve passing through all control points
 * Provides C1 continuity with automatic tangent calculation
 */
class SFCGAL_API CatmullRomCurve : public Curve {
public:
  enum BoundaryType : std::uint8_t {
    BOUNDARY_ZERO,    // Zero tangent at boundaries
    BOUNDARY_CLAMPED, // Use first/last segment direction
    BOUNDARY_CYCLIC,  // Closed curve
    BOUNDARY_NATURAL  // Natural spline boundary conditions
  };

  CatmullRomCurve();
  CatmullRomCurve(const std::vector<Point> &points, double tension = 0.5);
  CatmullRomCurve(const LineString &lineString, double tension = 0.5);
  CatmullRomCurve(const CatmullRomCurve &other) = default;
  auto
  operator=(const CatmullRomCurve &other) -> CatmullRomCurve & = default;
  ~CatmullRomCurve() override                                  = default;

  // Geometry interface
  [[nodiscard]] auto
  geometryTypeId() const -> GeometryType override;
  [[nodiscard]] auto
  geometryType() const -> std::string override;
  [[nodiscard]] auto
  clone() const -> CatmullRomCurve * override;
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

  // CatmullRomCurve specific
  [[nodiscard]] auto
  numPoints() const -> size_t
  {
    return _points.size();
  }
  [[nodiscard]] auto
  point(size_t index) const -> const Point &;
  auto
  point(size_t index) -> Point &;
  void
  setPoint(size_t index, const Point &point);
  void
  addPoint(const Point &point);
  void
  insertPoint(size_t index, const Point &point);
  void
  removePoint(size_t index);
  void
  clear();

  // Tension control (0.5 = standard Catmull-Rom, 0 = linear, 1 = tight)
  [[nodiscard]] auto
  tension() const -> double
  {
    return _tension;
  }
  void
  setTension(double tension);

  // Boundary conditions
  [[nodiscard]] auto
  boundaryType() const -> BoundaryType
  {
    return _boundaryType;
  }
  void
  setBoundaryType(BoundaryType type);

  // Conversion utilities
  [[nodiscard]] auto
  toHermiteCurve() const -> std::unique_ptr<HermiteCurve>;
  [[nodiscard]] auto
  toBSplineCurve() const -> std::unique_ptr<BSplineCurve>;

  // Parameterization options
  enum ParametrizationType : std::uint8_t {
    UNIFORM,      // Uniform parameter spacing
    CHORD_LENGTH, // Based on distance between points
    CENTRIPETAL   // Square root of chord length (smoother)
  };

  [[nodiscard]] auto
  parametrization() const -> ParametrizationType
  {
    return _parametrization;
  }
  void
  setParametrization(ParametrizationType type);

private:
  std::vector<Point>              _points;
  double                          _tension;
  BoundaryType                    _boundaryType;
  ParametrizationType             _parametrization;
  mutable std::vector<Kernel::FT> _parameterValues; // Cached parameter values

  /**
   * @brief Updates the cached parameter values based on current parametrization
   * type
   *
   * This method recalculates the parameter values used for curve evaluation
   * based on the selected parametrization type (uniform, chord length, or
   * centripetal). The parameter values are normalized to the range [0, 1].
   */
  void
  updateParameterValues() const;

  /**
   * @brief Finds the curve segment containing the given parameter value
   * @param parameter The parameter value to locate within the curve segments
   * @return The index of the segment containing the parameter
   *
   * Uses binary search to efficiently locate the appropriate segment for
   * parameter evaluation. Returns the left segment index for the interval
   * containing the parameter.
   */
  [[nodiscard]] auto
  findSegment(double parameter) const -> size_t;

  /**
   * @brief Performs Catmull-Rom interpolation between four control points
   * @param point0 First control point (phantom point before segment)
   * @param point1 Second control point (segment start)
   * @param point2 Third control point (segment end)
   * @param point3 Fourth control point (phantom point after segment)
   * @param parameter Local parameter value within the segment [0,1]
   * @return Interpolated point at the given parameter
   *
   * Applies the Catmull-Rom basis functions with tension control to compute
   * the interpolated position along the curve segment.
   */
  [[nodiscard]] auto
  catmullRomInterpolate(const Point &point0, const Point &point1,
                        const Point &point2, const Point &point3,
                        double parameter) const -> Point;

  /**
   * @brief Calculates the tangent vector at a specific control point
   * @param index Index of the control point
   * @return Tangent vector at the specified point
   *
   * Computes tangent vectors using neighboring points and tension parameter.
   * Handles special cases for first and last points appropriately.
   */
  [[nodiscard]] auto
  getTangent(size_t index) const -> Point;

  /**
   * @brief Retrieves boundary points for curve segments at curve endpoints
   * @param point0 Reference to store the phantom point before segment
   * @param point3 Reference to store the phantom point after segment
   * @param segment Index of the curve segment
   *
   * Generates phantom control points needed for Catmull-Rom interpolation
   * at curve boundaries based on the selected boundary type.
   */
  void
  getBoundaryPoints(Point &point0, Point &point3, size_t segment) const;

  /**
   * @brief Generates the phantom control point for the first curve segment
   * @return Phantom point to use as point0 for first segment interpolation
   *
   * Creates an appropriate phantom point based on the boundary type setting
   * to ensure smooth curve behavior at the start of the curve.
   */
  [[nodiscard]] auto
  getFirstSegmentBoundaryPoint() const -> Point;

  /**
   * @brief Generates the phantom control point for the last curve segment
   * @return Phantom point to use as point3 for last segment interpolation
   *
   * Creates an appropriate phantom point based on the boundary type setting
   * to ensure smooth curve behavior at the end of the curve.
   */
  [[nodiscard]] auto
  getLastSegmentBoundaryPoint() const -> Point;

  /**
   * @brief Calculates derivative for linear interpolation between two points
   * @return Constant derivative vector for linear segment
   *
   * Helper function that computes the derivative for the simple case where
   * the curve contains only two control points.
   */
  [[nodiscard]] auto
  calculateLinearDerivative() const -> Point;

  /**
   * @brief Retrieves the four control points needed for derivative calculation
   * @param segment Index of the curve segment
   * @param point0 Reference to store first control point
   * @param point1 Reference to store second control point
   * @param point2 Reference to store third control point
   * @param point3 Reference to store fourth control point
   *
   * Gathers the appropriate control points including phantom boundary points
   * needed for computing derivatives along the specified curve segment.
   */
  void
  getDerivativeControlPoints(size_t segment, Point &point0, Point &point1,
                             Point &point2, Point &point3) const;

  /**
   * @brief Computes the derivative basis functions for Catmull-Rom
   * interpolation
   * @param localParam Local parameter value within segment [0,1]
   * @return Array containing the four derivative basis function values
   *
   * Calculates the derivatives of the Catmull-Rom basis functions with respect
   * to the parameter, incorporating the tension setting for curve control.
   */
  [[nodiscard]] auto
  calculateDerivativeBasisFunctions(double localParam) const
      -> std::array<double, 4>;
};

} // namespace SFCGAL

#endif // SFCGAL_CATMULLROMCURVE_H
