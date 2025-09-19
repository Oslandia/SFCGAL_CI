// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_NURBSCURVE_H
#define SFCGAL_NURBSCURVE_H

#include "SFCGAL/Curve.h"
#include "SFCGAL/Kernel.h"
#include <cstdint>
#include <memory>
#include <optional>
#include <vector>

namespace SFCGAL {

/**
 * @brief Non-Uniform Rational B-Spline (NURbS) curve implementation
 *
 * Complete OGC/SQL-MM compliant implementation supporting XY, Z, M, ZM
 * dimensions. Weights affect XYZ coordinates in homogeneous space (x*w, y*w,
 * z*w, w). M coordinate is interpolated non-rationally through basis functions.
 *
 * Based on "The NURBS Book" by Piegl & Tiller algorithms and NURBS-Python/geomdl.
 * Uses Kernel::FT (Exact Predicates Exact Constructions) for numerical
 * robustness.
 *
 * ## Key Features
 * - **Complete EndCondition support**: CLAMPED, NATURAL, PERIODIC, TANGENT
 * - **All KnotMethod implementations**: UNIFORM, CHORD_LENGTH, CENTRIPETAL
 * - **Unified FitMethod interface**: INTERPOLATE (exact) vs APPROXIMATE
 * (smooth)
 * - **AutoCAD SPLINE compatibility**: Supports CAD-standard curve fitting
 * - **Mathematical operations**: evaluation, derivatives, arc length,
 * intersections
 * - **Advanced manipulations**: splitting, joining, offsetting, degree
 * elevation
 *
 * ## Usage Examples
 * @code
 * // Simple interpolation (exact fit through points)
 * auto curve1 = NURBSCurve::interpolateCurve(points, 3);
 *
 * // Natural end conditions for smoother curves
 * auto curve2 = NURBSCurve::interpolateCurve(points, 3,
 *     NURBSCurve::KnotMethod::CENTRIPETAL,
 *     NURBSCurve::EndCondition::NATURAL);
 *
 * // Unified interface with explicit control
 * auto curve3 = NURBSCurve::fitCurve(points, 3,
 *     NURBSCurve::FitMethod::INTERPOLATE,
 *     NURBSCurve::KnotMethod::CHORD_LENGTH,
 *     NURBSCurve::EndCondition::CLAMPED);
 * @endcode
 *
 * @ingroup public_api
 * @since 2.3.0
 */
class SFCGAL_API NURBSCurve : public Curve {
public:
  // Type aliases for clarity and consistency
  using FT           = Kernel::FT;
  using ControlPoint = Point;
  using Weight       = FT;
  using Knot         = FT;
  using Parameter    = FT;

  /**
   * @brief Knot vector generation methods for parameterization
   */
  enum class KnotMethod : std::uint8_t {
    UNIFORM,      ///< Equal knot spacing (fastest generation)
    CHORD_LENGTH, ///< Based on Euclidean distances between points
    CENTRIPETAL ///< Square root of chord length (often best shape preservation)
  };

  /**
   * @brief End conditions for curve fitting and interpolation
   *
   * Controls the behavior at curve endpoints to achieve different
   * mathematical and aesthetic properties.
   */
  enum class EndCondition : std::uint8_t {
    CLAMPED,  ///< Multiplicity degree+1 at curve ends (standard CAD behavior)
    NATURAL,  ///< Minimal curvature conditions (smoother, more organic curves)
    PERIODIC, ///< Closed periodic curve with C^(degree-1) continuity at
              ///< junction
    TANGENT ///< User-specified tangent constraints at ends (future enhancement)
  };

  /**
   * @brief Curve fitting methods
   *
   * Determines the mathematical approach to curve generation from input points.
   */
  enum class FitMethod : std::uint8_t {
    INTERPOLATE, ///< Pass exactly through all points (higher fidelity, may be
                 ///< less smooth)
    APPROXIMATE  ///< Least squares approximation within tolerance (smoother,
                 ///< fewer control points)
  };


  // Constructors

  /**
   * @brief Default constructor creating empty curve
   */
  NURBSCurve();

  /**
   * @brief Constructor with control points and uniform weights
   * @param controlPoints Control points defining curve shape
   * @param degree Polynomial degree of curve (must be < controlPoints.size())
   * @param knotMethod Method for generating knot vector
   */
  NURBSCurve(const std::vector<Point> &controlPoints, unsigned int degree,
             KnotMethod knotMethod = KnotMethod::UNIFORM);

  /**
   * @brief Constructor with control points and explicit weights
   * @param controlPoints Control points defining curve shape
   * @param weights Weight for each control point (must match size)
   * @param degree Polynomial degree of curve
   * @param knotMethod Method for generating knot vector
   */
  NURBSCurve(const std::vector<Point> &controlPoints,
             const std::vector<FT> &weights, unsigned int degree,
             KnotMethod knotMethod = KnotMethod::UNIFORM);

  /**
   * @brief Full constructor with all NURBS parameters
   * @param controlPoints Control points defining curve shape
   * @param weights Weights (empty for uniform weights = 1.0)
   * @param degree Polynomial degree
   * @param knotVector Custom knot vector (must satisfy m = n + p + 1)
   */
  NURBSCurve(const std::vector<Point> &controlPoints,
             const std::vector<FT> &weights, unsigned int degree,
             const std::vector<Knot> &knotVector);

  /**
   * @brief Copy constructor
   */
  NURBSCurve(const NURBSCurve &other) = default;

  /**
   * @brief Assignment operator
   */
  auto
  operator=(const NURBSCurve &other) -> NURBSCurve & = default;

  /**
   * @brief Destructor
   */
  ~NURBSCurve() override = default;

  // Factory methods for common curve types

  /**
   * @brief Create from Bezier control points
   * @param controlPoints Bezier control points
   * @return NURBS curve with clamped knot vector [0,...,0,1,...,1]
   */
  static auto
  fromBezier(const std::vector<Point> &controlPoints)
      -> std::unique_ptr<NURBSCurve>;

  /**
   * @brief Create from B-spline parameters
   * @param controlPoints B-spline control points
   * @param degree Spline degree
   * @param knots Knot vector
   * @return NURBS curve with uniform weights
   */
  static auto
  fromBSpline(const std::vector<Point> &controlPoints, unsigned int degree,
              const std::vector<Knot> &knots) -> std::unique_ptr<NURBSCurve>;

  /**
   * @brief Create uniform B-spline with standard parameterization
   * @param controlPoints Control points
   * @param degree Spline degree
   * @return NURBS curve with uniform knots and weights
   */
  static auto
  createBSpline(const std::vector<Point> &controlPoints, unsigned int degree)
      -> std::unique_ptr<NURBSCurve>;

  /**
   * @brief Create exact circular arc representation
   * @param center Arc center point
   * @param radius Arc radius (must be positive)
   * @param startAngle Start angle in radians
   * @param endAngle End angle in radians
   * @param normal Normal vector for 3D arcs (default Z-axis)
   * @return Exact NURBS representation of circular arc
   */
  static auto
  createCircularArc(const Point &center, FT radius, FT startAngle, FT endAngle,
                    const Point &normal = Point(0, 0, 1), KnotMethod knotMethod = KnotMethod::UNIFORM)
      -> std::unique_ptr<NURBSCurve>;

  /**
   * @brief Interpolate points exactly with specified continuity
   * @param points Points to interpolate through
   * @param degree Curve degree (will be reduced if too high)
   * @param knotMethod Parameterization method for stability
   * @param endCondition Boundary conditions
   * @return Interpolating NURBS curve
   */
  static auto
  interpolateCurve(const std::vector<Point> &points, unsigned int degree = 3,
                   KnotMethod   knotMethod   = KnotMethod::UNIFORM,
                   EndCondition endCondition = EndCondition::CLAMPED)
      -> std::unique_ptr<NURBSCurve>;

  /**
   * @brief Approximate points within tolerance using least squares
   * @param points Points to approximate
   * @param degree Target curve degree
   * @param tolerance Maximum allowed deviation
   * @param maxControlPoints Maximum control points to use
   * @return Approximating NURBS curve within tolerance
   */
  static auto
  approximateCurve(const std::vector<Point> &points, unsigned int degree,
                   FT tolerance, size_t maxControlPoints = 50)
      -> std::unique_ptr<NURBSCurve>;

  /**
   * @brief Unified curve fitting interface using FitMethod enum
   * @param points Points to fit curve through or approximate
   * @param degree Target curve degree
   * @param fitMethod Whether to interpolate exactly or approximate within
   * tolerance
   * @param knotMethod Parameterization method for stability
   * @param endCondition Boundary conditions for interpolation
   * @param tolerance Maximum deviation for approximation (ignored for
   * interpolation)
   * @param maxControlPoints Maximum control points for approximation (ignored
   * for interpolation)
   * @return NURBS curve fitting the specified points
   */
  static auto
  fitCurve(const std::vector<Point> &points, unsigned int degree = 3,
           FitMethod         fitMethod         = FitMethod::INTERPOLATE,
           KnotMethod        knotMethod        = KnotMethod::UNIFORM,
           EndCondition      endCondition      = EndCondition::CLAMPED,
           FT                tolerance         = FT(1e-6),
           size_t            maxControlPoints  = 50)
      -> std::unique_ptr<NURBSCurve>;

  // Geometry interface implementation

  void
  accept(GeometryVisitor &visitor) override;

  void
  accept(ConstGeometryVisitor &visitor) const override;

  [[nodiscard]] auto
  geometryTypeId() const -> GeometryType override;

  [[nodiscard]] auto
  geometryType() const -> std::string override;

  [[nodiscard]] auto
  clone() const -> NURBSCurve * override;

  [[nodiscard]] auto
  isEmpty() const -> bool override;

  [[nodiscard]] auto
  is3D() const -> bool override;

  [[nodiscard]] auto
  isMeasured() const -> bool override;

  auto
  dropZ() -> bool override;

  auto
  dropM() -> bool override;

  void
  swapXY() override;

  // Curve interface implementation

  [[nodiscard]] auto
  evaluate(Parameter parameter) const -> Point override;

  [[nodiscard]] auto
  derivative(Parameter parameter, unsigned int order = 1) const
      -> Point override;

  [[nodiscard]] auto
  tangent(Parameter parameter) const -> Point override;

  [[nodiscard]] auto
  normal(Parameter parameter) const -> Point override;

  [[nodiscard]] auto
  binormal(Parameter parameter) const -> Point override;

  [[nodiscard]] auto
  curvature(Parameter parameter) const -> FT override;

  [[nodiscard]] auto
  torsion(Parameter parameter) const -> FT override;

  [[nodiscard]] auto
  frenetFrame(Parameter parameter) const
      -> std::tuple<Point, Point, Point> override;

  [[nodiscard]] auto
  toLineString(unsigned int numSegments = 32) const
      -> std::unique_ptr<LineString> override;

  [[nodiscard]] auto
  toLineStringAdaptive(FT tolerance = FT(1e-3), unsigned int minSegments = 8,
                       unsigned int maxSegments = 256) const
      -> std::unique_ptr<LineString> override;

  /**
   * @brief Convert curve to LineString using arc-length-based sampling
   * @param numSegments Number of line segments to use
   * @return LineString with points sampled uniformly in arc-length space
   */
  [[nodiscard]] auto
  toLineStringArcLength(unsigned int numSegments = 32) const
      -> std::unique_ptr<LineString>;

  [[nodiscard]] auto
  parameterBounds() const -> std::pair<Parameter, Parameter> override;

  [[nodiscard]] auto
  isPeriodic() const -> bool override;

  [[nodiscard]] auto
  isPlanar(std::vector<FT> *plane = nullptr) const -> bool override;

  [[nodiscard]] auto
  isLinear() const -> bool override;

  [[nodiscard]] auto
  degree() const -> unsigned int override
  {
    return _degree;
  }

  [[nodiscard]] auto
  length(Parameter from = Parameter(-1), Parameter to = Parameter(-1),
         FT tolerance = FT(1e-6)) const -> FT override;

  [[nodiscard]] auto
  parameterAtLength(FT arcLength, FT tolerance = FT(1e-6)) const
      -> Parameter override;

  [[nodiscard]] auto
  reparameterizeByArcLength() const -> std::unique_ptr<Curve> override;

  [[nodiscard]] auto
  split(Parameter parameter) const
      -> std::pair<std::unique_ptr<Curve>, std::unique_ptr<Curve>> override;

  [[nodiscard]] auto
  subcurve(Parameter from, Parameter to) const
      -> std::unique_ptr<Curve> override;

  [[nodiscard]] auto
  reverse() const -> std::unique_ptr<Curve> override;

  [[nodiscard]] auto
  join(const Curve &other, Continuity continuity = Continuity::C0,
       FT tolerance = FT(1e-6)) const -> std::unique_ptr<Curve> override;

  [[nodiscard]] auto
  closestPoint(const Point &point, Parameter *outParameter = nullptr) const
      -> Point override;

  // NURBS-specific access and manipulation methods

  /**
   * @brief Get control points array
   * @return Const reference to control points vector
   */
  [[nodiscard]] auto
  controlPoints() const -> const std::vector<Point> &
  {
    return _controlPoints;
  }

  /**
   * @brief Set all control points with validation
   * @param controlPoints New control points
   * @throws Exception if dimensionally inconsistent or degree constraint
   * violated
   */
  void
  setControlPoints(const std::vector<Point> &controlPoints);

  /**
   * @brief Access control point for modification
   * @param index Control point index (must be valid)
   * @return Non-const reference to control point
   * @throws Exception if index out of bounds
   */
  [[nodiscard]] auto
  controlPointN(size_t index) -> Point & override;

  /**
   * @brief Access control point for reading
   * @param index Control point index (must be valid)
   * @return Const reference to control point
   * @throws Exception if index out of bounds
   */
  [[nodiscard]] auto
  controlPointN(size_t index) const -> const Point & override;

  /**
   * @brief Set individual control point with validation
   * @param index Control point index
   * @param point New control point value
   * @throws Exception if index invalid or dimension mismatch
   */
  void
  setControlPoint(size_t index, const Point &point);

  /**
   * @brief Get weights array
   * @return Const reference to weights vector
   */
  [[nodiscard]] auto
  weights() const -> const std::vector<FT> &
  {
    return _weights;
  }

  /**
   * @brief Set all weights with validation
   * @param weights New weights (must be positive and match control point count)
   * @throws Exception if invalid weights provided
   */
  void
  setWeights(const std::vector<FT> &weights);

  /**
   * @brief Get weight at specified index
   * @param index Weight index
   * @return Weight value (1.0 if using uniform weights)
   * @throws Exception if index out of bounds
   */
  [[nodiscard]] auto
  weight(size_t index) const -> FT override;

  /**
   * @brief Set weight at specified index
   * @param index Weight index
   * @param weight New weight value (must be positive)
   * @throws Exception if index invalid or weight non-positive
   */
  void
  setWeight(size_t index, FT weight);

  /**
   * @brief Check if curve uses non-uniform rational weights
   * @return true if weights vary significantly from uniform
   */
  [[nodiscard]] auto
  isRational() const -> bool override;

  /**
   * @brief Check if curve is in Bezier form
   * @return true if knot vector has Bezier structure [0,...,0,1,...,1]
   */
  [[nodiscard]] auto
  isBezier() const -> bool;

  /**
   * @brief Check if curve is B-spline (uniform weights)
   * @return true if all weights are effectively equal
   */
  [[nodiscard]] auto
  isBSpline() const -> bool;

  /**
   * @brief Get knot vector
   * @return Const reference to knot vector
   */
  [[nodiscard]] auto
  knotVector() const -> const std::vector<Knot> &
  {
    return _knotVector;
  }

  /**
   * @brief Set knot vector with full validation
   * @param knots New knot vector (must satisfy NURBS constraints)
   * @throws Exception if knot vector invalid for current curve parameters
   */
  void
  setKnotVector(const std::vector<Knot> &knots);

  /**
   * @brief Get knot value at index
   * @param index Knot index
   * @return Knot value
   * @throws Exception if index out of bounds
   */
  [[nodiscard]] auto
  knot(size_t index) const -> Knot;

  /**
   * @brief Count multiplicity of knot value
   * @param value Knot value to count
   * @param tolerance Tolerance for knot equality
   * @return Number of occurrences of knot value
   */
  [[nodiscard]] auto
  knotMultiplicity(Knot value, FT tolerance = FT(1e-10)) const -> unsigned int;

  /**
   * @brief Get curve order (degree + 1)
   * @return Curve order
   */
  [[nodiscard]] auto
  order() const -> unsigned int
  {
    return _degree + 1;
  }

  /**
   * @brief Get number of control points
   * @return Control point count
   */
  [[nodiscard]] auto
  numControlPoints() const -> size_t override
  {
    return _controlPoints.size();
  }

  /**
   * @brief Get number of knots in vector
   * @return Knot count
   */
  [[nodiscard]] auto
  numKnots() const -> size_t
  {
    return _knotVector.size();
  }

  // Advanced NURBS operations

  /**
   * @brief Insert knot into curve without changing shape (Oslo algorithm)
   * @param parameter Knot value to insert
   * @param times Number of insertions (limited by degree)
   * @return New curve with inserted knot
   * @throws Exception if parameter outside knot span or too many insertions
   */
  [[nodiscard]] auto
  insertKnot(Knot parameter, unsigned int times = 1) const
      -> std::unique_ptr<NURBSCurve>;

  /**
   * @brief Refine knot vector by adding multiple knots
   * @param newKnots Sorted knots to insert
   * @return Refined curve (shape unchanged)
   * @throws Exception if knots not properly ordered
   */
  [[nodiscard]] auto
  refineKnotVector(const std::vector<Knot> &newKnots) const
      -> std::unique_ptr<NURBSCurve>;

  /**
   * @brief Elevate degree of curve (increases smoothness)
   * @param times Number of degree elevations
   * @return Curve with elevated degree (shape unchanged)
   * @throws Exception if elevation would create too many control points
   */
  [[nodiscard]] auto
  elevateDegree(unsigned int times = 1) const -> std::unique_ptr<NURBSCurve>;

  /**
   * @brief Reduce degree if possible without shape change
   * @param tolerance Maximum allowed deviation
   * @return Degree-reduced curve or nullptr if reduction impossible within
   * tolerance
   */
  [[nodiscard]] auto
  reduceDegree(FT tolerance = FT(1e-6)) const -> std::unique_ptr<NURBSCurve>;

  /**
   * @brief Remove unnecessary knots while preserving shape
   * @param tolerance Maximum allowed deviation
   * @return Simplified curve with fewer knots
   */
  [[nodiscard]] auto
  removeKnots(FT tolerance = FT(1e-6)) const -> std::unique_ptr<NURBSCurve>;

  // Fitting and approximation data access

  /**
   * @brief Get original fit points if curve was created by fitting
   * @return Original fit points (empty if not applicable)
   */
  [[nodiscard]] auto
  fitPoints() const -> const std::vector<Point> &
  {
    return _fitPoints;
  }

  /**
   * @brief Get fitting tolerance used
   * @return Fit tolerance (0.0 for exact interpolation)
   */
  [[nodiscard]] auto
  fitTolerance() const -> FT
  {
    return _fitTolerance;
  }

  /**
   * @brief Get start tangent constraint if specified
   * @return Optional start tangent vector
   */
  [[nodiscard]] auto
  startTangent() const -> const std::optional<Point> &
  {
    return _startTangent;
  }

  /**
   * @brief Get end tangent constraint if specified
   * @return Optional end tangent vector
   */
  [[nodiscard]] auto
  endTangent() const -> const std::optional<Point> &
  {
    return _endTangent;
  }

  /**
   * @brief Comprehensive NURBS data validation
   * @return Pair of (is_valid, error_description)
   */
  [[nodiscard]] auto
  validateData() const -> std::pair<bool, std::string>;

  /**
   * @brief Get detailed curve statistics for debugging
   * @return Map of curve properties and values
   */
  [[nodiscard]] auto
  getCurveStatistics() const -> std::map<std::string, double>;

protected:
  // Core NURBS data with proper encapsulation
  std::vector<Point> _controlPoints; ///< Control points in curve space
  std::vector<FT>    _weights;       ///< Weights (empty = all 1.0)
  unsigned int       _degree{0};     ///< Polynomial degree
  std::vector<Knot>  _knotVector;    ///< Non-decreasing knot sequence

  // Optional fitting metadata for CAD compatibility
  std::vector<Point>   _fitPoints;           ///< Original points if fitted
  FT                   _fitTolerance{FT(0)}; ///< Fitting tolerance used
  std::optional<Point> _startTangent;        ///< Start tangent constraint
  std::optional<Point> _endTangent;          ///< End tangent constraint

  // Core algorithmic methods

  /**
   * @brief Find knot span containing parameter (Algorithm A2.1)
   * @param parameter Parameter value to locate
   * @return Knot span index
   */
  [[nodiscard]] auto
  findSpan(Parameter parameter) const -> size_t;

  /**
   * @brief Compute basis functions and derivatives (Algorithm A2.3)
   * @param span Knot span index
   * @param parameter Parameter value
   * @param maxDerivative Maximum derivative order needed
   * @return Matrix [derivative_order][basis_function] of values
   */
  [[nodiscard]] auto
  basisFunctionDerivatives(size_t span, Parameter parameter,
                           unsigned int maxDerivative) const
      -> std::vector<std::vector<FT>>;

  /**
   * @brief Evaluate curve using De Boor's algorithm for rational curves
   * @param span Knot span containing parameter
   * @param parameter Parameter value
   * @return Point on curve
   */
  [[nodiscard]] auto
  deBoorRational(size_t span, Parameter parameter) const -> Point;

  /**
   * @brief Compute curve derivatives up to specified order (Algorithm A4.2)
   * @param parameter Evaluation parameter
   * @param maxOrder Maximum derivative order
   * @return Vector of derivatives [C(u), C'(u), C''(u), ...]
   */
  [[nodiscard]] auto
  computeDerivatives(Parameter parameter, unsigned int maxOrder) const
      -> std::vector<Point>;

  // Utility methods for curve generation and validation

  /**
   * @brief Generate knot vector using specified method
   * @param points Reference points for parameterization
   * @param degree Curve degree
   * @param method Knot generation method
   * @return Generated knot vector
   */
  static auto
  generateKnotVector(const std::vector<Point> &points, unsigned int degree,
                     KnotMethod method) -> std::vector<Knot>;

  /**
   * @brief Compute parameter values for point sequence
   * @param points Points to parameterize
   * @param method Parameterization method
   * @return Normalized parameter values [0,1]
   */
  static auto
  computeParameters(const std::vector<Point> &points, KnotMethod method)
      -> std::vector<Parameter>;

  /**
   * @brief Check dimensional consistency across all control points
   * @return true if all points have same XYZ and M dimensions
   */
  [[nodiscard]] auto
  checkDimensionalConsistency() const -> bool;

  /**
   * @brief Validate complete NURBS data structure
   * @param controlPoints Control points to validate
   * @param weights Weights to validate
   * @param degree Degree to validate
   * @param knots Knot vector to validate
   * @return Pair of (is_valid, detailed_error_message)
   */
  static auto
  validateNURBSData(const std::vector<Point> &controlPoints,
                    const std::vector<FT> &weights, unsigned int degree,
                    const std::vector<Knot> &knots)
      -> std::pair<bool, std::string>;

  // Advanced mathematical operations

  /**
   * @brief Compute arc length using adaptive Simpson quadrature
   * @param startParam Start parameter
   * @param endParam End parameter
   * @param tolerance Integration tolerance
   * @return Arc length between parameters
   */
  [[nodiscard]] auto
  computeArcLength(Parameter startParam, Parameter endParam, FT tolerance) const
      -> FT;

  /**
   * @brief Compute speed function ||C'(t)|| for arc length integration
   * @param t Parameter value
   * @return Magnitude of derivative vector at parameter t
   */
  [[nodiscard]] auto
  speedFunction(Parameter t) const -> FT;

  /**
   * @brief Apply Simpson's rule over interval [a,b]
   * @param a Start parameter
   * @param b End parameter
   * @return Simpson approximation of integral
   */
  [[nodiscard]] auto
  simpsonRule(Parameter a, Parameter b) const -> FT;

  /**
   * @brief Adaptive Simpson quadrature for arc length computation
   * @param a Start parameter
   * @param b End parameter
   * @param tolerance Integration tolerance
   * @param maxDepth Maximum recursion depth (default 20)
   * @return Accurate arc length using adaptive subdivision
   */
  [[nodiscard]] auto
  adaptiveSimpsonArcLength(Parameter a, Parameter b, FT tolerance,
                           unsigned int maxDepth = 20) const -> FT;

  /**
   * @brief Find parameter corresponding to arc length using Newton-Raphson
   * @param targetLength Target arc length from start
   * @param tolerance Convergence tolerance
   * @return Parameter value at specified arc length
   */
  [[nodiscard]] auto
  findParameterByArcLength(FT targetLength, FT tolerance) const -> Parameter;

  /**
   * @brief Project point onto curve using Newton's method
   * @param point Point to project
   * @param tolerance Convergence tolerance
   * @param maxIterations Maximum iterations
   * @return Pair of (closest_point_on_curve, parameter_value)
   */
  [[nodiscard]] auto
  projectPointToCurve(const Point &point, FT tolerance = FT(1e-9),
                      unsigned int maxIterations = 50) const
      -> std::pair<Point, Parameter>;

  // Memory and performance optimization helpers

  /**
   * @brief Reserve appropriate capacity for internal vectors
   * @param expectedSize Expected final size for optimization
   */
  void
  reserveCapacity(size_t expectedSize);

  /**
   * @brief Clear all optional data to minimize memory usage
   */
  void
  clearOptionalData();

  // Helper methods for end condition implementations

  /**
   * @brief Generate knot vector appropriate for specified end condition
   */
  static auto
  generateKnotVectorForEndCondition(const std::vector<Parameter> &parameters,
                                    unsigned int                  degree,
                                    EndCondition                  endCondition)
      -> std::vector<Knot>;

  /**
   * @brief Generate knot vector for least-squares approximation
   * @param parameters Parameter values for data points
   * @param degree Curve degree
   * @param numControlPoints Number of control points (m in Piegl & Tiller)
   * @return Knot vector optimized for approximation
   */
  static auto
  generateApproximationKnotVector(const std::vector<Parameter> &parameters,
                                  unsigned int degree, size_t numControlPoints)
      -> std::vector<Knot>;

  /**
   * @brief Interpolate with clamped end conditions (original behavior)
   */
  static auto
  interpolateClampedCurve(const std::vector<Point>     &points,
                          const std::vector<Parameter> &parameters,
                          unsigned int degree, const std::vector<Knot> &knots)
      -> std::vector<Point>;

  /**
   * @brief Interpolate with natural end conditions (minimal curvature)
   */
  static auto
  interpolateNaturalCurve(const std::vector<Point>     &points,
                          const std::vector<Parameter> &parameters,
                          unsigned int degree, const std::vector<Knot> &knots)
      -> std::vector<Point>;

  /**
   * @brief Create periodic interpolating curve
   */
  static auto
  interpolatePeriodicCurve(const std::vector<Point> &points,
                           unsigned int degree, KnotMethod knotMethod)
      -> std::unique_ptr<NURBSCurve>;
};

} // namespace SFCGAL

#endif // SFCGAL_NURBSCURVE_H
