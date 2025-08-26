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
 * @brief Non-Uniform Rational B-Spline (NURBS) curve implementation
 *
 * OGC/SQL-MM compliance: supports XY, Z, M, ZM dimensions.
 * Weights affect XYZ in homogeneous coordinates (x*w, y*w, z*w, w).
 * M coordinate is interpolated non-rationally through basis functions.
 *
 * Based on "The NURBS Book" by Piegl & Tiller formulations.
 * Uses Kernel::FT (Epeck) for numerical robustness.
 *
 * @ingroup public_api
 * @since 2.3.0
 */
class SFCGAL_API NURBSCurve : public Curve {
public:
  // Type aliases for clarity
  using FT           = Kernel::FT;
  using ControlPoint = Point;
  using Weight       = FT;
  using Knot         = FT;
  using Parameter    = FT;

  /**
   * @brief Knot vector generation methods
   */
  enum class KnotMethod : std::uint8_t {
    UNIFORM,      ///< Uniform knot spacing
    CHORD_LENGTH, ///< Based on Euclidean distances
    CENTRIPETAL   ///< Square root of chord length
  };

  /**
   * @brief End conditions for curve fitting
   */
  enum class EndCondition : std::uint8_t {
    CLAMPED,  ///< Multiplicity degree+1 at ends
    NATURAL,  ///< Minimal curvature (cubic spline)
    PERIODIC, ///< Closed periodic curve
    TANGENT   ///< Imposed tangent constraints
  };

  /**
   * @brief Fitting methods
   */
  enum class FitMethod : std::uint8_t {
    INTERPOLATE, ///< Pass exactly through points
    APPROXIMATE  ///< Least squares fitting
  };

  // Constructors

  /**
   * @brief Default constructor creating empty curve
   */
  NURBSCurve();

  /**
   * @brief Constructor with control points and uniform weights
   * @param controlPoints Control points
   * @param degree Curve degree
   * @param knotMethod Method for generating knot vector
   */
  NURBSCurve(const std::vector<Point> &controlPoints, unsigned int degree,
             KnotMethod knotMethod = KnotMethod::CHORD_LENGTH);

  /**
   * @brief Constructor with control points and explicit weights
   * @param controlPoints Control points
   * @param weights Weight for each control point
   * @param degree Curve degree
   * @param knotMethod Method for generating knot vector
   */
  NURBSCurve(const std::vector<Point> &controlPoints,
             const std::vector<FT> &weights, unsigned int degree,
             KnotMethod knotMethod = KnotMethod::CHORD_LENGTH);

  /**
   * @brief Full constructor with all parameters
   * @param controlPoints Control points
   * @param weights Weights (empty for uniform)
   * @param degree Curve degree
   * @param knotVector Custom knot vector
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

  // Factory methods

  /**
   * @brief Create from Bezier control points
   * @param controlPoints Bezier control points
   * @return NURBS curve with Bezier knot vector
   */
  static auto
  fromBezier(const std::vector<Point> &controlPoints)
      -> std::unique_ptr<NURBSCurve>;

  /**
   * @brief Create from B-spline
   * @param controlPoints B-spline control points
   * @param degree Spline degree
   * @param knots Knot vector
   * @return NURBS curve with uniform weights
   */
  static auto
  fromBSpline(const std::vector<Point> &controlPoints, unsigned int degree,
              const std::vector<Knot> &knots) -> std::unique_ptr<NURBSCurve>;

  /**
   * @brief Create uniform B-spline
   * @param controlPoints Control points
   * @param degree Spline degree
   * @return NURBS curve with uniform knots and weights
   */
  static auto
  createBSpline(const std::vector<Point> &controlPoints, unsigned int degree)
      -> std::unique_ptr<NURBSCurve>;

  /**
   * @brief Create circular arc
   * @param center Arc center
   * @param radius Arc radius
   * @param startAngle Start angle in radians
   * @param endAngle End angle in radians
   * @param normal Normal vector (for 3D arcs)
   * @return Exact NURBS representation of arc
   */
  static auto
  createCircularArc(const Point &center, FT radius, FT startAngle, FT endAngle,
                    const Point &normal = Point(0, 0, 1))
      -> std::unique_ptr<NURBSCurve>;

  /**
   * @brief Interpolate points exactly (Piegl & Tiller Ch.9)
   * @param points Points to interpolate
   * @param degree Curve degree
   * @param knotMethod Parameterization method
   * @param endCondition Boundary conditions
   * @return Interpolating NURBS curve
   */
  static auto
  interpolateCurve(const std::vector<Point> &points, unsigned int degree = 3,
                   KnotMethod   knotMethod   = KnotMethod::CENTRIPETAL,
                   EndCondition endCondition = EndCondition::CLAMPED)
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

  [[nodiscard]] auto
  parameterBounds() const -> std::pair<Parameter, Parameter> override;

  [[nodiscard]] auto
  isClosed() const -> bool override;

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
  offset(FT distance) const -> std::unique_ptr<Curve> override;

  [[nodiscard]] auto
  closestPoint(const Point &point, Parameter *outParameter = nullptr) const
      -> Point override;

  [[nodiscard]] auto
  distance(const Point &point) const -> FT override;

  [[nodiscard]] auto
  hasSelfIntersections(std::vector<std::pair<Parameter, Parameter>>
                           *intersections = nullptr) const -> bool override;

  [[nodiscard]] auto
  intersect(const Curve &other, FT tolerance = FT(1e-6)) const
      -> std::vector<std::tuple<Point, Parameter, Parameter>> override;

  [[nodiscard]] auto
  boundingBox() const -> std::pair<Point, Point> override;

  // NURBS-specific access methods

  /**
   * @brief Get control points
   * @return Vector of control points
   */
  [[nodiscard]] auto
  controlPoints() const -> const std::vector<Point> &
  {
    return _controlPoints;
  }

  /**
   * @brief Set control points
   * @param controlPoints New control points
   */
  void
  setControlPoints(const std::vector<Point> &controlPoints);

  /**
   * @brief Get control point at index (non-const version)
   * @param index Control point index
   * @return Control point reference
   */
  [[nodiscard]] auto
  controlPointN(size_t index) -> Point &;

  /**
   * @brief Get control point at index
   * @param index Control point index
   * @return Control point
   */
  [[nodiscard]] auto
  controlPointN(size_t index) const -> const Point &;

  /**
   * @brief Set control point at index
   * @param index Control point index
   * @param point New control point value
   */
  void
  setControlPoint(size_t index, const Point &point);

  /**
   * @brief Get weights
   * @return Vector of weights
   */
  [[nodiscard]] auto
  weights() const -> const std::vector<FT> &
  {
    return _weights;
  }

  /**
   * @brief Set weights
   * @param weights New weights
   */
  void
  setWeights(const std::vector<FT> &weights);

  /**
   * @brief Get weight at index
   * @param index Weight index
   * @return Weight value
   */
  [[nodiscard]] auto
  weight(size_t index) const -> FT;

  /**
   * @brief Set weight at index
   * @param index Weight index
   * @param weight New weight value
   */
  void
  setWeight(size_t index, FT weight);

  /**
   * @brief Check if curve is rational (non-uniform weights)
   * @return true if weights are non-uniform
   */
  [[nodiscard]] auto
  isRational() const -> bool;

  /**
   * @brief Check if curve is a Bezier curve
   * @return true if knot vector has Bezier form
   */
  [[nodiscard]] auto
  isBezier() const -> bool;

  /**
   * @brief Check if curve is a B-spline (uniform weights)
   * @return true if all weights are equal
   */
  [[nodiscard]] auto
  isBSpline() const -> bool;

  /**
   * @brief Get knot vector
   * @return Knot vector
   */
  [[nodiscard]] auto
  knotVector() const -> const std::vector<Knot> &
  {
    return _knotVector;
  }

  /**
   * @brief Set knot vector
   * @param knots New knot vector
   */
  void
  setKnotVector(const std::vector<Knot> &knots);

  /**
   * @brief Get knot at index
   * @param index Knot index
   * @return Knot value
   */
  [[nodiscard]] auto
  knot(size_t index) const -> Knot;

  /**
   * @brief Get knot multiplicity
   * @param value Knot value to check
   * @return Multiplicity of the knot
   */
  [[nodiscard]] auto
  knotMultiplicity(Knot value) const -> unsigned int;

  /**
   * @brief Get order (degree + 1)
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
  numControlPoints() const -> size_t
  {
    return _controlPoints.size();
  }

  /**
   * @brief Get number of knots
   * @return Knot count
   */
  [[nodiscard]] auto
  numKnots() const -> size_t
  {
    return _knotVector.size();
  }

  // NURBS operations

  /**
   * @brief Insert a knot
   * @param parameter Knot value to insert
   * @param times Number of times to insert
   * @return New curve with inserted knot
   */
  auto
  insertKnot(Knot parameter, unsigned int times = 1)
      -> std::unique_ptr<NURBSCurve>;

  /**
   * @brief Refine knot vector
   * @param newKnots Knots to insert
   * @return Refined curve
   */
  auto
  refineKnotVector(const std::vector<Knot> &newKnots)
      -> std::unique_ptr<NURBSCurve>;

  /**
   * @brief Elevate degree
   * @param times Number of degree elevations
   * @return Curve with elevated degree
   */
  auto
  elevateDegree(unsigned int times = 1) -> std::unique_ptr<NURBSCurve>;

  /**
   * @brief Get fit points if curve was created by fitting
   * @return Original fit points (empty if not applicable)
   */
  [[nodiscard]] auto
  fitPoints() const -> const std::vector<Point> &
  {
    return _fitPoints;
  }

  /**
   * @brief Get fit tolerance if curve was created by fitting
   * @return Fit tolerance (0 if exact interpolation)
   */
  [[nodiscard]] auto
  fitTolerance() const -> FT
  {
    return _fitTolerance;
  }

  /**
   * @brief Validate NURBS data consistency
   * @return <true if data is valid, the invalidity raeson>
   */
  [[nodiscard]] auto
  validateData() const -> std::pair<bool, std::string>;

protected:
  // Core NURBS data
  std::vector<Point> _controlPoints; ///< Control points (XYZM)
  std::vector<FT>    _weights;       ///< Weights (empty = uniform)
  unsigned int       _degree{0};     ///< Polynomial degree
  std::vector<Knot>  _knotVector;    ///< Non-decreasing knot vector

  // Optional fit data (AutoCAD compatibility)
  std::vector<Point>   _fitPoints;           ///< Original fit points
  FT                   _fitTolerance{FT(0)}; ///< Fit tolerance
  std::optional<Point> _startTangent;        ///< Start tangent constraint
  std::optional<Point> _endTangent;          ///< End tangent constraint

  // Helper methods

  /**
   * @brief Find knot span for parameter (Piegl-Tiller Algorithm A2.1)
   * @param parameter Parameter value
   * @return Knot span index
   */
  [[nodiscard]] auto
  findSpan(Parameter parameter) const -> size_t;

  /**
   * @brief Compute basis function derivatives
   * @param span Knot span
   * @param parameter Parameter value
   * @param order Derivative order
   * @return Matrix of derivatives
   */
  [[nodiscard]] auto
  basisFunctionDerivatives(size_t span, Parameter parameter,
                           unsigned int order) const
      -> std::vector<std::vector<FT>>;

  /**
   * @brief Generate knot vector by method
   * @param points Reference points
   * @param degree Curve degree
   * @param method Generation method
   * @return Generated knot vector
   */
  static auto
  generateKnotVector(const std::vector<Point> &points, unsigned int degree,
                     KnotMethod method) -> std::vector<Knot>;

  /**
   * @brief Check dimensional consistency of control points
   * @return true if all points have same dimensions
   */
  [[nodiscard]] auto
  checkDimensionalConsistency() const -> bool;

  /**
   * @brief De Boor's algorithm for rational evaluation
   * @param span Knot span
   * @param parameter Parameter value
   * @return Evaluated point
   */
  [[nodiscard]] auto
  deBoorRational(size_t span, Parameter parameter) const -> Point;

  /**
   * @brief Compute parameter values for point sequence
   * @param points Points to parameterize
   * @param method Parameterization method
   * @return Parameter values
   */
  static auto
  computeParameters(const std::vector<Point> &points, KnotMethod method)
      -> std::vector<Parameter>;

  /**
   * @brief Compute all derivatives up to order (NURBS Book A4.2)
   * @param parameter Evaluation parameter
   * @param maxOrder Maximum derivative order
   * @return Vector [C(u), C'(u), C''(u), ..., C^(maxOrder)(u)]
   */
  auto
  derivativesAt(Parameter parameter, unsigned int maxOrder) const
      -> std::vector<Point>;

  /**
   * @brief Setup collocation matrix for interpolation
   * @param parameters Parameter values
   * @param degree Curve degree
   * @param knots Knot vector
   * @return Collocation matrix
   */
  static auto
  setupCollocationMatrix(const std::vector<Parameter> &parameters,
                         unsigned int degree, const std::vector<Knot> &knots)
      -> std::vector<std::vector<FT>>;

  /**
   * @brief Solve linear system for interpolation
   * @param matrix Coefficient matrix
   * @param rightHandSide RHS vectors (one per dimension)
   * @return Solution vectors
   */
  static auto
  solveLinearSystem(const std::vector<std::vector<FT>> &matrix,
                    const std::vector<std::vector<FT>> &rightHandSide)
      -> std::vector<std::vector<FT>>;
};

} // namespace SFCGAL

#endif // SFCGAL_NURBSCURVE_H
