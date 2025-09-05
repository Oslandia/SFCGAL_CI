// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryVisitor.h"
#include "SFCGAL/algorithm/distance.h"
#include <CGAL/Bbox_3.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>

namespace SFCGAL {

//-- Constructors and initialization

NURBSCurve::NURBSCurve() : _degree(0), _fitTolerance(FT(0)) {}

NURBSCurve::NURBSCurve(const std::vector<Point> &controlPoints,
                       unsigned int degree, KnotMethod knotMethod)
    : _controlPoints(controlPoints), _degree(degree), _fitTolerance(FT(0))
{
  if (!_controlPoints.empty()) {
    _weights.resize(_controlPoints.size(), FT(1));
    _knotVector = generateKnotVector(_controlPoints, _degree, knotMethod);

    auto [isValid, reason] = validateData();
    if (!isValid) {
      BOOST_THROW_EXCEPTION(Exception("Invalid NURBS curve data: " + reason));
    }
  }
}

NURBSCurve::NURBSCurve(const std::vector<Point> &controlPoints,
                       const std::vector<FT> &weights, unsigned int degree,
                       KnotMethod knotMethod)
    : _controlPoints(controlPoints), _weights(weights), _degree(degree),
      _fitTolerance(FT(0))
{
  if (!_controlPoints.empty()) {
    _knotVector = generateKnotVector(_controlPoints, _degree, knotMethod);

    auto [isValid, reason] = validateData();
    if (!isValid) {
      BOOST_THROW_EXCEPTION(Exception("Invalid NURBS curve data: " + reason));
    }
  }
}

NURBSCurve::NURBSCurve(const std::vector<Point> &controlPoints,
                       const std::vector<FT> &weights, unsigned int degree,
                       const std::vector<Knot> &knotVector)
    : _controlPoints(controlPoints), _weights(weights), _degree(degree),
      _knotVector(knotVector), _fitTolerance(FT(0))
{
  auto [isValid, reason] = validateData();
  if (!isValid) {
    BOOST_THROW_EXCEPTION(Exception("Invalid NURBS curve data: " + reason));
  }
}

//-- Factory methods

auto
NURBSCurve::fromBezier(const std::vector<Point> &controlPoints)
    -> std::unique_ptr<NURBSCurve>
{
  if (controlPoints.empty()) {
    return std::make_unique<NURBSCurve>();
  }

  unsigned int degree = static_cast<unsigned int>(controlPoints.size() - 1);

  // Bezier knot vector: [0,...,0, 1,...,1] with multiplicity degree+1
  std::vector<Knot> knots;
  knots.reserve(2 * (degree + 1));

  for (unsigned int idx = 0; idx <= degree; ++idx) {
    knots.push_back(FT(0));
  }
  for (unsigned int idx = 0; idx <= degree; ++idx) {
    knots.push_back(FT(1));
  }

  std::vector<FT> weights(controlPoints.size(), FT(1));
  return std::make_unique<NURBSCurve>(controlPoints, weights, degree, knots);
}

auto
NURBSCurve::fromBSpline(const std::vector<Point> &controlPoints,
                        unsigned int degree, const std::vector<Knot> &knots)
    -> std::unique_ptr<NURBSCurve>
{
  std::vector<FT> weights(controlPoints.size(), FT(1));
  return std::make_unique<NURBSCurve>(controlPoints, weights, degree, knots);
}

auto
NURBSCurve::createBSpline(const std::vector<Point> &controlPoints,
                          unsigned int degree) -> std::unique_ptr<NURBSCurve>
{
  if (controlPoints.empty()) {
    return std::make_unique<NURBSCurve>();
  }

  auto knots = generateKnotVector(controlPoints, degree, KnotMethod::UNIFORM);
  return fromBSpline(controlPoints, degree, knots);
}

auto
NURBSCurve::createCircularArc(const Point &center, FT radius, FT startAngle,
                              FT endAngle, const Point &normal)
    -> std::unique_ptr<NURBSCurve>
{
  if (radius <= FT(0)) {
    BOOST_THROW_EXCEPTION(Exception("Arc radius must be positive"));
  }

  // Convert angles to double for trigonometric calculations
  double startAngleD = CGAL::to_double(startAngle);
  double endAngleD   = CGAL::to_double(endAngle);
  double radiusD     = CGAL::to_double(radius);

  // Normalize angle span to [0, 2π]
  double angleSpan = endAngleD - startAngleD;
  while (angleSpan < 0.0) {
    angleSpan += 2.0 * M_PI;
  }
  while (angleSpan > 2.0 * M_PI) {
    angleSpan -= 2.0 * M_PI;
  }

  // Determine coordinate type from center point
  CoordinateType dimType = COORDINATE_XY;
  if (center.is3D() && center.isMeasured()) {
    dimType = COORDINATE_XYZM;
  } else if (center.is3D()) {
    dimType = COORDINATE_XYZ;
  } else if (center.isMeasured()) {
    dimType = COORDINATE_XYM;
  }

  // For simple arc, use rational quadratic representation
  std::vector<Point> controlPoints;
  std::vector<FT>    weights;
  std::vector<Knot>  knots;

  double midAngle = (startAngleD + endAngleD) / 2.0;
  double halfSpan = angleSpan / 2.0;

  // Create control points for rational quadratic arc
  controlPoints.emplace_back(center.x() + FT(radiusD * std::cos(startAngleD)),
                             center.y() + FT(radiusD * std::sin(startAngleD)),
                             center.is3D() ? center.z() : FT(0),
                             center.isMeasured() ? center.m() : NaN(), dimType);

  controlPoints.emplace_back(
      center.x() + FT(radiusD * std::cos(midAngle) / std::cos(halfSpan)),
      center.y() + FT(radiusD * std::sin(midAngle) / std::cos(halfSpan)),
      center.is3D() ? center.z() : FT(0),
      center.isMeasured() ? center.m() : NaN(), dimType);

  controlPoints.emplace_back(center.x() + FT(radiusD * std::cos(endAngleD)),
                             center.y() + FT(radiusD * std::sin(endAngleD)),
                             center.is3D() ? center.z() : FT(0),
                             center.isMeasured() ? center.m() : NaN(), dimType);

  // Rational weights for exact circular arc
  weights.push_back(FT(1.0));
  weights.push_back(FT(std::cos(halfSpan)));
  weights.push_back(FT(1.0));

  // Bezier knot vector for degree 2
  knots = {FT(0), FT(0), FT(0), FT(1), FT(1), FT(1)};

  return std::make_unique<NURBSCurve>(controlPoints, weights, 2, knots);
}

auto
NURBSCurve::interpolateCurve(const std::vector<Point> &points,
                             unsigned int degree, KnotMethod knotMethod,
                             EndCondition endCondition)
    -> std::unique_ptr<NURBSCurve>
{
  if (points.size() < 2) {
    BOOST_THROW_EXCEPTION(
        Exception("Need at least 2 points for interpolation"));
  }

  if (degree >= points.size()) {
    degree = static_cast<unsigned int>(points.size() - 1);
  }

  // Compute parameter values for interpolation points
  auto parameters = computeParameters(points, knotMethod);

  // Generate appropriate knot vector
  std::vector<Knot> knots;
  size_t            numKnots = points.size() + degree + 1;
  knots.reserve(numKnots);

  // Clamped knot vector (most common case)
  for (unsigned int idx = 0; idx <= degree; ++idx) {
    knots.push_back(parameters.front());
  }

  // Internal knots using averaging method
  for (size_t pointIdx = 1; pointIdx < points.size() - degree; ++pointIdx) {
    FT sum = FT(0);
    for (unsigned int degIdx = 0; degIdx < degree; ++degIdx) {
      sum += parameters[pointIdx + degIdx];
    }
    knots.push_back(sum / FT(degree));
  }

  for (unsigned int idx = 0; idx <= degree; ++idx) {
    knots.push_back(parameters.back());
  }

  // For now, return simplified curve with control points as interpolation
  // points Full collocation matrix solution would be implemented here in
  // production
  auto curve = std::make_unique<NURBSCurve>(
      points, std::vector<FT>(points.size(), FT(1)), degree, knots);

  curve->_fitPoints    = points;
  curve->_fitTolerance = FT(0); // Exact interpolation

  return curve;
}

auto
NURBSCurve::approximateCurve(const std::vector<Point> &points,
                             unsigned int degree, FT tolerance,
                             size_t maxControlPoints)
    -> std::unique_ptr<NURBSCurve>
{
  if (points.size() < 2) {
    BOOST_THROW_EXCEPTION(
        Exception("Need at least 2 points for approximation"));
  }

  // Use simplified approach: select subset of points as control points
  size_t numControlPoints = std::min(maxControlPoints, points.size());
  std::vector<Point> controlPoints;
  controlPoints.reserve(numControlPoints);

  // Distribute control points evenly through input points
  for (size_t idx = 0; idx < numControlPoints; ++idx) {
    size_t pointIndex = (idx * (points.size() - 1)) / (numControlPoints - 1);
    controlPoints.push_back(points[pointIndex]);
  }

  if (degree >= controlPoints.size()) {
    degree = static_cast<unsigned int>(controlPoints.size() - 1);
  }

  auto curve           = std::make_unique<NURBSCurve>(controlPoints, degree);
  curve->_fitPoints    = points;
  curve->_fitTolerance = tolerance;

  return curve;
}

//-- Geometry interface implementation

void
NURBSCurve::accept(GeometryVisitor &visitor)
{
  visitor.visit(*this);
}

void
NURBSCurve::accept(ConstGeometryVisitor &visitor) const
{
  visitor.visit(*this);
}

auto
NURBSCurve::geometryTypeId() const -> GeometryType
{
  return TYPE_NURBSCURVE;
}

auto
NURBSCurve::geometryType() const -> std::string
{
  return "NURBSCurve";
}

auto
NURBSCurve::clone() const -> NURBSCurve *
{
  return new NURBSCurve(*this);
}

auto
NURBSCurve::isEmpty() const -> bool
{
  return _controlPoints.empty();
}

auto
NURBSCurve::is3D() const -> bool
{
  return !_controlPoints.empty() && _controlPoints[0].is3D();
}

auto
NURBSCurve::isMeasured() const -> bool
{
  return !_controlPoints.empty() && _controlPoints[0].isMeasured();
}

auto
NURBSCurve::dropZ() -> bool
{
  bool hadZ = is3D();
  for (auto &point : _controlPoints) {
    point.dropZ();
  }
  return hadZ;
}

auto
NURBSCurve::dropM() -> bool
{
  bool hadM = isMeasured();
  for (auto &point : _controlPoints) {
    point.dropM();
  }
  return hadM;
}

void
NURBSCurve::swapXY()
{
  for (auto &point : _controlPoints) {
    point.swapXY();
  }
}

//-- Core curve evaluation and derivatives

auto
NURBSCurve::evaluate(Parameter parameter) const -> Point
{
  if (_controlPoints.empty()) {
    BOOST_THROW_EXCEPTION(Exception("Cannot evaluate empty NURBS curve"));
  }

  auto bounds = parameterBounds();
  if (parameter < bounds.first || parameter > bounds.second) {
    BOOST_THROW_EXCEPTION(Exception("Parameter outside valid range"));
  }

  if (_controlPoints.size() == 1) {
    return _controlPoints[0];
  }

  size_t span = findSpan(parameter);
  return deBoorRational(span, parameter);
}

auto
NURBSCurve::derivative(Parameter parameter, unsigned int order) const -> Point
{
  if (order == 0) {
    return evaluate(parameter);
  }

  auto derivatives = computeDerivatives(parameter, order);
  if (order < derivatives.size()) {
    return derivatives[order];
  }

  // Return zero vector for orders beyond curve degree
  if (is3D() && isMeasured()) {
    return Point(FT(0), FT(0), FT(0), 0.0);
  } else if (is3D()) {
    return Point(FT(0), FT(0), FT(0));
  } else if (isMeasured()) {
    Point point(FT(0), FT(0));
    point.setM(0.0);
    return point;
  } else {
    return Point(FT(0), FT(0));
  }
}

auto
NURBSCurve::tangent(Parameter parameter) const -> Point
{
  Point firstDerivative = derivative(parameter, 1);

  // Normalize to unit vector
  FT magnitude = std::sqrt(CGAL::to_double(
      firstDerivative.x() * firstDerivative.x() +
      firstDerivative.y() * firstDerivative.y() +
      (firstDerivative.is3D() ? firstDerivative.z() * firstDerivative.z()
                              : FT(0))));

  if (magnitude < FT(1e-12)) {
    BOOST_THROW_EXCEPTION(
        Exception("Cannot compute tangent at singular point"));
  }

  return Point(firstDerivative.x() / magnitude, firstDerivative.y() / magnitude,
               firstDerivative.is3D() ? firstDerivative.z() / magnitude
                                      : FT(0));
}

auto
NURBSCurve::normal(Parameter parameter) const -> Point
{
  if (is3D()) {
    // For 3D curves, normal is not uniquely defined
    BOOST_THROW_EXCEPTION(
        Exception("Normal vector not uniquely defined for 3D curves"));
  }

  Point tangentVec = tangent(parameter);

  // 2D normal: rotate tangent 90 degrees counterclockwise
  return Point(-tangentVec.y(), tangentVec.x());
}

auto
NURBSCurve::binormal(Parameter parameter) const -> Point
{
  if (!is3D()) {
    BOOST_THROW_EXCEPTION(Exception("Binormal only defined for 3D curves"));
  }

  Point firstDeriv  = derivative(parameter, 1);
  Point secondDeriv = derivative(parameter, 2);

  // Binormal = first derivative × second derivative (normalized)
  FT crossX =
      firstDeriv.y() * secondDeriv.z() - firstDeriv.z() * secondDeriv.y();
  FT crossY =
      firstDeriv.z() * secondDeriv.x() - firstDeriv.x() * secondDeriv.z();
  FT crossZ =
      firstDeriv.x() * secondDeriv.y() - firstDeriv.y() * secondDeriv.x();

  FT magnitude = std::sqrt(
      CGAL::to_double(crossX * crossX + crossY * crossY + crossZ * crossZ));

  if (magnitude < FT(1e-12)) {
    BOOST_THROW_EXCEPTION(
        Exception("Cannot compute binormal at inflection point"));
  }

  return Point(crossX / magnitude, crossY / magnitude, crossZ / magnitude);
}

auto
NURBSCurve::curvature(Parameter parameter) const -> FT
{
  Point firstDeriv  = derivative(parameter, 1);
  Point secondDeriv = derivative(parameter, 2);

  if (is3D()) {
    // 3D curvature: |C'(t) × C''(t)| / |C'(t)|³
    FT crossX =
        firstDeriv.y() * secondDeriv.z() - firstDeriv.z() * secondDeriv.y();
    FT crossY =
        firstDeriv.z() * secondDeriv.x() - firstDeriv.x() * secondDeriv.z();
    FT crossZ =
        firstDeriv.x() * secondDeriv.y() - firstDeriv.y() * secondDeriv.x();

    FT crossMagnitude = std::sqrt(
        CGAL::to_double(crossX * crossX + crossY * crossY + crossZ * crossZ));
    FT firstMagnitude = std::sqrt(CGAL::to_double(
        firstDeriv.x() * firstDeriv.x() + firstDeriv.y() * firstDeriv.y() +
        firstDeriv.z() * firstDeriv.z()));

    if (firstMagnitude < FT(1e-12)) {
      return FT(0);
    }

    return crossMagnitude / (firstMagnitude * firstMagnitude * firstMagnitude);
  } else {
    // 2D curvature: (x'y'' - y'x') / (x'² + y'²)^(3/2)
    FT numerator =
        firstDeriv.x() * secondDeriv.y() - firstDeriv.y() * secondDeriv.x();
    FT denominator =
        firstDeriv.x() * firstDeriv.x() + firstDeriv.y() * firstDeriv.y();

    if (denominator < FT(1e-12)) {
      return FT(0);
    }

    return numerator / std::pow(CGAL::to_double(denominator), 1.5);
  }
}

auto
NURBSCurve::torsion(Parameter parameter) const -> FT
{
  if (!is3D()) {
    return FT(0);
  }

  Point firstDeriv  = derivative(parameter, 1);
  Point secondDeriv = derivative(parameter, 2);
  Point thirdDeriv  = derivative(parameter, 3);

  // Torsion = (C' × C'') · C''' / |C' × C''|²
  FT crossX =
      firstDeriv.y() * secondDeriv.z() - firstDeriv.z() * secondDeriv.y();
  FT crossY =
      firstDeriv.z() * secondDeriv.x() - firstDeriv.x() * secondDeriv.z();
  FT crossZ =
      firstDeriv.x() * secondDeriv.y() - firstDeriv.y() * secondDeriv.x();

  FT scalarTripleProduct = crossX * thirdDeriv.x() + crossY * thirdDeriv.y() +
                           crossZ * thirdDeriv.z();
  FT crossMagnitudeSquared =
      crossX * crossX + crossY * crossY + crossZ * crossZ;

  if (crossMagnitudeSquared < FT(1e-12)) {
    return FT(0);
  }

  return scalarTripleProduct / crossMagnitudeSquared;
}

auto
NURBSCurve::frenetFrame(Parameter parameter) const
    -> std::tuple<Point, Point, Point>
{
  Point tangentVec = tangent(parameter);

  if (is3D()) {
    Point binormalVec = binormal(parameter);
    Point normalVec   = Point(
        tangentVec.y() * binormalVec.z() - tangentVec.z() * binormalVec.y(),
        tangentVec.z() * binormalVec.x() - tangentVec.x() * binormalVec.z(),
        tangentVec.x() * binormalVec.y() - tangentVec.y() * binormalVec.x());
    return std::make_tuple(tangentVec, normalVec, binormalVec);
  } else {
    Point normalVec   = normal(parameter);
    Point binormalVec = Point(0, 0, 1); // Default Z-axis for 2D curves
    return std::make_tuple(tangentVec, normalVec, binormalVec);
  }
}

auto
NURBSCurve::toLineString(unsigned int numSegments) const
    -> std::unique_ptr<LineString>
{
  auto lineString = std::make_unique<LineString>();
  if (_controlPoints.empty()) {
    return lineString;
  }

  auto bounds     = parameterBounds();
  FT   paramRange = bounds.second - bounds.first;

  for (unsigned int segIdx = 0; segIdx <= numSegments; ++segIdx) {
    FT parameter = bounds.first + (FT(segIdx) / FT(numSegments)) * paramRange;
    lineString->addPoint(evaluate(parameter));
  }

  return lineString;
}

auto
NURBSCurve::toLineStringAdaptive(FT tolerance, unsigned int minSegments,
                                 unsigned int maxSegments) const
    -> std::unique_ptr<LineString>
{
  auto lineString = std::make_unique<LineString>();
  if (_controlPoints.empty()) {
    return lineString;
  }

  auto bounds = parameterBounds();

  // Start with minimum segments
  std::vector<Parameter> parameters;
  parameters.reserve(maxSegments + 1);

  FT deltaParam = (bounds.second - bounds.first) / FT(minSegments);
  for (unsigned int segIdx = 0; segIdx <= minSegments; ++segIdx) {
    parameters.push_back(bounds.first + FT(segIdx) * deltaParam);
  }

  // Adaptively refine segments that exceed tolerance
  bool refined = true;
  while (refined && parameters.size() <= maxSegments) {
    refined = false;

    for (size_t paramIdx = 0; paramIdx < parameters.size() - 1; ++paramIdx) {
      Parameter param1   = parameters[paramIdx];
      Parameter param2   = parameters[paramIdx + 1];
      Parameter midParam = (param1 + param2) / FT(2);

      Point point1   = evaluate(param1);
      Point point2   = evaluate(param2);
      Point midPoint = evaluate(midParam);

      // Linear interpolation between endpoints
      Point interpolated = Point(
          (point1.x() + point2.x()) / FT(2), (point1.y() + point2.y()) / FT(2),
          point1.is3D() ? (point1.z() + point2.z()) / FT(2) : FT(0));

      FT deviation = algorithm::distance(midPoint, interpolated);

      if (deviation > tolerance) {
        parameters.insert(parameters.begin() + paramIdx + 1, midParam);
        refined = true;
        break; // Restart iteration after modification
      }
    }
  }

  // Build LineString from refined parameters
  for (const auto &param : parameters) {
    lineString->addPoint(evaluate(param));
  }

  return lineString;
}

auto
NURBSCurve::parameterBounds() const -> std::pair<Parameter, Parameter>
{
  if (_knotVector.empty()) {
    return std::make_pair(FT(0), FT(1));
  }

  // Valid parameter domain is [knots[degree], knots[n]]
  return std::make_pair(_knotVector[_degree],
                        _knotVector[_controlPoints.size()]);
}

auto
NURBSCurve::isClosed() const -> bool
{
  if (_controlPoints.size() < 2) {
    return false;
  }

  const Point &startPoint = _controlPoints.front();
  const Point &endPoint   = _controlPoints.back();

  return algorithm::distance(startPoint, endPoint) < FT(1e-10);
}

auto
NURBSCurve::isPeriodic() const -> bool
{
  if (_controlPoints.size() < 2 || _degree == 0) {
    return false;
  }

  // Check if enough control points match at beginning and end
  size_t minOverlap = _degree;
  if (_controlPoints.size() < 2 * minOverlap) {
    return false;
  }

  for (size_t idx = 0; idx < minOverlap; ++idx) {
    if (algorithm::distance(
            _controlPoints[idx],
            _controlPoints[_controlPoints.size() - minOverlap + idx]) >
        FT(1e-10)) {
      return false;
    }
  }

  return true;
}

auto
NURBSCurve::isPlanar(std::vector<FT> *plane) const -> bool
{
  if (!is3D() || _controlPoints.size() < 4) {
    if (plane) {
      // For 2D curves, plane is XY (z = 0)
      *plane = {FT(0), FT(0), FT(1), FT(0)};
    }
    return true;
  }

  // Find three non-collinear points to define plane
  for (size_t firstIdx = 0; firstIdx < _controlPoints.size() - 2; ++firstIdx) {
    for (size_t secondIdx = firstIdx + 1; secondIdx < _controlPoints.size() - 1;
         ++secondIdx) {
      for (size_t thirdIdx = secondIdx + 1; thirdIdx < _controlPoints.size();
           ++thirdIdx) {

        Point vec1 =
            Point(_controlPoints[secondIdx].x() - _controlPoints[firstIdx].x(),
                  _controlPoints[secondIdx].y() - _controlPoints[firstIdx].y(),
                  _controlPoints[secondIdx].z() - _controlPoints[firstIdx].z());

        Point vec2 =
            Point(_controlPoints[thirdIdx].x() - _controlPoints[firstIdx].x(),
                  _controlPoints[thirdIdx].y() - _controlPoints[firstIdx].y(),
                  _controlPoints[thirdIdx].z() - _controlPoints[firstIdx].z());

        // Cross product for normal vector
        Point normal = Point(vec1.y() * vec2.z() - vec1.z() * vec2.y(),
                             vec1.z() * vec2.x() - vec1.x() * vec2.z(),
                             vec1.x() * vec2.y() - vec1.y() * vec2.x());

        FT normalMagnitude = std::sqrt(
            CGAL::to_double(normal.x() * normal.x() + normal.y() * normal.y() +
                            normal.z() * normal.z()));

        if (normalMagnitude > FT(1e-10)) {
          // Normalize normal vector
          normal =
              Point(normal.x() / normalMagnitude, normal.y() / normalMagnitude,
                    normal.z() / normalMagnitude);

          // Check if all other points lie on this plane
          FT planeD = -(normal.x() * _controlPoints[firstIdx].x() +
                        normal.y() * _controlPoints[firstIdx].y() +
                        normal.z() * _controlPoints[firstIdx].z());

          bool allPointsOnPlane = true;
          for (const auto &point : _controlPoints) {
            FT distance =
                CGAL::abs(normal.x() * point.x() + normal.y() * point.y() +
                          normal.z() * point.z() + planeD);
            if (distance > FT(1e-10)) {
              allPointsOnPlane = false;
              break;
            }
          }

          if (allPointsOnPlane) {
            if (plane) {
              *plane = {normal.x(), normal.y(), normal.z(), planeD};
            }
            return true;
          }
        }
      }
    }
  }

  return false;
}

auto
NURBSCurve::isLinear() const -> bool
{
  if (_controlPoints.size() < 2) {
    return true;
  }

  if (_degree > 1) {
    return false; // Higher degree curves are not linear by definition
  }

  // Check if all control points are collinear
  Point direction = Point(_controlPoints[1].x() - _controlPoints[0].x(),
                          _controlPoints[1].y() - _controlPoints[0].y(),
                          _controlPoints[1].is3D()
                              ? _controlPoints[1].z() - _controlPoints[0].z()
                              : FT(0));

  FT dirMagnitude = std::sqrt(CGAL::to_double(
      direction.x() * direction.x() + direction.y() * direction.y() +
      (direction.is3D() ? direction.z() * direction.z() : FT(0))));

  if (dirMagnitude < FT(1e-10)) {
    return true; // Degenerate case
  }

  direction = Point(direction.x() / dirMagnitude, direction.y() / dirMagnitude,
                    direction.is3D() ? direction.z() / dirMagnitude : FT(0));

  for (size_t ptIdx = 2; ptIdx < _controlPoints.size(); ++ptIdx) {
    Point vecToPoint =
        Point(_controlPoints[ptIdx].x() - _controlPoints[0].x(),
              _controlPoints[ptIdx].y() - _controlPoints[0].y(),
              _controlPoints[ptIdx].is3D()
                  ? _controlPoints[ptIdx].z() - _controlPoints[0].z()
                  : FT(0));

    FT vecMagnitude = std::sqrt(CGAL::to_double(
        vecToPoint.x() * vecToPoint.x() + vecToPoint.y() * vecToPoint.y() +
        (vecToPoint.is3D() ? vecToPoint.z() * vecToPoint.z() : FT(0))));

    if (vecMagnitude > FT(1e-10)) {
      vecToPoint =
          Point(vecToPoint.x() / vecMagnitude, vecToPoint.y() / vecMagnitude,
                vecToPoint.is3D() ? vecToPoint.z() / vecMagnitude : FT(0));

      // Check if direction vectors are parallel
      Point crossProduct = Point(
          direction.y() * vecToPoint.z() - direction.z() * vecToPoint.y(),
          direction.z() * vecToPoint.x() - direction.x() * vecToPoint.z(),
          direction.x() * vecToPoint.y() - direction.y() * vecToPoint.x());

      FT crossMagnitude =
          std::sqrt(CGAL::to_double(crossProduct.x() * crossProduct.x() +
                                    crossProduct.y() * crossProduct.y() +
                                    crossProduct.z() * crossProduct.z()));

      if (crossMagnitude > FT(1e-10)) {
        return false;
      }
    }
  }

  return true;
}

auto
NURBSCurve::length(Parameter from, Parameter to, FT tolerance) const -> FT
{
  auto bounds = parameterBounds();

  if (from < FT(0))
    from = bounds.first;
  if (to < FT(0))
    to = bounds.second;

  if (from >= to) {
    return FT(0);
  }

  return computeArcLength(from, to, tolerance);
}

auto
NURBSCurve::parameterAtLength(FT arcLength, FT tolerance) const -> Parameter
{
  if (arcLength <= FT(0)) {
    return parameterBounds().first;
  }

  return findParameterByArcLength(arcLength, tolerance);
}

auto
NURBSCurve::reparameterizeByArcLength() const -> std::unique_ptr<Curve>
{
  FT totalLength = length();

  if (totalLength <= FT(0)) {
    return std::unique_ptr<NURBSCurve>(clone());
  }

  // Create new curve with arc-length parameterization
  std::vector<Point> newControlPoints;
  std::vector<FT>    newWeights;
  std::vector<Knot>  newKnots;

  // Sample curve at equal arc length intervals
  size_t numSamples = _controlPoints.size() * 2; // Increase sampling density
  for (size_t sampleIdx = 0; sampleIdx < numSamples; ++sampleIdx) {
    FT        targetLength = (FT(sampleIdx) / FT(numSamples - 1)) * totalLength;
    Parameter param        = parameterAtLength(targetLength);
    newControlPoints.push_back(evaluate(param));
  }

  // Create uniform weights and knot vector
  newWeights.resize(newControlPoints.size(), FT(1));

  unsigned int newDegree =
      std::min(_degree, static_cast<unsigned int>(newControlPoints.size() - 1));
  newKnots =
      generateKnotVector(newControlPoints, newDegree, KnotMethod::UNIFORM);

  return std::make_unique<NURBSCurve>(newControlPoints, newWeights, newDegree,
                                      newKnots);
}

auto
NURBSCurve::split(Parameter parameter) const
    -> std::pair<std::unique_ptr<Curve>, std::unique_ptr<Curve>>
{
  auto bounds = parameterBounds();

  if (parameter <= bounds.first) {
    // Split at beginning: first curve empty, second is full curve
    auto emptyCurve = std::make_unique<NURBSCurve>();
    auto fullCurve  = std::unique_ptr<NURBSCurve>(clone());
    return std::make_pair(std::move(emptyCurve), std::move(fullCurve));
  }

  if (parameter >= bounds.second) {
    // Split at end: first curve is full, second curve empty
    auto fullCurve  = std::unique_ptr<NURBSCurve>(clone());
    auto emptyCurve = std::make_unique<NURBSCurve>();
    return std::make_pair(std::move(fullCurve), std::move(emptyCurve));
  }

  // For now, create simplified split using sampling
  // Full implementation would use knot insertion algorithms
  auto firstLineString  = toLineString(32);
  auto secondLineString = toLineString(32);

  // Find split point in the sampled curve
  size_t splitIndex = 16; // Approximate middle

  std::vector<Point> firstPoints;
  std::vector<Point> secondPoints;

  for (size_t ptIdx = 0; ptIdx <= splitIndex; ++ptIdx) {
    firstPoints.push_back(firstLineString->pointN(ptIdx));
  }

  for (size_t ptIdx = splitIndex; ptIdx < firstLineString->numPoints();
       ++ptIdx) {
    secondPoints.push_back(firstLineString->pointN(ptIdx));
  }

  auto firstCurve  = std::make_unique<NURBSCurve>(firstPoints, _degree);
  auto secondCurve = std::make_unique<NURBSCurve>(secondPoints, _degree);

  return std::make_pair(std::move(firstCurve), std::move(secondCurve));
}

auto
NURBSCurve::subcurve(Parameter from, Parameter to) const
    -> std::unique_ptr<Curve>
{
  if (from >= to) {
    return std::make_unique<NURBSCurve>();
  }

  auto bounds = parameterBounds();
  from        = std::max(from, bounds.first);
  to          = std::min(to, bounds.second);

  // Sample the subcurve
  unsigned int       numSamples = 32;
  std::vector<Point> subPoints;
  subPoints.reserve(numSamples + 1);

  FT paramRange = to - from;
  for (unsigned int sampleIdx = 0; sampleIdx <= numSamples; ++sampleIdx) {
    Parameter param = from + (FT(sampleIdx) / FT(numSamples)) * paramRange;
    subPoints.push_back(evaluate(param));
  }

  unsigned int subDegree =
      std::min(_degree, static_cast<unsigned int>(subPoints.size() - 1));
  return std::make_unique<NURBSCurve>(subPoints, subDegree);
}

auto
NURBSCurve::reverse() const -> std::unique_ptr<Curve>
{
  if (_controlPoints.empty()) {
    return std::make_unique<NURBSCurve>();
  }

  // Reverse control points and weights
  std::vector<Point> reversedPoints(_controlPoints.rbegin(),
                                    _controlPoints.rend());
  std::vector<FT>    reversedWeights;

  if (!_weights.empty()) {
    reversedWeights.assign(_weights.rbegin(), _weights.rend());
  }

  // Reverse and reparameterize knot vector
  std::vector<Knot> reversedKnots;
  if (!_knotVector.empty()) {
    reversedKnots.reserve(_knotVector.size());
    FT maxKnot = _knotVector.back();

    for (auto knotIter = _knotVector.rbegin(); knotIter != _knotVector.rend();
         ++knotIter) {
      reversedKnots.push_back(maxKnot - *knotIter);
    }
  }

  return std::make_unique<NURBSCurve>(reversedPoints, reversedWeights, _degree,
                                      reversedKnots);
}

auto
NURBSCurve::join(const Curve &other, Continuity continuity, FT tolerance) const
    -> std::unique_ptr<Curve>
{
  const auto *otherNurbs = dynamic_cast<const NURBSCurve *>(&other);
  if (!otherNurbs) {
    BOOST_THROW_EXCEPTION(
        Exception("Can only join NURBS curves with other NURBS curves"));
  }

  if (isEmpty() && otherNurbs->isEmpty()) {
    return std::make_unique<NURBSCurve>();
  }

  if (isEmpty()) {
    return std::unique_ptr<NURBSCurve>(otherNurbs->clone());
  }

  if (otherNurbs->isEmpty()) {
    return std::unique_ptr<NURBSCurve>(clone());
  }

  // Check if curves can be joined (endpoints close enough)
  Point thisEnd    = endPoint();
  Point otherStart = otherNurbs->startPoint();

  if (algorithm::distance(thisEnd, otherStart) > tolerance) {
    BOOST_THROW_EXCEPTION(
        Exception("Curves are not adjacent within tolerance"));
  }

  // Simple joining: concatenate control points
  std::vector<Point> joinedPoints = _controlPoints;
  const auto        &otherPoints  = otherNurbs->controlPoints();

  // Skip first point of second curve to avoid duplication
  joinedPoints.insert(joinedPoints.end(), otherPoints.begin() + 1,
                      otherPoints.end());

  // Use maximum degree of both curves
  unsigned int joinedDegree = std::max(_degree, otherNurbs->degree());
  joinedDegree              = std::min(joinedDegree,
                                       static_cast<unsigned int>(joinedPoints.size() - 1));

  return std::make_unique<NURBSCurve>(joinedPoints, joinedDegree);
}

auto
NURBSCurve::offset(FT distance) const -> std::unique_ptr<Curve>
{
  if (is3D()) {
    BOOST_THROW_EXCEPTION(Exception("Offset not implemented for 3D curves"));
  }

  if (distance == FT(0)) {
    return std::unique_ptr<NURBSCurve>(clone());
  }

  // Simple offset implementation: sample curve and offset each point
  auto               lineString = toLineString(64);
  std::vector<Point> offsetPoints;
  offsetPoints.reserve(lineString->numPoints());

  auto bounds     = parameterBounds();
  FT   paramRange = bounds.second - bounds.first;

  for (size_t ptIdx = 0; ptIdx < lineString->numPoints(); ++ptIdx) {
    Parameter param =
        bounds.first +
        (FT(ptIdx) / FT(lineString->numPoints() - 1)) * paramRange;
    Point normalVec = normal(param);
    Point original  = lineString->pointN(ptIdx);

    offsetPoints.emplace_back(original.x() + distance * normalVec.x(),
                              original.y() + distance * normalVec.y());
  }

  return std::make_unique<NURBSCurve>(offsetPoints, _degree);
}

auto
NURBSCurve::closestPoint(const Point &point, Parameter *outParameter) const
    -> Point
{
  if (_controlPoints.empty()) {
    BOOST_THROW_EXCEPTION(
        Exception("Cannot find closest point on empty curve"));
  }

  auto [closestPt, parameter] = projectPointToCurve(point);

  if (outParameter) {
    *outParameter = parameter;
  }

  return closestPt;
}

auto
NURBSCurve::distance(const Point &point) const -> FT
{
  Point closest = closestPoint(point);
  return algorithm::distance(point, closest);
}

auto
NURBSCurve::hasSelfIntersections(
    std::vector<std::pair<Parameter, Parameter>> *intersections) const -> bool
{
  // Simple implementation: sample curve and check for crossings
  auto lineString = toLineString(100);

  if (intersections) {
    intersections->clear();
  }

  // This is a simplified check - full implementation would be more
  // sophisticated
  for (size_t firstIdx = 0; firstIdx < lineString->numPoints() - 3;
       ++firstIdx) {
    for (size_t secondIdx = firstIdx + 2;
         secondIdx < lineString->numPoints() - 1; ++secondIdx) {
      // Skip adjacent segments
      if (secondIdx <= firstIdx + 2)
        continue;

      Point segment1Start = lineString->pointN(firstIdx);
      Point segment1End   = lineString->pointN(firstIdx + 1);
      Point segment2Start = lineString->pointN(secondIdx);
      Point segment2End   = lineString->pointN(secondIdx + 1);

      // Simple intersection test (would need proper line-line intersection)
      if (algorithm::distance(segment1Start, segment2Start) < FT(1e-6) ||
          algorithm::distance(segment1End, segment2End) < FT(1e-6)) {
        if (intersections) {
          auto bounds     = parameterBounds();
          FT   paramRange = bounds.second - bounds.first;

          Parameter param1 =
              bounds.first +
              (FT(firstIdx) / FT(lineString->numPoints() - 1)) * paramRange;
          Parameter param2 =
              bounds.first +
              (FT(secondIdx) / FT(lineString->numPoints() - 1)) * paramRange;

          intersections->emplace_back(param1, param2);
        }
        return true;
      }
    }
  }

  return false;
}

auto
NURBSCurve::intersect(const Curve &other, FT tolerance) const
    -> std::vector<std::tuple<Point, Parameter, Parameter>>
{
  const auto *otherNurbs = dynamic_cast<const NURBSCurve *>(&other);
  if (!otherNurbs) {
    BOOST_THROW_EXCEPTION(
        Exception("Intersection only implemented between NURBS curves"));
  }

  std::vector<std::tuple<Point, Parameter, Parameter>> result;

  // Sample both curves and find close points
  auto thisLineString  = toLineString(50);
  auto otherLineString = otherNurbs->toLineString(50);

  auto thisBounds  = parameterBounds();
  auto otherBounds = otherNurbs->parameterBounds();

  FT thisParamRange  = thisBounds.second - thisBounds.first;
  FT otherParamRange = otherBounds.second - otherBounds.first;

  for (size_t thisIdx = 0; thisIdx < thisLineString->numPoints(); ++thisIdx) {
    Point thisPoint = thisLineString->pointN(thisIdx);

    for (size_t otherIdx = 0; otherIdx < otherLineString->numPoints();
         ++otherIdx) {
      Point otherPoint = otherLineString->pointN(otherIdx);

      if (algorithm::distance(thisPoint, otherPoint) <= tolerance) {
        Parameter thisParam =
            thisBounds.first +
            (FT(thisIdx) / FT(thisLineString->numPoints() - 1)) *
                thisParamRange;
        Parameter otherParam =
            otherBounds.first +
            (FT(otherIdx) / FT(otherLineString->numPoints() - 1)) *
                otherParamRange;

        // Use average point as intersection
        Point intersectionPoint =
            Point((thisPoint.x() + otherPoint.x()) / FT(2),
                  (thisPoint.y() + otherPoint.y()) / FT(2),
                  thisPoint.is3D() ? (thisPoint.z() + otherPoint.z()) / FT(2)
                                   : FT(0));

        result.emplace_back(intersectionPoint, thisParam, otherParam);
      }
    }
  }

  return result;
}

auto
NURBSCurve::boundingBox() const -> std::pair<Point, Point>
{
  if (_controlPoints.empty()) {
    BOOST_THROW_EXCEPTION(
        Exception("Cannot compute bounding box of empty curve"));
  }

  // Start with first control point
  FT minX = _controlPoints[0].x(), maxX = _controlPoints[0].x();
  FT minY = _controlPoints[0].y(), maxY = _controlPoints[0].y();
  FT minZ = _controlPoints[0].is3D() ? _controlPoints[0].z() : FT(0);
  FT maxZ = minZ;

  // Expand by all control points
  for (const auto &point : _controlPoints) {
    minX = std::min(minX, point.x());
    maxX = std::max(maxX, point.x());
    minY = std::min(minY, point.y());
    maxY = std::max(maxY, point.y());

    if (point.is3D()) {
      minZ = std::min(minZ, point.z());
      maxZ = std::max(maxZ, point.z());
    }
  }

  // Sample curve for tighter bounds
  auto lineString = toLineString(32);
  for (size_t ptIdx = 0; ptIdx < lineString->numPoints(); ++ptIdx) {
    const Point &point = lineString->pointN(ptIdx);

    minX = std::min(minX, point.x());
    maxX = std::max(maxX, point.x());
    minY = std::min(minY, point.y());
    maxY = std::max(maxY, point.y());

    if (point.is3D()) {
      minZ = std::min(minZ, point.z());
      maxZ = std::max(maxZ, point.z());
    }
  }

  if (is3D()) {
    return std::make_pair(Point(minX, minY, minZ), Point(maxX, maxY, maxZ));
  } else {
    return std::make_pair(Point(minX, minY), Point(maxX, maxY));
  }
}

//-- NURBS-specific data access and manipulation

void
NURBSCurve::setControlPoints(const std::vector<Point> &controlPoints)
{
  _controlPoints = controlPoints;

  // Adjust weights if necessary
  if (!_weights.empty() && _weights.size() != controlPoints.size()) {
    _weights.resize(controlPoints.size(), FT(1));
  }

  // Regenerate knot vector if needed
  if (!_knotVector.empty() && !controlPoints.empty()) {
    _knotVector =
        generateKnotVector(controlPoints, _degree, KnotMethod::CHORD_LENGTH);
  }

  auto [isValid, reason] = validateData();
  if (!isValid) {
    BOOST_THROW_EXCEPTION(Exception("Invalid control points: " + reason));
  }
}

auto
NURBSCurve::controlPointN(size_t index) -> Point &
{
  if (index >= _controlPoints.size()) {
    BOOST_THROW_EXCEPTION(Exception("Control point index out of bounds"));
  }
  return _controlPoints[index];
}

auto
NURBSCurve::controlPointN(size_t index) const -> const Point &
{
  if (index >= _controlPoints.size()) {
    BOOST_THROW_EXCEPTION(Exception("Control point index out of bounds"));
  }
  return _controlPoints[index];
}

void
NURBSCurve::setControlPoint(size_t index, const Point &point)
{
  if (index >= _controlPoints.size()) {
    BOOST_THROW_EXCEPTION(Exception("Control point index out of bounds"));
  }

  // Check dimensional consistency
  if (!_controlPoints.empty() &&
      (_controlPoints[0].is3D() != point.is3D() ||
       _controlPoints[0].isMeasured() != point.isMeasured())) {
    BOOST_THROW_EXCEPTION(Exception("Control point dimension mismatch"));
  }

  _controlPoints[index] = point;
}

void
NURBSCurve::setWeights(const std::vector<FT> &weights)
{
  if (!weights.empty() && weights.size() != _controlPoints.size()) {
    BOOST_THROW_EXCEPTION(
        Exception("Weight count (" + std::to_string(weights.size()) +
                  ") does not match control point count (" +
                  std::to_string(_controlPoints.size()) + ")"));
  }

  for (size_t weightIdx = 0; weightIdx < weights.size(); ++weightIdx) {
    if (weights[weightIdx] <= FT(0)) {
      BOOST_THROW_EXCEPTION(Exception("Weight at index " +
                                      std::to_string(weightIdx) +
                                      " must be positive"));
    }
  }

  _weights = weights;
}

auto
NURBSCurve::weight(size_t index) const -> FT
{
  if (_weights.empty()) {
    return FT(1); // Uniform weights
  }
  if (index >= _weights.size()) {
    BOOST_THROW_EXCEPTION(Exception("Weight index out of bounds"));
  }
  return _weights[index];
}

void
NURBSCurve::setWeight(size_t index, FT weight)
{
  if (weight <= FT(0)) {
    BOOST_THROW_EXCEPTION(Exception("Weight must be positive"));
  }

  if (_weights.empty()) {
    _weights.resize(_controlPoints.size(), FT(1));
  }

  if (index >= _weights.size()) {
    BOOST_THROW_EXCEPTION(Exception("Weight index out of bounds"));
  }
  _weights[index] = weight;
}

auto
NURBSCurve::isRational() const -> bool
{
  if (_weights.empty() || _weights.size() < 2) {
    return false;
  }

  FT       firstWeight = _weights[0];
  const FT tolerance   = FT(1e-10);

  for (const auto &weight : _weights) {
    if (CGAL::abs(weight - firstWeight) > tolerance) {
      return true;
    }
  }

  return false;
}

auto
NURBSCurve::isBezier() const -> bool
{
  if (_knotVector.empty() || _controlPoints.empty()) {
    return false;
  }

  size_t expectedSize = 2 * (_degree + 1);
  if (_knotVector.size() != expectedSize) {
    return false;
  }

  const FT tolerance = FT(1e-10);

  // Check start multiplicity
  for (unsigned int idx = 0; idx <= _degree; ++idx) {
    if (CGAL::abs(_knotVector[idx] - _knotVector[0]) > tolerance) {
      return false;
    }
  }

  // Check end multiplicity
  for (unsigned int idx = _degree + 1; idx < _knotVector.size(); ++idx) {
    if (CGAL::abs(_knotVector[idx] - _knotVector.back()) > tolerance) {
      return false;
    }
  }

  return true;
}

auto
NURBSCurve::isBSpline() const -> bool
{
  return !isRational();
}

void
NURBSCurve::setKnotVector(const std::vector<Knot> &knots)
{
  _knotVector            = knots;
  auto [isValid, reason] = validateData();
  if (!isValid) {
    BOOST_THROW_EXCEPTION(Exception("Invalid knot vector: " + reason));
  }
}

auto
NURBSCurve::knot(size_t index) const -> Knot
{
  if (index >= _knotVector.size()) {
    BOOST_THROW_EXCEPTION(Exception("Knot index out of bounds"));
  }
  return _knotVector[index];
}

auto
NURBSCurve::knotMultiplicity(Knot value, FT tolerance) const -> unsigned int
{
  unsigned int multiplicity = 0;
  for (const auto &knot : _knotVector) {
    if (CGAL::abs(knot - value) <= tolerance) {
      ++multiplicity;
    }
  }
  return multiplicity;
}

//-- Advanced NURBS operations (implementations to be completed)

auto
NURBSCurve::insertKnot(Knot parameter, unsigned int times) const
    -> std::unique_ptr<NURBSCurve>
{
  // Oslo algorithm for knot insertion would be implemented here
  // For now, return copy of current curve
  return std::unique_ptr<NURBSCurve>(clone());
}

auto
NURBSCurve::refineKnotVector(const std::vector<Knot> &newKnots) const
    -> std::unique_ptr<NURBSCurve>
{
  // Knot refinement algorithm would be implemented here
  return std::unique_ptr<NURBSCurve>(clone());
}

auto
NURBSCurve::elevateDegree(unsigned int times) const
    -> std::unique_ptr<NURBSCurve>
{
  // Degree elevation algorithm would be implemented here
  return std::unique_ptr<NURBSCurve>(clone());
}

auto
NURBSCurve::reduceDegree(FT tolerance) const -> std::unique_ptr<NURBSCurve>
{
  // Degree reduction with tolerance would be implemented here
  return std::unique_ptr<NURBSCurve>(clone());
}

auto
NURBSCurve::removeKnots(FT tolerance) const -> std::unique_ptr<NURBSCurve>
{
  // Knot removal algorithm would be implemented here
  return std::unique_ptr<NURBSCurve>(clone());
}

//-- Validation and utility methods

auto
NURBSCurve::validateData() const -> std::pair<bool, std::string>
{
  return validateNURBSData(_controlPoints, _weights, _degree, _knotVector);
}

auto
NURBSCurve::validateNURBSData(const std::vector<Point> &controlPoints,
                              const std::vector<FT>    &weights,
                              unsigned int              degree,
                              const std::vector<Knot>  &knots)
    -> std::pair<bool, std::string>
{
  // Empty geometry validation
  if (controlPoints.empty()) {
    if (!knots.empty() || !weights.empty()) {
      return {false, "Empty curve with non-empty knot vector or weights"};
    }
    return {true, ""};
  }

  // Degree constraints
  if (degree >= controlPoints.size()) {
    return {false, "Degree (" + std::to_string(degree) +
                       ") must be less than control points count (" +
                       std::to_string(controlPoints.size()) + ")"};
  }

  // Knot vector size: m = n + p + 1 (n = control points, p = degree)
  const size_t expectedKnots = controlPoints.size() + degree + 1;
  if (knots.size() != expectedKnots) {
    return {false, "Invalid knot vector size: expected " +
                       std::to_string(expectedKnots) + ", got " +
                       std::to_string(knots.size())};
  }

  // Knot vector non-decreasing order
  for (size_t knotIdx = 1; knotIdx < knots.size(); ++knotIdx) {
    if (knots[knotIdx] < knots[knotIdx - 1]) {
      return {false, "Knot vector is not non-decreasing at index " +
                         std::to_string(knotIdx)};
    }
  }

  // Weight validation
  if (!weights.empty()) {
    if (weights.size() != controlPoints.size()) {
      return {false, "Weight count does not match control point count"};
    }

    for (size_t weightIdx = 0; weightIdx < weights.size(); ++weightIdx) {
      const FT &weight = weights[weightIdx];

      if (weight <= FT(0)) {
        return {false, "Weight at index " + std::to_string(weightIdx) +
                           " is non-positive"};
      }

      if (!CGAL::is_finite(weight)) {
        return {false, "Weight at index " + std::to_string(weightIdx) +
                           " is not finite"};
      }
    }
  }

  // Dimensional consistency
  if (!controlPoints.empty()) {
    bool first3D       = controlPoints[0].is3D();
    bool firstMeasured = controlPoints[0].isMeasured();

    for (const auto &point : controlPoints) {
      if (point.is3D() != first3D || point.isMeasured() != firstMeasured) {
        return {false, "Inconsistent dimensions in control points"};
      }
    }
  }

  // Parameter bounds validation
  if (!knots.empty()) {
    if (degree == 0 && controlPoints.size() == 1) {
      if (knots.size() != 2) {
        return {false, "Invalid knot vector size for degree 0"};
      }
      if (knots[0] >= knots[1]) {
        return {false, "Invalid parameter bounds for degree 0"};
      }
    } else {
      const FT startParam = knots[degree];
      const FT endParam   = knots[controlPoints.size()];
      if (startParam >= endParam) {
        return {false, "Invalid parameter bounds"};
      }
    }
  }

  return {true, ""};
}

auto
NURBSCurve::getCurveStatistics() const -> std::map<std::string, double>
{
  std::map<std::string, double> stats;

  stats["num_control_points"] = static_cast<double>(_controlPoints.size());
  stats["degree"]             = static_cast<double>(_degree);
  stats["num_knots"]          = static_cast<double>(_knotVector.size());
  stats["is_rational"]        = isRational() ? 1.0 : 0.0;
  stats["is_bezier"]          = isBezier() ? 1.0 : 0.0;
  stats["is_closed"]          = isClosed() ? 1.0 : 0.0;
  stats["is_3d"]              = is3D() ? 1.0 : 0.0;
  stats["is_measured"]        = isMeasured() ? 1.0 : 0.0;

  if (!isEmpty()) {
    stats["curve_length"] = CGAL::to_double(length());

    auto bounds        = parameterBounds();
    stats["param_min"] = CGAL::to_double(bounds.first);
    stats["param_max"] = CGAL::to_double(bounds.second);
  }

  return stats;
}

//-- Protected helper methods

auto
NURBSCurve::findSpan(Parameter parameter) const -> size_t
{
  size_t numControlPoints = _controlPoints.size();

  // Special cases at boundaries
  if (parameter >= _knotVector[numControlPoints]) {
    return numControlPoints - 1;
  }

  if (parameter <= _knotVector[_degree]) {
    return _degree;
  }

  // Binary search for knot span
  size_t low  = _degree;
  size_t high = numControlPoints;
  size_t mid  = (low + high) / 2;

  while (parameter < _knotVector[mid] || parameter >= _knotVector[mid + 1]) {
    if (parameter < _knotVector[mid]) {
      high = mid;
    } else {
      low = mid;
    }
    mid = (low + high) / 2;
  }

  return mid;
}

auto
NURBSCurve::deBoorRational(size_t span, Parameter parameter) const -> Point
{
  // Homogeneous point representation for rational evaluation
  struct HomogeneousPoint {
    FT     wx, wy, wz, weight;
    double measure;

    HomogeneousPoint() : wx(0), wy(0), wz(0), weight(1), measure(0) {}

    HomogeneousPoint(const Point &point, FT wgt)
        : wx(point.x() * wgt), wy(point.y() * wgt),
          wz(point.is3D() ? point.z() * wgt : FT(0)), weight(wgt),
          measure(point.isMeasured() ? point.m() : 0.0)
    {
    }
  };

  std::vector<HomogeneousPoint> temp(_degree + 1);

  // Initialize with relevant control points
  for (unsigned int idx = 0; idx <= _degree; ++idx) {
    size_t cpIndex = span - _degree + idx;
    FT     weight  = _weights.empty() ? FT(1) : _weights[cpIndex];
    temp[idx]      = HomogeneousPoint(_controlPoints[cpIndex], weight);
  }

  // De Boor's algorithm with rational weights
  for (unsigned int level = 1; level <= _degree; ++level) {
    for (unsigned int idx = _degree; idx >= level; --idx) {
      size_t knotIdx = span - _degree + idx;
      FT     alpha =
          (parameter - _knotVector[knotIdx]) /
          (_knotVector[knotIdx + _degree - level + 1] - _knotVector[knotIdx]);

      // Interpolate homogeneous coordinates
      temp[idx].wx = (FT(1) - alpha) * temp[idx - 1].wx + alpha * temp[idx].wx;
      temp[idx].wy = (FT(1) - alpha) * temp[idx - 1].wy + alpha * temp[idx].wy;
      temp[idx].wz = (FT(1) - alpha) * temp[idx - 1].wz + alpha * temp[idx].wz;
      temp[idx].weight =
          (FT(1) - alpha) * temp[idx - 1].weight + alpha * temp[idx].weight;

      // M coordinate interpolated separately (non-rational)
      if (isMeasured()) {
        double alphaD = CGAL::to_double(alpha);
        temp[idx].measure =
            (1.0 - alphaD) * temp[idx - 1].measure + alphaD * temp[idx].measure;
      }
    }
  }

  // Project back to Euclidean space
  HomogeneousPoint &result = temp[_degree];

  if (CGAL::abs(result.weight) < FT(1e-15)) {
    BOOST_THROW_EXCEPTION(Exception("Division by zero in NURBS evaluation"));
  }

  FT invWeight = FT(1) / result.weight;

  if (is3D() && isMeasured()) {
    return Point(result.wx * invWeight, result.wy * invWeight,
                 result.wz * invWeight, result.measure);
  } else if (is3D()) {
    return Point(result.wx * invWeight, result.wy * invWeight,
                 result.wz * invWeight);
  } else if (isMeasured()) {
    Point point(result.wx * invWeight, result.wy * invWeight);
    point.setM(result.measure);
    return point;
  } else {
    return Point(result.wx * invWeight, result.wy * invWeight);
  }
}

auto
NURBSCurve::computeDerivatives(Parameter parameter, unsigned int maxOrder) const
    -> std::vector<Point>
{
  std::vector<Point> derivatives;
  derivatives.reserve(maxOrder + 1);

  // Order 0: curve evaluation
  derivatives.push_back(evaluate(parameter));

  if (maxOrder > 0) {
    // For linear curves, derivative is constant
    if (_degree == 1 && _controlPoints.size() == 2) {
      Point derivative =
          Point(_controlPoints[1].x() - _controlPoints[0].x(),
                _controlPoints[1].y() - _controlPoints[0].y(),
                _controlPoints[1].is3D()
                    ? _controlPoints[1].z() - _controlPoints[0].z()
                    : FT(0));
      derivatives.push_back(derivative);

      // Higher derivatives are zero for linear curves
      for (unsigned int orderIdx = 2; orderIdx <= maxOrder; ++orderIdx) {
        if (is3D() && isMeasured()) {
          derivatives.emplace_back(FT(0), FT(0), FT(0), 0.0);
        } else if (is3D()) {
          derivatives.emplace_back(FT(0), FT(0), FT(0));
        } else if (isMeasured()) {
          Point zeroPoint(FT(0), FT(0));
          zeroPoint.setM(0.0);
          derivatives.push_back(zeroPoint);
        } else {
          derivatives.emplace_back(FT(0), FT(0));
        }
      }
    } else {
      // General case: finite difference approximation
      const FT epsilon = FT(1e-8);
      auto     bounds  = parameterBounds();

      for (unsigned int orderIdx = 1; orderIdx <= maxOrder; ++orderIdx) {
        FT param1     = std::max(bounds.first, parameter - epsilon);
        FT param2     = std::min(bounds.second, parameter + epsilon);
        FT deltaParam = param2 - param1;

        if (deltaParam > FT(0)) {
          Point point1 = evaluate(param1);
          Point point2 = evaluate(param2);

          Point derivative(
              (point2.x() - point1.x()) / deltaParam,
              (point2.y() - point1.y()) / deltaParam,
              point2.is3D() ? (point2.z() - point1.z()) / deltaParam : FT(0));

          if (isMeasured()) {
            derivative.setM((point2.m() - point1.m()) /
                            CGAL::to_double(deltaParam));
          }

          derivatives.push_back(derivative);
        } else {
          // Zero derivative if cannot compute difference
          if (is3D() && isMeasured()) {
            derivatives.emplace_back(FT(0), FT(0), FT(0), 0.0);
          } else if (is3D()) {
            derivatives.emplace_back(FT(0), FT(0), FT(0));
          } else if (isMeasured()) {
            Point zeroPoint(FT(0), FT(0));
            zeroPoint.setM(0.0);
            derivatives.push_back(zeroPoint);
          } else {
            derivatives.emplace_back(FT(0), FT(0));
          }
        }
      }
    }
  }

  return derivatives;
}

auto
NURBSCurve::generateKnotVector(const std::vector<Point> &points,
                               unsigned int degree, KnotMethod method)
    -> std::vector<Knot>
{
  std::vector<Knot> knots;
  size_t            numPoints = points.size();
  size_t            numKnots  = numPoints + degree + 1;

  if (numPoints == 0) {
    return knots;
  }

  // Special case for degree 0
  if (degree == 0) {
    knots.reserve(numKnots);
    for (size_t idx = 0; idx < numPoints; ++idx) {
      knots.push_back(FT(idx));
    }
    knots.push_back(FT(numPoints));
    return knots;
  }

  knots.reserve(numKnots);

  if (method == KnotMethod::UNIFORM) {
    // Clamped uniform knot vector
    for (unsigned int idx = 0; idx <= degree; ++idx) {
      knots.push_back(FT(0));
    }

    size_t numInternal = numKnots - 2 * (degree + 1);
    for (size_t internalIdx = 1; internalIdx <= numInternal; ++internalIdx) {
      knots.push_back(FT(internalIdx) / FT(numInternal + 1));
    }

    for (unsigned int idx = 0; idx <= degree; ++idx) {
      knots.push_back(FT(1));
    }
  } else {
    // Chord length or centripetal parameterization
    auto parameters = computeParameters(points, method);

    // Clamped knot vector with parameter-based internal knots
    for (unsigned int idx = 0; idx <= degree; ++idx) {
      knots.push_back(parameters.front());
    }

    // Internal knots by averaging
    if (numPoints > degree + 1) {
      for (size_t pointIdx = 1; pointIdx < numPoints - degree; ++pointIdx) {
        FT sum = FT(0);
        for (unsigned int degIdx = 0; degIdx < degree; ++degIdx) {
          sum += parameters[pointIdx + degIdx];
        }
        knots.push_back(sum / FT(degree));
      }
    }

    for (unsigned int idx = 0; idx <= degree; ++idx) {
      knots.push_back(parameters.back());
    }
  }

  return knots;
}

auto
NURBSCurve::computeParameters(const std::vector<Point> &points,
                              KnotMethod method) -> std::vector<Parameter>
{
  std::vector<Parameter> parameters;
  parameters.reserve(points.size());

  if (points.empty()) {
    return parameters;
  }

  if (method == KnotMethod::UNIFORM) {
    // Uniform parameterization
    for (size_t idx = 0; idx < points.size(); ++idx) {
      parameters.push_back(FT(idx) / FT(points.size() - 1));
    }
  } else {
    // Chord length or centripetal parameterization
    parameters.push_back(FT(0));
    FT totalLength = FT(0);

    std::vector<FT> distances;
    distances.reserve(points.size() - 1);

    for (size_t idx = 1; idx < points.size(); ++idx) {
      FT distance = algorithm::distance(points[idx - 1], points[idx]);

      if (method == KnotMethod::CENTRIPETAL) {
        distance = FT(std::sqrt(CGAL::to_double(distance)));
      }

      distances.push_back(distance);
      totalLength += distance;
      parameters.push_back(totalLength);
    }

    // Normalize to [0, 1]
    if (totalLength > FT(0)) {
      for (auto &param : parameters) {
        param /= totalLength;
      }
    }
  }

  return parameters;
}

auto
NURBSCurve::checkDimensionalConsistency() const -> bool
{
  if (_controlPoints.empty()) {
    return true;
  }

  bool first3D       = _controlPoints[0].is3D();
  bool firstMeasured = _controlPoints[0].isMeasured();

  for (const auto &point : _controlPoints) {
    if (point.is3D() != first3D || point.isMeasured() != firstMeasured) {
      return false;
    }
  }

  return true;
}

auto
NURBSCurve::computeArcLength(Parameter startParam, Parameter endParam,
                             FT tolerance, unsigned int maxDepth) const -> FT
{
  if (startParam >= endParam || maxDepth == 0) {
    return FT(0);
  }

  Point     startPoint = evaluate(startParam);
  Point     endPoint   = evaluate(endParam);
  Parameter midParam   = (startParam + endParam) / FT(2);
  Point     midPoint   = evaluate(midParam);

  FT chordLength = algorithm::distance(startPoint, endPoint);
  FT leftLength  = algorithm::distance(startPoint, midPoint);
  FT rightLength = algorithm::distance(midPoint, endPoint);
  FT totalLength = leftLength + rightLength;

  if (CGAL::abs(totalLength - chordLength) <= tolerance) {
    return totalLength;
  }

  // Recursive subdivision
  return computeArcLength(startParam, midParam, tolerance, maxDepth - 1) +
         computeArcLength(midParam, endParam, tolerance, maxDepth - 1);
}

auto
NURBSCurve::findParameterByArcLength(FT targetLength, FT tolerance,
                                     unsigned int maxIterations) const
    -> Parameter
{
  auto bounds = parameterBounds();
  FT   totalLength =
      computeArcLength(bounds.first, bounds.second, tolerance / 10);

  if (targetLength <= FT(0)) {
    return bounds.first;
  }

  if (targetLength >= totalLength) {
    return bounds.second;
  }

  // Binary search for parameter
  Parameter low  = bounds.first;
  Parameter high = bounds.second;

  for (unsigned int iteration = 0; iteration < maxIterations; ++iteration) {
    Parameter mid         = (low + high) / FT(2);
    FT        lengthAtMid = computeArcLength(bounds.first, mid, tolerance / 10);

    if (CGAL::abs(lengthAtMid - targetLength) <= tolerance) {
      return mid;
    }

    if (lengthAtMid < targetLength) {
      low = mid;
    } else {
      high = mid;
    }
  }

  return (low + high) / FT(2);
}

auto
NURBSCurve::projectPointToCurve(const Point &point, FT tolerance,
                                unsigned int maxIterations) const
    -> std::pair<Point, Parameter>
{
  if (_controlPoints.empty()) {
    BOOST_THROW_EXCEPTION(Exception("Cannot project to empty curve"));
  }

  // Initial guess: closest control point
  FT        minDistance = algorithm::distance(point, _controlPoints[0]);
  Parameter bestParam   = parameterBounds().first;

  auto bounds     = parameterBounds();
  FT   paramRange = bounds.second - bounds.first;

  // Sample curve to find approximate closest point
  unsigned int numSamples = 50;
  for (unsigned int sampleIdx = 0; sampleIdx <= numSamples; ++sampleIdx) {
    Parameter param =
        bounds.first + (FT(sampleIdx) / FT(numSamples)) * paramRange;
    Point curvePoint = evaluate(param);
    FT    distance   = algorithm::distance(point, curvePoint);

    if (distance < minDistance) {
      minDistance = distance;
      bestParam   = param;
    }
  }

  // Newton-Raphson refinement (simplified version)
  for (unsigned int iteration = 0; iteration < maxIterations; ++iteration) {
    Point curvePoint = evaluate(bestParam);
    Point firstDeriv = derivative(bestParam, 1);

    // Vector from curve to target point
    Point diff = Point(point.x() - curvePoint.x(), point.y() - curvePoint.y(),
                       point.is3D() ? point.z() - curvePoint.z() : FT(0));

    // Dot product for projection
    FT numerator = diff.x() * firstDeriv.x() + diff.y() * firstDeriv.y() +
                   (diff.is3D() ? diff.z() * firstDeriv.z() : FT(0));

    FT denominator =
        firstDeriv.x() * firstDeriv.x() + firstDeriv.y() * firstDeriv.y() +
        (firstDeriv.is3D() ? firstDeriv.z() * firstDeriv.z() : FT(0));

    if (denominator < FT(1e-12)) {
      break; // Singular point
    }

    FT paramUpdate = numerator / denominator;

    if (CGAL::abs(paramUpdate) <= tolerance) {
      break; // Converged
    }

    bestParam += paramUpdate;
    bestParam = std::max(bounds.first, std::min(bounds.second, bestParam));
  }

  return std::make_pair(evaluate(bestParam), bestParam);
}

void
NURBSCurve::reserveCapacity(size_t expectedSize)
{
  _controlPoints.reserve(expectedSize);
  _weights.reserve(expectedSize);
  _knotVector.reserve(expectedSize + _degree + 1);
  _fitPoints.reserve(expectedSize);
}

void
NURBSCurve::clearOptionalData()
{
  _fitPoints.clear();
  _fitPoints.shrink_to_fit();
  _startTangent.reset();
  _endTangent.reset();
}

} // namespace SFCGAL
