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

  double startAngleDouble = CGAL::to_double(startAngle);
  double endAngleDouble   = CGAL::to_double(endAngle);
  double radiusDouble     = CGAL::to_double(radius);

  // Normalize angle range
  double angleSpan = endAngleDouble - startAngleDouble;
  while (angleSpan < 0.0) {
    angleSpan += 2.0 * M_PI;
  }
  while (angleSpan > 2.0 * M_PI) {
    angleSpan -= 2.0 * M_PI;
  }

  CoordinateType dimType = COORDINATE_XY;
  if (center.is3D() && center.isMeasured()) {
    dimType = COORDINATE_XYZM;
  } else if (center.is3D()) {
    dimType = COORDINATE_XYZ;
  } else if (center.isMeasured()) {
    dimType = COORDINATE_XYM;
  }

  std::vector<Point> controlPoints;
  std::vector<FT>    weights;

  // For arcs <= π, use standard quadratic rational representation
  if (angleSpan <= M_PI + 1e-10) {
    double halfSpan    = angleSpan / 2.0;
    double midAngle    = (startAngleDouble + endAngleDouble) / 2.0;
    double cosHalfSpan = std::cos(halfSpan);

    // Start point
    controlPoints.emplace_back(
        center.x() + FT(radiusDouble * std::cos(startAngleDouble)),
        center.y() + FT(radiusDouble * std::sin(startAngleDouble)),
        center.is3D() ? center.z() : FT(0),
        center.isMeasured() ? center.m() : NaN(), dimType);
    weights.push_back(FT(1.0));

    // Middle control point (off the circle)
    controlPoints.emplace_back(
        center.x() + FT(radiusDouble * std::cos(midAngle) / cosHalfSpan),
        center.y() + FT(radiusDouble * std::sin(midAngle) / cosHalfSpan),
        center.is3D() ? center.z() : FT(0),
        center.isMeasured() ? center.m() : NaN(), dimType);
    weights.push_back(FT(cosHalfSpan));

    // End point
    controlPoints.emplace_back(
        center.x() + FT(radiusDouble * std::cos(endAngleDouble)),
        center.y() + FT(radiusDouble * std::sin(endAngleDouble)),
        center.is3D() ? center.z() : FT(0),
        center.isMeasured() ? center.m() : NaN(), dimType);
    weights.push_back(FT(1.0));

    std::vector<Knot> knots = {FT(0), FT(0), FT(0), FT(1), FT(1), FT(1)};
    return std::make_unique<NURBSCurve>(controlPoints, weights, 2, knots);
  }

  // For arcs > π, split into two segments
  double midAngle = startAngleDouble + angleSpan / 2.0;

  // First segment: start to mid
  double halfSpan1    = (midAngle - startAngleDouble) / 2.0;
  double midAngle1    = (startAngleDouble + midAngle) / 2.0;
  double cosHalfSpan1 = std::cos(halfSpan1);

  controlPoints.emplace_back(
      center.x() + FT(radiusDouble * std::cos(startAngleDouble)),
      center.y() + FT(radiusDouble * std::sin(startAngleDouble)),
      center.is3D() ? center.z() : FT(0),
      center.isMeasured() ? center.m() : NaN(), dimType);
  weights.push_back(FT(1.0));

  controlPoints.emplace_back(
      center.x() + FT(radiusDouble * std::cos(midAngle1) / cosHalfSpan1),
      center.y() + FT(radiusDouble * std::sin(midAngle1) / cosHalfSpan1),
      center.is3D() ? center.z() : FT(0),
      center.isMeasured() ? center.m() : NaN(), dimType);
  weights.push_back(FT(cosHalfSpan1));

  controlPoints.emplace_back(center.x() + FT(radiusDouble * std::cos(midAngle)),
                             center.y() + FT(radiusDouble * std::sin(midAngle)),
                             center.is3D() ? center.z() : FT(0),
                             center.isMeasured() ? center.m() : NaN(), dimType);
  weights.push_back(FT(1.0));

  // Second segment: mid to end
  double halfSpan2    = (endAngleDouble - midAngle) / 2.0;
  double midAngle2    = (midAngle + endAngleDouble) / 2.0;
  double cosHalfSpan2 = std::cos(halfSpan2);

  controlPoints.emplace_back(
      center.x() + FT(radiusDouble * std::cos(midAngle2) / cosHalfSpan2),
      center.y() + FT(radiusDouble * std::sin(midAngle2) / cosHalfSpan2),
      center.is3D() ? center.z() : FT(0),
      center.isMeasured() ? center.m() : NaN(), dimType);
  weights.push_back(FT(cosHalfSpan2));

  controlPoints.emplace_back(
      center.x() + FT(radiusDouble * std::cos(endAngleDouble)),
      center.y() + FT(radiusDouble * std::sin(endAngleDouble)),
      center.is3D() ? center.z() : FT(0),
      center.isMeasured() ? center.m() : NaN(), dimType);
  weights.push_back(FT(1.0));

  // Knot vector for two Bézier segments
  std::vector<Knot> knots = {FT(0),   FT(0), FT(0), FT(0.5),
                             FT(0.5), FT(1), FT(1), FT(1)};

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

  auto parameters = computeParameters(points, knotMethod);

  std::vector<Knot> knots;
  size_t            numKnots = points.size() + degree + 1;
  knots.reserve(numKnots);

  for (unsigned int idx = 0; idx <= degree; ++idx) {
    knots.push_back(parameters.front());
  }

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

  auto curve = std::make_unique<NURBSCurve>(
      points, std::vector<FT>(points.size(), FT(1)), degree, knots);

  curve->_fitPoints    = points;
  curve->_fitTolerance = FT(0);

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

  size_t numControlPoints = std::min(maxControlPoints, points.size());
  std::vector<Point> controlPoints;
  controlPoints.reserve(numControlPoints);

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
    BOOST_THROW_EXCEPTION(
        Exception("Normal vector not uniquely defined for 3D curves"));
  }

  Point tangentVec = tangent(parameter);
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
    Point binormalVec = Point(0, 0, 1);
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

  std::vector<Parameter> parameters;
  parameters.reserve(maxSegments + 1);

  FT deltaParam = (bounds.second - bounds.first) / FT(minSegments);
  for (unsigned int segIdx = 0; segIdx <= minSegments; ++segIdx) {
    parameters.push_back(bounds.first + FT(segIdx) * deltaParam);
  }

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

      Point interpolated = Point(
          (point1.x() + point2.x()) / FT(2), (point1.y() + point2.y()) / FT(2),
          point1.is3D() ? (point1.z() + point2.z()) / FT(2) : FT(0));

      FT deviation = algorithm::distance(midPoint, interpolated);

      if (deviation > tolerance) {
        parameters.insert(parameters.begin() + paramIdx + 1, midParam);
        refined = true;
        break;
      }
    }
  }

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
      *plane = {FT(0), FT(0), FT(1), FT(0)};
    }
    return true;
  }

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

        Point normal = Point(vec1.y() * vec2.z() - vec1.z() * vec2.y(),
                             vec1.z() * vec2.x() - vec1.x() * vec2.z(),
                             vec1.x() * vec2.y() - vec1.y() * vec2.x());

        FT normalMagnitude = std::sqrt(
            CGAL::to_double(normal.x() * normal.x() + normal.y() * normal.y() +
                            normal.z() * normal.z()));

        if (normalMagnitude > FT(1e-10)) {
          normal =
              Point(normal.x() / normalMagnitude, normal.y() / normalMagnitude,
                    normal.z() / normalMagnitude);

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
    return false;
  }

  Point direction = Point(_controlPoints[1].x() - _controlPoints[0].x(),
                          _controlPoints[1].y() - _controlPoints[0].y(),
                          _controlPoints[1].is3D()
                              ? _controlPoints[1].z() - _controlPoints[0].z()
                              : FT(0));

  FT dirMagnitude = std::sqrt(CGAL::to_double(
      direction.x() * direction.x() + direction.y() * direction.y() +
      (direction.is3D() ? direction.z() * direction.z() : FT(0))));

  if (dirMagnitude < FT(1e-10)) {
    return true;
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
  if (isEmpty()) {
    BOOST_THROW_EXCEPTION(Exception("Cannot compute length of empty curve"));
  }

  auto bounds = parameterBounds();

  if (from < FT(0))
    from = bounds.first;
  if (to < FT(0))
    to = bounds.second;

  if (from >= to) {
    return FT(0);
  }

  tolerance = CGAL::abs(tolerance);
  if (tolerance <= FT(0)) {
    tolerance = FT(1e-6);
  }

  return computeArcLength(from, to, tolerance);
}

auto
NURBSCurve::parameterAtLength(FT arcLength, FT tolerance) const -> Parameter
{
  if (arcLength <= FT(0)) {
    return parameterBounds().first;
  }

  FT totalLength = length();
  if (arcLength >= totalLength) {
    return parameterBounds().second;
  }

  return findParameterByArcLength(arcLength, tolerance);
}

auto
NURBSCurve::reparameterizeByArcLength() const -> std::unique_ptr<Curve>
{
  if (isEmpty()) {
    return std::make_unique<NURBSCurve>();
  }

  FT totalLength = length();

  if (totalLength <= FT(1e-12)) {
    return std::unique_ptr<NURBSCurve>(clone());
  }

  // Use a reasonable number of samples for reparameterization
  size_t numSamples =
      std::min(static_cast<size_t>(32),
               std::max(static_cast<size_t>(8), _controlPoints.size()));

  std::vector<Point> newControlPoints;
  newControlPoints.reserve(numSamples);

  for (size_t idx = 0; idx < numSamples; ++idx) {
    FT        targetLength = (FT(idx) / FT(numSamples - 1)) * totalLength;
    Parameter param        = parameterAtLength(targetLength, FT(1e-8));
    newControlPoints.push_back(evaluate(param));
  }

  // Use a reasonable degree for the reparameterized curve
  unsigned int newDegree =
      std::min(3u, static_cast<unsigned int>(newControlPoints.size() - 1));

  return std::make_unique<NURBSCurve>(newControlPoints, newDegree,
                                      KnotMethod::UNIFORM);
}

auto
NURBSCurve::split(Parameter parameter) const
    -> std::pair<std::unique_ptr<Curve>, std::unique_ptr<Curve>>
{
  if (isEmpty()) {
    return std::make_pair(std::make_unique<NURBSCurve>(),
                          std::make_unique<NURBSCurve>());
  }

  auto bounds = parameterBounds();

  parameter = std::max(bounds.first, std::min(bounds.second, parameter));

  if (parameter <= bounds.first) {
    auto emptyCurve = std::make_unique<NURBSCurve>();
    auto fullCurve  = std::unique_ptr<NURBSCurve>(clone());
    return std::make_pair(std::move(emptyCurve), std::move(fullCurve));
  }

  if (parameter >= bounds.second) {
    auto fullCurve  = std::unique_ptr<NURBSCurve>(clone());
    auto emptyCurve = std::make_unique<NURBSCurve>();
    return std::make_pair(std::move(fullCurve), std::move(emptyCurve));
  }

  auto firstCurve  = subcurve(bounds.first, parameter);
  auto secondCurve = subcurve(parameter, bounds.second);

  auto firstNurbs = std::unique_ptr<NURBSCurve>(
      dynamic_cast<NURBSCurve *>(firstCurve.release()));
  auto secondNurbs = std::unique_ptr<NURBSCurve>(
      dynamic_cast<NURBSCurve *>(secondCurve.release()));

  if (!firstNurbs) {
    firstNurbs = std::make_unique<NURBSCurve>();
  }
  if (!secondNurbs) {
    secondNurbs = std::make_unique<NURBSCurve>();
  }

  return std::make_pair(std::move(firstNurbs), std::move(secondNurbs));
}

auto
NURBSCurve::subcurve(Parameter from, Parameter to) const
    -> std::unique_ptr<Curve>
{
  if (isEmpty()) {
    return std::make_unique<NURBSCurve>();
  }

  if (from >= to) {
    return std::make_unique<NURBSCurve>();
  }

  auto bounds = parameterBounds();
  from        = std::max(from, bounds.first);
  to          = std::min(to, bounds.second);

  if (from >= to) {
    return std::make_unique<NURBSCurve>();
  }

  unsigned int numSamples =
      std::max(16u, static_cast<unsigned int>(_controlPoints.size()));

  std::vector<Point> subPoints;
  subPoints.reserve(numSamples + 1);

  FT paramRange = to - from;

  for (unsigned int sampleIdx = 0; sampleIdx <= numSamples; ++sampleIdx) {
    Parameter param = from + (FT(sampleIdx) / FT(numSamples)) * paramRange;
    Point     evaluatedPoint = evaluate(param);
    subPoints.push_back(evaluatedPoint);
  }

  unsigned int subDegree =
      std::min(_degree, static_cast<unsigned int>(subPoints.size() - 1));

  auto subKnots =
      generateKnotVector(subPoints, subDegree, KnotMethod::CHORD_LENGTH);

  std::vector<FT> subWeights;
  if (isRational()) {
    subWeights.resize(subPoints.size(), FT(1));
  }

  if (subWeights.empty()) {
    return std::make_unique<NURBSCurve>(subPoints, subDegree,
                                        KnotMethod::CHORD_LENGTH);
  } else {
    return std::make_unique<NURBSCurve>(subPoints, subWeights, subDegree,
                                        subKnots);
  }
}

auto
NURBSCurve::reverse() const -> std::unique_ptr<Curve>
{
  if (_controlPoints.empty()) {
    return std::make_unique<NURBSCurve>();
  }

  std::vector<Point> reversedPoints(_controlPoints.rbegin(),
                                    _controlPoints.rend());
  std::vector<FT>    reversedWeights;

  if (!_weights.empty()) {
    reversedWeights.assign(_weights.rbegin(), _weights.rend());
  }

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

  Point thisEnd    = endPoint();
  Point otherStart = otherNurbs->startPoint();

  if (algorithm::distance(thisEnd, otherStart) > tolerance) {
    BOOST_THROW_EXCEPTION(
        Exception("Curves are not adjacent within tolerance"));
  }

  std::vector<Point> joinedPoints = _controlPoints;
  const auto        &otherPoints  = otherNurbs->controlPoints();

  joinedPoints.insert(joinedPoints.end(), otherPoints.begin() + 1,
                      otherPoints.end());

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
  auto lineString = toLineString(100);

  if (intersections) {
    intersections->clear();
  }

  for (size_t firstIdx = 0; firstIdx < lineString->numPoints() - 3;
       ++firstIdx) {
    for (size_t secondIdx = firstIdx + 2;
         secondIdx < lineString->numPoints() - 1; ++secondIdx) {
      if (secondIdx <= firstIdx + 2)
        continue;

      Point segment1Start = lineString->pointN(firstIdx);
      Point segment1End   = lineString->pointN(firstIdx + 1);
      Point segment2Start = lineString->pointN(secondIdx);
      Point segment2End   = lineString->pointN(secondIdx + 1);

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

  FT minX = _controlPoints[0].x(), maxX = _controlPoints[0].x();
  FT minY = _controlPoints[0].y(), maxY = _controlPoints[0].y();
  FT minZ = _controlPoints[0].is3D() ? _controlPoints[0].z() : FT(0);
  FT maxZ = minZ;

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

  if (!_weights.empty() && _weights.size() != controlPoints.size()) {
    _weights.resize(controlPoints.size(), FT(1));
  }

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
    return FT(1);
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

  for (unsigned int idx = 0; idx <= _degree; ++idx) {
    if (CGAL::abs(_knotVector[idx] - _knotVector[0]) > tolerance) {
      return false;
    }
  }

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

//-- Advanced NURBS operations

auto
NURBSCurve::insertKnot(Knot parameter, unsigned int times) const
    -> std::unique_ptr<NURBSCurve>
{
  return std::unique_ptr<NURBSCurve>(clone());
}

auto
NURBSCurve::refineKnotVector(const std::vector<Knot> &newKnots) const
    -> std::unique_ptr<NURBSCurve>
{
  return std::unique_ptr<NURBSCurve>(clone());
}

auto
NURBSCurve::elevateDegree(unsigned int times) const
    -> std::unique_ptr<NURBSCurve>
{
  return std::unique_ptr<NURBSCurve>(clone());
}

auto
NURBSCurve::reduceDegree(FT tolerance) const -> std::unique_ptr<NURBSCurve>
{
  return std::unique_ptr<NURBSCurve>(clone());
}

auto
NURBSCurve::removeKnots(FT tolerance) const -> std::unique_ptr<NURBSCurve>
{
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
  if (controlPoints.empty()) {
    if (!knots.empty() || !weights.empty()) {
      return {false, "Empty curve with non-empty knot vector or weights"};
    }
    return {true, ""};
  }

  for (const auto &point : controlPoints) {
    if ((controlPoints[0].is3D() != point.is3D() ||
         controlPoints[0].isMeasured() != point.isMeasured())) {
      return {false, "Control point dimension mismatch"};
    }
  }

  if (degree >= controlPoints.size()) {
    return {false, "Degree must be less than control points count"};
  }

  const size_t expectedKnots = controlPoints.size() + degree + 1;
  if (knots.size() != expectedKnots) {
    return {false, "Invalid knot vector size"};
  }

  for (size_t idx = 1; idx < knots.size(); ++idx) {
    if (knots[idx] < knots[idx - 1]) {
      return {false, "Knot vector is not non-decreasing"};
    }
  }

  if (!weights.empty()) {
    if (weights.size() != controlPoints.size()) {
      return {false, "Weight count does not match control point count"};
    }

    for (const auto &weight : weights) {
      if (weight <= FT(0) || !CGAL::is_finite(weight)) {
        return {false, "Invalid weight value"};
      }
    }
  }

  if (!knots.empty() && controlPoints.size() > 0) {
    FT startParam = knots[degree];
    FT endParam   = knots[controlPoints.size()];

    if (startParam > endParam) {
      return {false, "Invalid parameter bounds"};
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

  if (parameter >= _knotVector[numControlPoints]) {
    return numControlPoints - 1;
  }

  if (parameter <= _knotVector[_degree]) {
    return _degree;
  }

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
  struct HomogeneousPoint {
    FT     weightedX, weightedY, weightedZ, weight;
    double measure;

    HomogeneousPoint()
        : weightedX(0), weightedY(0), weightedZ(0), weight(1), measure(0)
    {
    }

    HomogeneousPoint(const Point &point, FT pointWeight)
        : weightedX(point.x() * pointWeight),
          weightedY(point.y() * pointWeight),
          weightedZ(point.is3D() ? point.z() * pointWeight : FT(0)),
          weight(pointWeight), measure(point.isMeasured() ? point.m() : 0.0)
    {
    }
  };

  std::vector<HomogeneousPoint> temp(_degree + 1);

  for (unsigned int controlIdx = 0; controlIdx <= _degree; ++controlIdx) {
    size_t cpIndex = span - _degree + controlIdx;
    if (cpIndex >= _controlPoints.size()) {
      cpIndex = _controlPoints.size() - 1;
    }
    FT pointWeight   = _weights.empty() ? FT(1) : _weights[cpIndex];
    temp[controlIdx] = HomogeneousPoint(_controlPoints[cpIndex], pointWeight);
  }

  for (unsigned int level = 1; level <= _degree; ++level) {
    for (unsigned int basisIdx = _degree; basisIdx >= level; --basisIdx) {
      size_t knotIdx = span - _degree + basisIdx;

      if (knotIdx + _degree - level + 1 >= _knotVector.size() ||
          knotIdx >= _knotVector.size()) {
        continue;
      }

      FT denominator =
          _knotVector[knotIdx + _degree - level + 1] - _knotVector[knotIdx];
      FT alpha = FT(0);

      if (CGAL::abs(denominator) > FT(1e-15)) {
        alpha = (parameter - _knotVector[knotIdx]) / denominator;
      }

      alpha = std::max(FT(0), std::min(FT(1), alpha));

      temp[basisIdx].weightedX =
          (FT(1) - alpha) * temp[basisIdx - 1].weightedX +
          alpha * temp[basisIdx].weightedX;
      temp[basisIdx].weightedY =
          (FT(1) - alpha) * temp[basisIdx - 1].weightedY +
          alpha * temp[basisIdx].weightedY;
      temp[basisIdx].weightedZ =
          (FT(1) - alpha) * temp[basisIdx - 1].weightedZ +
          alpha * temp[basisIdx].weightedZ;
      temp[basisIdx].weight = (FT(1) - alpha) * temp[basisIdx - 1].weight +
                              alpha * temp[basisIdx].weight;

      if (isMeasured()) {
        double alphaDouble = CGAL::to_double(alpha);
        temp[basisIdx].measure =
            (1.0 - alphaDouble) * temp[basisIdx - 1].measure +
            alphaDouble * temp[basisIdx].measure;
      }
    }
  }

  HomogeneousPoint &result = temp[_degree];

  if (CGAL::abs(result.weight) < FT(1e-15)) {
    BOOST_THROW_EXCEPTION(Exception("Division by zero in NURBS evaluation"));
  }

  FT invWeight = FT(1) / result.weight;

  if (is3D() && isMeasured()) {
    return Point(result.weightedX * invWeight, result.weightedY * invWeight,
                 result.weightedZ * invWeight, result.measure);
  } else if (is3D()) {
    return Point(result.weightedX * invWeight, result.weightedY * invWeight,
                 result.weightedZ * invWeight);
  } else if (isMeasured()) {
    Point point(result.weightedX * invWeight, result.weightedY * invWeight);
    point.setM(result.measure);
    return point;
  } else {
    return Point(result.weightedX * invWeight, result.weightedY * invWeight);
  }
}

auto
NURBSCurve::computeDerivatives(Parameter parameter, unsigned int maxOrder) const
    -> std::vector<Point>
{
  std::vector<Point> derivatives;
  derivatives.reserve(maxOrder + 1);

  derivatives.push_back(evaluate(parameter));

  if (maxOrder > 0) {
    if (_degree == 1 && _controlPoints.size() == 2) {
      Point derivative =
          Point(_controlPoints[1].x() - _controlPoints[0].x(),
                _controlPoints[1].y() - _controlPoints[0].y(),
                _controlPoints[1].is3D()
                    ? _controlPoints[1].z() - _controlPoints[0].z()
                    : FT(0));
      derivatives.push_back(derivative);

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

  if (degree == 0) {
    knots.reserve(numKnots);
    for (size_t pointIdx = 0; pointIdx < numPoints; ++pointIdx) {
      knots.push_back(FT(pointIdx));
    }
    knots.push_back(FT(numPoints));
    return knots;
  }

  knots.reserve(numKnots);

  if (method == KnotMethod::UNIFORM) {
    for (unsigned int degreeIdx = 0; degreeIdx <= degree; ++degreeIdx) {
      knots.push_back(FT(0));
    }

    size_t numInternal = numKnots - 2 * (degree + 1);
    for (size_t internalIdx = 1; internalIdx <= numInternal; ++internalIdx) {
      knots.push_back(FT(internalIdx) / FT(numInternal + 1));
    }

    for (unsigned int degreeIdx = 0; degreeIdx <= degree; ++degreeIdx) {
      knots.push_back(FT(1));
    }
  } else {
    auto parameters = computeParameters(points, method);

    for (unsigned int degreeIdx = 0; degreeIdx <= degree; ++degreeIdx) {
      knots.push_back(parameters.front());
    }

    if (numPoints > degree + 1) {
      for (size_t pointIdx = 1; pointIdx < numPoints - degree; ++pointIdx) {
        FT sum = FT(0);
        for (unsigned int degreeIdx = 0; degreeIdx < degree; ++degreeIdx) {
          sum += parameters[pointIdx + degreeIdx];
        }
        knots.push_back(sum / FT(degree));
      }
    }

    for (unsigned int degreeIdx = 0; degreeIdx <= degree; ++degreeIdx) {
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
    for (size_t pointIdx = 0; pointIdx < points.size(); ++pointIdx) {
      parameters.push_back(FT(pointIdx) / FT(points.size() - 1));
    }
  } else {
    parameters.push_back(FT(0));
    FT totalLength = FT(0);

    std::vector<FT> distances;
    distances.reserve(points.size() - 1);

    for (size_t pointIdx = 1; pointIdx < points.size(); ++pointIdx) {
      FT distance = algorithm::distance(points[pointIdx - 1], points[pointIdx]);

      if (method == KnotMethod::CENTRIPETAL) {
        distance = FT(std::sqrt(CGAL::to_double(distance)));
      }

      distances.push_back(distance);
      totalLength += distance;
      parameters.push_back(totalLength);
    }

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
                             FT tolerance) const -> FT
{
  if (startParam >= endParam) {
    return FT(0);
  }

  // Use simple trapezoidal rule with fixed subdivision
  const unsigned int numSegments = 32; // Fixed number of segments for stability

  FT paramStep   = (endParam - startParam) / FT(numSegments);
  FT totalLength = FT(0);

  Point prevPoint = evaluate(startParam);

  for (unsigned int segIdx = 1; segIdx <= numSegments; ++segIdx) {
    Parameter currentParam = startParam + FT(segIdx) * paramStep;
    Point     currentPoint = evaluate(currentParam);

    FT segmentLength = algorithm::distance(prevPoint, currentPoint);
    totalLength += segmentLength;

    prevPoint = currentPoint;
  }

  return totalLength;
}

auto
NURBSCurve::findParameterByArcLength(FT targetLength, FT tolerance) const
    -> Parameter
{
  auto bounds = parameterBounds();

  // Use binary search with fixed iterations to avoid infinite loops
  Parameter low  = bounds.first;
  Parameter high = bounds.second;

  const unsigned int maxIterations =
      50; // Limit iterations to prevent infinite loops

  for (unsigned int iter = 0; iter < maxIterations; ++iter) {
    Parameter mid         = (low + high) / FT(2);
    FT        lengthAtMid = computeArcLength(bounds.first, mid, tolerance);

    if (CGAL::abs(lengthAtMid - targetLength) <= tolerance) {
      return mid;
    }

    if (lengthAtMid < targetLength) {
      low = mid;
    } else {
      high = mid;
    }

    // Check for convergence
    if (high - low < FT(1e-12)) {
      break;
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

  FT        minDistance = algorithm::distance(point, _controlPoints[0]);
  Parameter bestParam   = parameterBounds().first;

  auto bounds     = parameterBounds();
  FT   paramRange = bounds.second - bounds.first;

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

  for (unsigned int iteration = 0; iteration < maxIterations; ++iteration) {
    Point curvePoint = evaluate(bestParam);
    Point firstDeriv = derivative(bestParam, 1);

    Point diff = Point(point.x() - curvePoint.x(), point.y() - curvePoint.y(),
                       point.is3D() ? point.z() - curvePoint.z() : FT(0));

    FT numerator = diff.x() * firstDeriv.x() + diff.y() * firstDeriv.y() +
                   (diff.is3D() ? diff.z() * firstDeriv.z() : FT(0));

    FT denominator =
        firstDeriv.x() * firstDeriv.x() + firstDeriv.y() * firstDeriv.y() +
        (firstDeriv.is3D() ? firstDeriv.z() * firstDeriv.z() : FT(0));

    if (denominator < FT(1e-12)) {
      break;
    }

    FT paramUpdate = numerator / denominator;

    if (CGAL::abs(paramUpdate) <= tolerance) {
      break;
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
