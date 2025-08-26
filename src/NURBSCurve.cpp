#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/Exception.h"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace SFCGAL {

NURBSCurve::NURBSCurve() = default;

NURBSCurve::NURBSCurve(const std::vector<Point> &controlPoints,
                       unsigned int              degree)
    : BSplineCurve(controlPoints, degree)
{
  // Initialize uniform weights
  _weights.resize(controlPoints.size(), 1.0);
}

NURBSCurve::NURBSCurve(const std::vector<Point> &controlPoints,
                       unsigned int degree, const std::vector<double> &weights)
    : BSplineCurve(controlPoints, degree), _weights(weights)
{
  if (!validateWeights()) {
    throw InappropriateGeometryException("Invalid weights for NURBS curve");
  }
}

NURBSCurve::NURBSCurve(const std::vector<Point> &controlPoints,
                       unsigned int degree, const std::vector<double> &weights,
                       const std::vector<double> &knotVector)
    : BSplineCurve(controlPoints, degree, knotVector), _weights(weights)
{
  if (!validateWeights()) {
    throw InappropriateGeometryException("Invalid weights for NURBS curve");
  }
}

NURBSCurve::NURBSCurve(const NURBSCurve &other) = default;

auto
NURBSCurve::operator=(const NURBSCurve &other) -> NURBSCurve &
{
  if (this != &other) {
    BSplineCurve::operator=(other);
    _weights = other._weights;
  }
  return *this;
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
NURBSCurve::evaluate(double parameter) const -> Point
{
  if (_controlPoints.empty()) {
    throw InappropriateGeometryException("Cannot evaluate empty NURBS curve");
  }

  auto bounds = parameterBounds();
  if (parameter < bounds.first || parameter > bounds.second) {
    throw InappropriateGeometryException("Parameter t is outside valid range");
  }

  // If all weights are uniform, fall back to B-spline evaluation
  if (!isRational()) {
    return BSplineCurve::evaluate(parameter);
  }

  return evaluateRational(parameter);
}

auto
NURBSCurve::evaluateRational(double parameter) const -> Point
{
  // De Boor's algorithm with weights (rational version)
  size_t span = findSpan(parameter);

  // Prepare homogeneous control points
  struct HomogeneousPoint {
    double weightedX;
    double weightedY;
    double weightedZ;
    double weight;

    HomogeneousPoint() : weightedX(0), weightedY(0), weightedZ(0), weight(0) {}

    HomogeneousPoint(const Point &point, double pointWeight, bool is3D)
    {
      weightedX = CGAL::to_double(point.x()) * pointWeight;
      weightedY = CGAL::to_double(point.y()) * pointWeight;
      weightedZ = is3D ? CGAL::to_double(point.z()) * pointWeight : 0.0;
      weight    = pointWeight;
    }

    [[nodiscard]] auto
    lerp(const HomogeneousPoint &other, double alpha) const -> HomogeneousPoint
    {
      HomogeneousPoint result;
      result.weightedX = ((1 - alpha) * weightedX) + (alpha * other.weightedX);
      result.weightedY = ((1 - alpha) * weightedY) + (alpha * other.weightedY);
      result.weightedZ = ((1 - alpha) * weightedZ) + (alpha * other.weightedZ);
      result.weight    = ((1 - alpha) * weight) + (alpha * other.weight);
      return result;
    }
  };

  std::vector<HomogeneousPoint> temp(_degree + 1);

  // Initialize with weighted control points
  for (size_t pointIndex = 0; pointIndex <= _degree; ++pointIndex) {
    size_t idx = span - _degree + pointIndex;
    temp[pointIndex] =
        HomogeneousPoint(_controlPoints[idx], _weights[idx], is3D());
  }

  // Apply de Boor's algorithm in homogeneous space
  for (size_t level = 1; level <= _degree; ++level) {
    for (size_t pointIndex = _degree; pointIndex >= level; --pointIndex) {
      size_t knotIdx = span - _degree + pointIndex;
      double alpha =
          (parameter - _knotVector[knotIdx]) /
          (_knotVector[knotIdx + _degree - level + 1] - _knotVector[knotIdx]);
      temp[pointIndex] = temp[pointIndex - 1].lerp(temp[pointIndex], alpha);
    }
  }

  // Project back to Euclidean space
  if (std::abs(temp[_degree].weight) < 1e-10) {
    throw InappropriateGeometryException(
        "Division by zero in NURBS evaluation");
  }

  double coordX = temp[_degree].weightedX / temp[_degree].weight;
  double coordY = temp[_degree].weightedY / temp[_degree].weight;

  if (is3D()) {
    double coordZ = temp[_degree].weightedZ / temp[_degree].weight;
    if (isMeasured()) {
      // Interpolate M values linearly
      double segmentT = (parameter - _knotVector[span]) /
                        (_knotVector[span + 1] - _knotVector[span]);
      size_t idx1    = span - _degree;
      size_t idx2    = std::min(idx1 + 1, _controlPoints.size() - 1);
      double measure = ((1 - segmentT) * _controlPoints[idx1].m()) +
                       (segmentT * _controlPoints[idx2].m());
      return {coordX, coordY, coordZ, measure};
    }
    return {coordX, coordY, coordZ};
  }
  if (isMeasured()) {
    double segmentT = (parameter - _knotVector[span]) /
                      (_knotVector[span + 1] - _knotVector[span]);
    size_t idx1    = span - _degree;
    size_t idx2    = std::min(idx1 + 1, _controlPoints.size() - 1);
    double measure = ((1 - segmentT) * _controlPoints[idx1].m()) +
                     (segmentT * _controlPoints[idx2].m());
    auto resultPoint = Point(coordX, coordY);
    resultPoint.setM(measure);
    return resultPoint;
  }

  return {coordX, coordY};
}

auto
NURBSCurve::derivative(double parameter, unsigned int order) const -> Point
{
  if (order == 0) {
    return evaluate(parameter);
  }

  if (!isRational()) {
    return BSplineCurve::derivative(parameter, order);
  }

  // For rational NURBS, use the quotient rule for derivatives
  if (order == 1) {
    return computeRationalDerivative(parameter);
  }

  // Higher order derivatives using finite differences
  const double stepSize = 1e-8;
  auto         bounds   = parameterBounds();

  Point derivative1;
  Point derivative2;

  if (parameter - stepSize >= bounds.first) {
    derivative1 = derivative(parameter - stepSize, order - 1);
  } else {
    derivative1 = derivative(parameter, order - 1);
  }

  if (parameter + stepSize <= bounds.second) {
    derivative2 = derivative(parameter + stepSize, order - 1);
  } else {
    derivative2 = derivative(parameter, order - 1);
  }

  double deltaX =
      CGAL::to_double(derivative2.x() - derivative1.x()) / (2 * stepSize);
  double deltaY =
      CGAL::to_double(derivative2.y() - derivative1.y()) / (2 * stepSize);

  if (is3D()) {
    double deltaZ =
        CGAL::to_double(derivative2.z() - derivative1.z()) / (2 * stepSize);
    return {deltaX, deltaY, deltaZ};
  }

  return {deltaX, deltaY};
}

auto
NURBSCurve::computeRationalDerivative(double parameter) const -> Point
{
  // Compute both the curve and weight function derivatives
  size_t span = findSpan(parameter);

  // Compute basis functions and their derivatives
  std::vector<double> basisFunctions(_degree + 1);
  std::vector<double> derivativeBasisFunctions(_degree + 1);

  computeBasisFunctions(span, parameter, basisFunctions);
  computeBasisDerivatives(span, parameter, 1, derivativeBasisFunctions);

  // Compute weighted position and its derivative
  double weightedX           = 0;
  double weightedY           = 0;
  double weightedZ           = 0;
  double weight              = 0;
  double derivativeWeightedX = 0;
  double derivativeWeightedY = 0;
  double derivativeWeightedZ = 0;
  double derivativeWeight    = 0;

  for (size_t pointIndex = 0; pointIndex <= _degree; ++pointIndex) {
    size_t idx         = span - _degree + pointIndex;
    double pointWeight = _weights[idx];

    weightedX += (basisFunctions[pointIndex] *
                  CGAL::to_double(_controlPoints[idx].x())) *
                 pointWeight;
    weightedY += (basisFunctions[pointIndex] *
                  CGAL::to_double(_controlPoints[idx].y())) *
                 pointWeight;
    weight += basisFunctions[pointIndex] * pointWeight;

    derivativeWeightedX += (derivativeBasisFunctions[pointIndex] *
                            CGAL::to_double(_controlPoints[idx].x())) *
                           pointWeight;
    derivativeWeightedY += (derivativeBasisFunctions[pointIndex] *
                            CGAL::to_double(_controlPoints[idx].y())) *
                           pointWeight;
    derivativeWeight += derivativeBasisFunctions[pointIndex] * pointWeight;

    if (is3D()) {
      weightedZ += (basisFunctions[pointIndex] *
                    CGAL::to_double(_controlPoints[idx].z())) *
                   pointWeight;
      derivativeWeightedZ += (derivativeBasisFunctions[pointIndex] *
                              CGAL::to_double(_controlPoints[idx].z())) *
                             pointWeight;
    }
  }

  if (std::abs(weight) < 1e-10) {
    throw InappropriateGeometryException(
        "Division by zero in NURBS derivative");
  }

  // Apply quotient rule
  double weightSquared = weight * weight;
  double deltaX =
      ((derivativeWeightedX * weight) - (weightedX * derivativeWeight)) /
      weightSquared;
  double deltaY =
      ((derivativeWeightedY * weight) - (weightedY * derivativeWeight)) /
      weightSquared;

  if (is3D()) {
    double deltaZ =
        ((derivativeWeightedZ * weight) - (weightedZ * derivativeWeight)) /
        weightSquared;
    return {deltaX, deltaY, deltaZ};
  }

  return {deltaX, deltaY};
}

void
NURBSCurve::computeBasisFunctions(size_t span, double parameter,
                                  std::vector<double> &basisFunctions) const
{
  basisFunctions.resize(_degree + 1);
  std::vector<double> left(_degree + 1);
  std::vector<double> right(_degree + 1);

  basisFunctions[0] = 1.0;

  for (size_t degree = 1; degree <= _degree; ++degree) {
    left[degree]  = parameter - _knotVector[span + 1 - degree];
    right[degree] = _knotVector[span + degree] - parameter;

    double saved = 0.0;

    for (size_t level = 0; level < degree; ++level) {
      double temp =
          basisFunctions[level] / (right[level + 1] + left[degree - level]);
      basisFunctions[level] = saved + (right[level + 1] * temp);
      saved                 = left[degree - level] * temp;
    }

    basisFunctions[degree] = saved;
  }
}

void
NURBSCurve::computeBasisDerivatives(
    size_t span, double parameter, unsigned int order,
    std::vector<double> &derivativeBasisFunctions) const
{
  // Simplified version for first derivative only
  if (order != 1) {
    throw InappropriateGeometryException(
        "Only first derivative is implemented");
  }

  derivativeBasisFunctions.resize(_degree + 1, 0.0);

  if (_degree == 0) {
    return;
  }

  // Compute N_{i,p-1}
  std::vector<double> basisLeft(_degree);

  // Recursively compute lower degree basis functions
  std::vector<double> left(_degree);
  std::vector<double> right(_degree);

  basisLeft[0] = 1.0;

  for (size_t degree = 1; degree < _degree; ++degree) {
    left[degree]  = parameter - _knotVector[span + 1 - degree];
    right[degree] = _knotVector[span + degree] - parameter;

    double saved = 0.0;

    for (size_t level = 0; level < degree; ++level) {
      double temp =
          basisLeft[level] / (right[level + 1] + left[degree - level]);
      basisLeft[level] = saved + (right[level + 1] * temp);
      saved            = left[degree - level] * temp;
    }

    basisLeft[degree] = saved;
  }

  // Compute derivatives
  for (size_t pointIndex = 0; pointIndex <= _degree; ++pointIndex) {
    if (pointIndex > 0) {
      double denomLeft = _knotVector[span + pointIndex] -
                         _knotVector[span + pointIndex - _degree];
      if (std::abs(denomLeft) > 1e-10) {
        derivativeBasisFunctions[pointIndex] +=
            (static_cast<double>(_degree) * basisLeft[pointIndex - 1]) /
            denomLeft;
      }
    }

    if (pointIndex < _degree) {
      double denomRight = _knotVector[span + pointIndex + 1] -
                          _knotVector[span + pointIndex - _degree + 1];
      if (std::abs(denomRight) > 1e-10) {
        derivativeBasisFunctions[pointIndex] -=
            (static_cast<double>(_degree) * basisLeft[pointIndex]) / denomRight;
      }
    }
  }
}

auto
NURBSCurve::isValid() const -> bool
{
  if (!BSplineCurve::isValid()) {
    return false;
  }
  return validateWeights();
}

auto
NURBSCurve::weights() const -> const std::vector<double> &
{
  return _weights;
}

void
NURBSCurve::setWeights(const std::vector<double> &weights)
{
  _weights = weights;
  if (!validateWeights()) {
    throw InappropriateGeometryException("Invalid weights for NURBS curve");
  }
}

auto
NURBSCurve::weight(size_t index) const -> double
{
  BOOST_ASSERT(index < _weights.size());
  return _weights[index];
}

void
NURBSCurve::setWeight(size_t index, double weight)
{
  BOOST_ASSERT(index < _weights.size());
  if (weight <= 0.0) {
    throw InappropriateGeometryException("NURBS weights must be positive");
  }
  _weights[index] = weight;
}

auto
NURBSCurve::isRational() const -> bool
{
  if (_weights.empty()) {
    return false;
  }

  // Check if all weights are equal (within tolerance)
  double       firstWeight = _weights[0];
  const double tolerance   = 1e-10;

  for (size_t index = 1; index < _weights.size(); ++index) {
    if (std::abs(_weights[index] - firstWeight) > tolerance) {
      return true;
    }
  }

  return false;
}

auto
NURBSCurve::homogeneousControlPoints() const -> std::vector<std::vector<double>>
{
  std::vector<std::vector<double>> result;

  for (size_t index = 0; index < _controlPoints.size(); ++index) {
    std::vector<double> point;
    double pointWeight = (index < _weights.size()) ? _weights[index] : 1.0;

    point.push_back(CGAL::to_double(_controlPoints[index].x()) * pointWeight);
    point.push_back(CGAL::to_double(_controlPoints[index].y()) * pointWeight);

    if (is3D()) {
      point.push_back(CGAL::to_double(_controlPoints[index].z()) * pointWeight);
    }

    point.push_back(pointWeight);
    result.push_back(point);
  }

  return result;
}

void
NURBSCurve::insertKnot(double parameter)
{
  auto bounds = parameterBounds();
  if (parameter < bounds.first || parameter > bounds.second) {
    throw InappropriateGeometryException("Knot value outside valid range");
  }

  size_t span = findSpan(parameter);

  // New control points and weights
  std::vector<Point>  newControlPoints;
  std::vector<double> newWeights;

  // Copy control points before insertion point
  for (size_t index = 0; index <= span - _degree; ++index) {
    newControlPoints.push_back(_controlPoints[index]);
    newWeights.push_back(_weights[index]);
  }

  // Compute new control points affected by knot insertion
  std::vector<Point>  temp(_degree);
  std::vector<double> tempWeights(_degree);

  for (size_t index = 0; index < _degree; ++index) {
    size_t idx   = span - _degree + index + 1;
    double alpha = (parameter - _knotVector[idx]) /
                   (_knotVector[idx + _degree] - _knotVector[idx]);

    double weightedX1 =
        CGAL::to_double(_controlPoints[idx - 1].x()) * _weights[idx - 1];
    double weightedY1 =
        CGAL::to_double(_controlPoints[idx - 1].y()) * _weights[idx - 1];
    double weightedX2 =
        CGAL::to_double(_controlPoints[idx].x()) * _weights[idx];
    double weightedY2 =
        CGAL::to_double(_controlPoints[idx].y()) * _weights[idx];
    double weight1 = _weights[idx - 1];
    double weight2 = _weights[idx];

    double newWeight = ((1 - alpha) * weight1) + (alpha * weight2);
    double newX =
        (((1 - alpha) * weightedX1) + (alpha * weightedX2)) / newWeight;
    double newY =
        (((1 - alpha) * weightedY1) + (alpha * weightedY2)) / newWeight;

    if (is3D()) {
      double weightedZ1 =
          CGAL::to_double(_controlPoints[idx - 1].z()) * _weights[idx - 1];
      double weightedZ2 =
          CGAL::to_double(_controlPoints[idx].z()) * _weights[idx];
      double newZ =
          (((1 - alpha) * weightedZ1) + (alpha * weightedZ2)) / newWeight;
      temp[index] = Point(newX, newY, newZ);
    } else {
      temp[index] = Point(newX, newY);
    }

    tempWeights[index] = newWeight;
  }

  // Add new control points
  for (size_t index = 0; index < _degree; ++index) {
    newControlPoints.push_back(temp[index]);
    newWeights.push_back(tempWeights[index]);
  }

  // Copy remaining control points
  for (size_t index = span; index < _controlPoints.size(); ++index) {
    newControlPoints.push_back(_controlPoints[index]);
    newWeights.push_back(_weights[index]);
  }

  // Update knot vector
  std::vector<double> newKnotVector;
  for (size_t index = 0; index <= span; ++index) {
    newKnotVector.push_back(_knotVector[index]);
  }
  newKnotVector.push_back(parameter);
  for (size_t index = span + 1; index < _knotVector.size(); ++index) {
    newKnotVector.push_back(_knotVector[index]);
  }

  // Update curve
  _controlPoints = newControlPoints;
  _weights       = newWeights;
  _knotVector    = newKnotVector;
}

void
NURBSCurve::refineKnots(const std::vector<double> &newKnots)
{
  for (double knot : newKnots) {
    insertKnot(knot);
  }
}

auto
NURBSCurve::createCircularArc(const Point &center, double radius,
                              double startAngle, double endAngle, bool is3D)
    -> std::unique_ptr<NURBSCurve>
{
  // Create a circular arc using NURBS representation
  double angleSpan = endAngle - startAngle;
  while (angleSpan < 0) {
    angleSpan += 2 * M_PI;
  }
  while (angleSpan > 2 * M_PI) {
    angleSpan -= 2 * M_PI;
  }

  int numSegments;
  if (angleSpan <= M_PI) {
    numSegments = 1;
  } else if (angleSpan <= 1.5 * M_PI) {
    numSegments = 2;
  } else {
    numSegments = 3;
  }

  std::vector<Point>  controlPoints;
  std::vector<double> weights;

  double segmentAngle = angleSpan / static_cast<double>(numSegments);
  double weightValue  = std::cos(segmentAngle / 2);

  for (int index = 0; index <= 2 * numSegments; ++index) {
    double angle =
        startAngle + ((static_cast<double>(index) * segmentAngle) / 2);
    double coordX = CGAL::to_double(center.x()) + (radius * std::cos(angle));
    double coordY = CGAL::to_double(center.y()) + (radius * std::sin(angle));

    if (is3D) {
      controlPoints.emplace_back(coordX, coordY, CGAL::to_double(center.z()));
    } else {
      controlPoints.emplace_back(coordX, coordY);
    }

    // Weights: 1 for endpoints, w for intermediate points
    weights.emplace_back((index % 2 == 0) ? 1.0 : weightValue);
  }

  // Create knot vector for degree 2
  std::vector<double> knots;
  knots.emplace_back(0.0);
  knots.emplace_back(0.0);
  knots.emplace_back(0.0);

  for (int segmentIndex = 1; segmentIndex < numSegments; ++segmentIndex) {
    double knotValue = static_cast<double>(segmentIndex) / numSegments;
    knots.emplace_back(knotValue);
    knots.emplace_back(knotValue);
  }

  knots.emplace_back(1.0);
  knots.emplace_back(1.0);
  knots.emplace_back(1.0);

  return std::make_unique<NURBSCurve>(controlPoints, 2, weights, knots);
}

auto
NURBSCurve::validateWeights() const -> bool
{
  if (_weights.empty() && !_controlPoints.empty()) {
    return false;
  }

  if (_weights.size() != _controlPoints.size()) {
    return false;
  }

  // All weights must be positive
  return std::all_of(_weights.begin(), _weights.end(),
                     [](double weightValue) { return weightValue > 0.0; });
}

} // namespace SFCGAL
