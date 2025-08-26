#include "SFCGAL/BSplineCurve.h"
#include "SFCGAL/Exception.h"
#include <algorithm>
#include <cmath>

namespace SFCGAL {

BSplineCurve::BSplineCurve() : _degree(0) {}

BSplineCurve::BSplineCurve(const std::vector<Point> &controlPoints,
                           unsigned int              degree)
    : _controlPoints(controlPoints), _degree(degree)
{
  if (!_controlPoints.empty() && _degree >= _controlPoints.size()) {
    throw InappropriateGeometryException(
        "Degree must be less than the number of control points");
  }
  generateUniformKnotVector();
}

BSplineCurve::BSplineCurve(const std::vector<Point>  &controlPoints,
                           unsigned int               degree,
                           const std::vector<double> &knotVector)
    : _controlPoints(controlPoints), _degree(degree), _knotVector(knotVector)
{
  if (!_controlPoints.empty() && _degree >= _controlPoints.size()) {
    throw InappropriateGeometryException(
        "Degree must be less than the number of control points");
  }
  if (!validateKnotVector()) {
    throw InappropriateGeometryException("Invalid knot vector");
  }
}

BSplineCurve::BSplineCurve(const BSplineCurve &other) = default;

auto
BSplineCurve::operator=(const BSplineCurve &other) -> BSplineCurve & = default;

auto
BSplineCurve::geometryTypeId() const -> GeometryType
{
  return TYPE_BSPLINECURVE;
}

auto
BSplineCurve::geometryType() const -> std::string
{
  return "BSplineCurve";
}

auto
BSplineCurve::clone() const -> BSplineCurve *
{
  return new BSplineCurve(*this);
}

auto
BSplineCurve::isEmpty() const -> bool
{
  return _controlPoints.empty();
}

auto
BSplineCurve::is3D() const -> bool
{
  return !_controlPoints.empty() && _controlPoints[0].is3D();
}

auto
BSplineCurve::isMeasured() const -> bool
{
  return !_controlPoints.empty() && _controlPoints[0].isMeasured();
}

auto
BSplineCurve::degree() const -> unsigned int
{
  return _degree;
}

auto
BSplineCurve::numControlPoints() const -> size_t
{
  return _controlPoints.size();
}

auto
BSplineCurve::controlPointAt(size_t index) const -> const Point &
{
  BOOST_ASSERT(index < _controlPoints.size());
  return _controlPoints[index];
}

auto
BSplineCurve::controlPointAt(size_t index) -> Point &
{
  BOOST_ASSERT(index < _controlPoints.size());
  return _controlPoints[index];
}

void
BSplineCurve::setControlPoint(size_t index, const Point &point)
{
  BOOST_ASSERT(index < _controlPoints.size());
  _controlPoints[index] = point;
}

auto
BSplineCurve::controlPoints() const -> std::vector<Point>
{
  return _controlPoints;
}

void
BSplineCurve::generateUniformKnotVector()
{
  if (_controlPoints.empty()) {
    _knotVector.clear();
    return;
  }

  size_t pointCount = _controlPoints.size();
  size_t knotCount  = pointCount + _degree + 1;
  _knotVector.resize(knotCount);

  // Set multiplicity at the beginning
  for (size_t index = 0; index <= _degree; ++index) {
    _knotVector[index] = 0.0;
  }

  // Set internal knots uniformly
  size_t segments = pointCount - _degree;
  for (size_t index = _degree + 1; index < knotCount - _degree - 1; ++index) {
    _knotVector[index] =
        static_cast<double>(index - _degree) / static_cast<double>(segments);
  }

  // Set multiplicity at the end
  for (size_t index = knotCount - _degree - 1; index < knotCount; ++index) {
    _knotVector[index] = 1.0;
  }
}

auto
BSplineCurve::findSpan(double parameter) const -> size_t
{
  size_t pointCount = _controlPoints.size() - 1;

  // Special cases
  if (parameter >= _knotVector[pointCount + 1]) {
    return pointCount;
  }
  if (parameter <= _knotVector[_degree]) {
    return _degree;
  }

  // Binary search
  size_t low  = _degree;
  size_t high = pointCount + 1;
  auto   mid  = (low + high) / 2;

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
BSplineCurve::evaluate(double parameter) const -> Point
{
  if (_controlPoints.empty()) {
    throw InappropriateGeometryException(
        "Cannot evaluate empty B-spline curve");
  }

  auto bounds = parameterBounds();
  if (parameter < bounds.first || parameter > bounds.second) {
    throw InappropriateGeometryException("Parameter t is outside valid range");
  }

  // De Boor's algorithm
  size_t             span = findSpan(parameter);
  std::vector<Point> temp(_degree + 1);

  // Initialize with relevant control points
  for (size_t pointIndex = 0; pointIndex <= _degree; ++pointIndex) {
    temp[pointIndex] = _controlPoints[span - _degree + pointIndex];
  }

  // Apply de Boor's algorithm
  for (size_t level = 1; level <= _degree; ++level) {
    for (size_t pointIndex = _degree; pointIndex >= level; --pointIndex) {
      auto alpha = (parameter - _knotVector[span - _degree + pointIndex]) /
                   (_knotVector[span + pointIndex - level + 1] -
                    _knotVector[span - _degree + pointIndex]);

      auto coordX = ((1 - alpha) * CGAL::to_double(temp[pointIndex - 1].x())) +
                    (alpha * CGAL::to_double(temp[pointIndex].x()));
      auto coordY = ((1 - alpha) * CGAL::to_double(temp[pointIndex - 1].y())) +
                    (alpha * CGAL::to_double(temp[pointIndex].y()));

      if (is3D()) {
        auto coordZ =
            ((1 - alpha) * CGAL::to_double(temp[pointIndex - 1].z())) +
            (alpha * CGAL::to_double(temp[pointIndex].z()));
        if (isMeasured()) {
          auto measure = ((1 - alpha) * temp[pointIndex - 1].m()) +
                         (alpha * temp[pointIndex].m());
          temp[pointIndex] = Point(coordX, coordY, coordZ, measure);
        } else {
          temp[pointIndex] = Point(coordX, coordY, coordZ);
        }
      } else if (isMeasured()) {
        temp[pointIndex] = Point(coordX, coordY);
        auto measure     = ((1 - alpha) * temp[pointIndex - 1].m()) +
                       (alpha * temp[pointIndex].m());
        temp[pointIndex].setM(measure);
      } else {
        temp[pointIndex] = Point(coordX, coordY);
      }
    }
  }

  return temp[_degree];
}

auto
BSplineCurve::derivative(double parameter, unsigned int order) const -> Point
{
  if (order == 0) {
    return evaluate(parameter);
  }

  if (_degree < order) {
    if (is3D()) {
      return {0, 0, 0};
    }
    return {0, 0};
  }

  // Create derivative B-spline curve
  std::vector<Point> derivativePoints;
  for (size_t index = 0; index < _controlPoints.size() - 1; ++index) {
    auto factor = static_cast<double>(_degree) /
                  (_knotVector[index + _degree + 1] - _knotVector[index + 1]);

    auto deltaX = CGAL::to_double(_controlPoints[index + 1].x() -
                                  _controlPoints[index].x()) *
                  factor;
    auto deltaY = CGAL::to_double(_controlPoints[index + 1].y() -
                                  _controlPoints[index].y()) *
                  factor;

    if (is3D()) {
      auto deltaZ = CGAL::to_double(_controlPoints[index + 1].z() -
                                    _controlPoints[index].z()) *
                    factor;
      derivativePoints.emplace_back(deltaX, deltaY, deltaZ);
    } else {
      derivativePoints.emplace_back(deltaX, deltaY);
    }
  }

  // Create derivative knot vector
  std::vector<double> derivativeKnots(_knotVector.begin() + 1,
                                      _knotVector.end() - 1);

  BSplineCurve derivativeCurve(derivativePoints, _degree - 1, derivativeKnots);

  if (order == 1) {
    return derivativeCurve.evaluate(parameter);
  }
  return derivativeCurve.derivative(parameter, order - 1);
}

auto
BSplineCurve::toLineString(unsigned int numSegments) const
    -> std::unique_ptr<LineString>
{
  if (_controlPoints.empty()) {
    return std::make_unique<LineString>();
  }

  if (numSegments < 2) {
    throw InappropriateGeometryException(
        "Number of segments must be at least 2");
  }

  auto lineString = std::make_unique<LineString>();
  auto bounds     = parameterBounds();

  for (unsigned int index = 0; index <= numSegments; ++index) {
    auto parameter = bounds.first + ((bounds.second - bounds.first) *
                                     static_cast<double>(index) /
                                     static_cast<double>(numSegments));
    lineString->addPoint(evaluate(parameter));
  }

  return lineString;
}

auto
BSplineCurve::parameterBounds() const -> std::pair<double, double>
{
  if (_knotVector.empty()) {
    return std::make_pair(0.0, 1.0);
  }
  return std::make_pair(_knotVector[_degree],
                        _knotVector[_controlPoints.size()]);
}

auto
BSplineCurve::isValid() const -> bool
{
  if (_controlPoints.empty()) {
    return true;
  }

  // Check degree validity
  if (_degree >= _controlPoints.size()) {
    return false;
  }

  // Check dimensional consistency
  bool first3D = _controlPoints[0].is3D();
  bool firstM  = _controlPoints[0].isMeasured();

  for (size_t index = 1; index < _controlPoints.size(); ++index) {
    if (_controlPoints[index].is3D() != first3D) {
      return false;
    }
    if (_controlPoints[index].isMeasured() != firstM) {
      return false;
    }
  }

  return validateKnotVector();
}

auto
BSplineCurve::knotVector() const -> const std::vector<double> &
{
  return _knotVector;
}

void
BSplineCurve::setKnotVector(const std::vector<double> &knots)
{
  _knotVector = knots;
  if (!validateKnotVector()) {
    throw InappropriateGeometryException("Invalid knot vector");
  }
}

auto
BSplineCurve::validateKnotVector() const -> bool
{
  if (_controlPoints.empty()) {
    return _knotVector.empty();
  }

  size_t expectedSize = _controlPoints.size() + _degree + 1;
  if (_knotVector.size() != expectedSize) {
    return false;
  }

  // Check non-decreasing order
  for (size_t index = 1; index < _knotVector.size(); ++index) {
    if (_knotVector[index] < _knotVector[index - 1]) {
      return false;
    }
  }

  return true;
}

} // namespace SFCGAL
