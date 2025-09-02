#include "SFCGAL/BSplineCurve.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryVisitor.h"
#include "SFCGAL/config.h"
#include <CGAL/number_utils.h>
#include <algorithm>
#include <cmath>

namespace SFCGAL {

BSplineCurve::BSplineCurve() : _degree(0) {}

BSplineCurve::BSplineCurve(const std::vector<Point> &controlPoints,
                           unsigned int              degree)
    : _controlPoints(controlPoints), _degree(degree)
{
  if (!_controlPoints.empty() && _degree >= _controlPoints.size()) {
    BOOST_THROW_EXCEPTION(
        Exception("Degree must be less than the number of control points"));
  }
  generateUniformKnotVector();
}

BSplineCurve::BSplineCurve(const std::vector<Point>  &controlPoints,
                           unsigned int               degree,
                           const std::vector<double> &knotVector)
    : _controlPoints(controlPoints), _degree(degree), _knotVector(knotVector)
{
  if (!_controlPoints.empty() && _degree >= _controlPoints.size()) {
    BOOST_THROW_EXCEPTION(
        Exception("Degree must be less than the number of control points"));
  }
  if (!validateKnotVector()) {
    BOOST_THROW_EXCEPTION(Exception("Invalid knot vector"));
  }
}

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

void
BSplineCurve::accept(GeometryVisitor &visitor)
{
  return visitor.visit(*this);
}

void
BSplineCurve::accept(ConstGeometryVisitor &visitor) const
{
  return visitor.visit(*this);
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

void
BSplineCurve::setDegree(unsigned int newDegree)
{
  if (!_controlPoints.empty() && newDegree >= _controlPoints.size()) {
    BOOST_THROW_EXCEPTION(
        Exception("Degree must be less than the number of control points"));
  }

  _degree = newDegree;

  // Regenerate knot vector with new degree
  if (!_knotVector.empty()) {
    generateUniformKnotVector();
  }
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

auto
BSplineCurve::lerp(const Point &point1, const Point &point2, double parameter)
    -> Point
{
  double coordX = (CGAL::to_double(point1.x()) * (1 - parameter)) +
                  (CGAL::to_double(point2.x()) * parameter);
  double coordY = (CGAL::to_double(point1.y()) * (1 - parameter)) +
                  (CGAL::to_double(point2.y()) * parameter);

  bool has3D = point1.is3D() && point2.is3D();
  bool hasM  = point1.isMeasured() && point2.isMeasured();

  CoordinateType coordType;
  double         coordZ  = 0.0;
  double         measure = 0.0;

  if (has3D && hasM) {
    coordType = COORDINATE_XYZM;
    coordZ    = (CGAL::to_double(point1.z()) * (1 - parameter)) +
             (CGAL::to_double(point2.z()) * parameter);
    measure = (point1.m() * (1 - parameter)) + (point2.m() * parameter);
  } else if (has3D) {
    coordType = COORDINATE_XYZ;
    coordZ    = (CGAL::to_double(point1.z()) * (1 - parameter)) +
             (CGAL::to_double(point2.z()) * parameter);
  } else if (hasM) {
    coordType = COORDINATE_XYM;
    measure   = (point1.m() * (1 - parameter)) + (point2.m() * parameter);
  } else {
    coordType = COORDINATE_XY;
  }

  return {coordX, coordY, coordZ, measure, coordType};
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
    BOOST_THROW_EXCEPTION(Exception("Cannot evaluate empty B-spline curve"));
  }

  auto bounds = parameterBounds();
  if (parameter < bounds.first || parameter > bounds.second) {
    BOOST_THROW_EXCEPTION(Exception("Parameter t is outside valid range"));
  }

  // De Boor's algorithm
  size_t             span = findSpan(parameter);
  std::vector<Point> temp(_degree + 1);

  // Initialize with relevant control points
  for (size_t pointIndex = 0; pointIndex <= _degree; ++pointIndex) {
    temp[pointIndex] = _controlPoints[span - _degree + pointIndex];
  }

  // Apply de Boor's algorithm using lerp
  for (size_t level = 1; level <= _degree; ++level) {
    for (size_t pointIndex = _degree; pointIndex >= level; --pointIndex) {
      auto alpha = (parameter - _knotVector[span - _degree + pointIndex]) /
                   (_knotVector[span + pointIndex - level + 1] -
                    _knotVector[span - _degree + pointIndex]);

      temp[pointIndex] = lerp(temp[pointIndex - 1], temp[pointIndex], alpha);
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

  CoordinateType coordType;
  if (is3D() && isMeasured()) {
    coordType = COORDINATE_XYZM;
  } else if (is3D()) {
    coordType = COORDINATE_XYZ;
  } else if (isMeasured()) {
    coordType = COORDINATE_XYM;
  } else {
    coordType = COORDINATE_XY;
  }

  if (_degree < order) {
    return Point(0.0, 0.0, 0.0, 0.0, coordType);
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

    double deltaZ = 0.0;
    double deltaM = 0.0;

    if (is3D()) {
      deltaZ = CGAL::to_double(_controlPoints[index + 1].z() -
                               _controlPoints[index].z()) *
               factor;
    }

    if (isMeasured()) {
      deltaM =
          (_controlPoints[index + 1].m() - _controlPoints[index].m()) * factor;
    }

    derivativePoints.emplace_back(deltaX, deltaY, deltaZ, deltaM, coordType);
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
    BOOST_THROW_EXCEPTION(Exception("Number of segments must be at least 2"));
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
    BOOST_THROW_EXCEPTION(Exception("Invalid knot vector"));
  }
}

void
BSplineCurve::addControlPoint(const Point &point)
{
  if (!_controlPoints.empty()) {
    // Ensure dimensional consistency
    bool needsZ = _controlPoints[0].is3D();
    bool needsM = _controlPoints[0].isMeasured();

    if (point.is3D() != needsZ || point.isMeasured() != needsM) {
      BOOST_THROW_EXCEPTION(
          Exception("Control point dimensions must be consistent"));
    }
  }

  _controlPoints.push_back(point);

  // Regenerate knot vector if uniform
  if (!_knotVector.empty()) {
    generateUniformKnotVector();
  }
}

void
BSplineCurve::removeControlPoint(size_t index)
{
  BOOST_ASSERT(index < _controlPoints.size());
  _controlPoints.erase(_controlPoints.begin() + static_cast<ptrdiff_t>(index));

  // Regenerate knot vector if uniform
  if (!_knotVector.empty()) {
    generateUniformKnotVector();
  }
}

void
BSplineCurve::clear()
{
  _controlPoints.clear();
  _knotVector.clear();
  _degree = 0;
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

auto
BSplineCurve::evaluateWithBasisFunctions(double parameter) const
    -> std::pair<Point, std::vector<double>>
{
  if (_controlPoints.empty()) {
    BOOST_THROW_EXCEPTION(Exception("Cannot evaluate empty B-spline curve"));
  }

  auto bounds = parameterBounds();
  if (parameter < bounds.first || parameter > bounds.second) {
    BOOST_THROW_EXCEPTION(Exception("Parameter t is outside valid range"));
  }

  size_t span = findSpan(parameter);

  // Compute basis functions
  std::vector<double> basisFunctions(_degree + 1);
  std::vector<double> left(_degree + 1);
  std::vector<double> right(_degree + 1);

  basisFunctions[0] = 1.0;

  for (size_t j = 1; j <= _degree; ++j) {
    left[j]      = parameter - _knotVector[span + 1 - j];
    right[j]     = _knotVector[span + j] - parameter;
    double saved = 0.0;

    for (size_t r = 0; r < j; ++r) {
      double temp       = basisFunctions[r] / (right[r + 1] + left[j - r]);
      basisFunctions[r] = saved + right[r + 1] * temp;
      saved             = left[j - r] * temp;
    }
    basisFunctions[j] = saved;
  }

  // Compute curve point
  CoordinateType coordType;
  if (is3D() && isMeasured()) {
    coordType = COORDINATE_XYZM;
  } else if (is3D()) {
    coordType = COORDINATE_XYZ;
  } else if (isMeasured()) {
    coordType = COORDINATE_XYM;
  } else {
    coordType = COORDINATE_XY;
  }

  double coordX = 0.0, coordY = 0.0, coordZ = 0.0, measure = 0.0;

  for (size_t i = 0; i <= _degree; ++i) {
    const Point &cp = _controlPoints[span - _degree + i];
    coordX += CGAL::to_double(cp.x()) * basisFunctions[i];
    coordY += CGAL::to_double(cp.y()) * basisFunctions[i];

    if (is3D()) {
      coordZ += CGAL::to_double(cp.z()) * basisFunctions[i];
    }
    if (isMeasured()) {
      measure += cp.m() * basisFunctions[i];
    }
  }

  return std::make_pair(Point(coordX, coordY, coordZ, measure, coordType),
                        basisFunctions);
}

} // namespace SFCGAL
