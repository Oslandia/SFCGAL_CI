#include "SFCGAL/CatmullRomCurve.h"
#include "SFCGAL/BSplineCurve.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/HermiteCurve.h"
#include <algorithm>
#include <cmath>

namespace SFCGAL {

CatmullRomCurve::CatmullRomCurve()
    : _tension(0.5), _boundaryType(BOUNDARY_CLAMPED),
      _parametrization(CENTRIPETAL)
{
}

CatmullRomCurve::CatmullRomCurve(const std::vector<Point> &points,
                                 double                    tension)
    : _points(points), _tension(tension), _boundaryType(BOUNDARY_CLAMPED),
      _parametrization(CENTRIPETAL)
{
  updateParameterValues();
}

CatmullRomCurve::CatmullRomCurve(const LineString &lineString, double tension)
    : _tension(tension), _boundaryType(BOUNDARY_CLAMPED),
      _parametrization(CENTRIPETAL)
{
  for (size_t index = 0; index < lineString.numPoints(); ++index) {
    _points.push_back(lineString.pointN(index));
  }
  updateParameterValues();
}

auto
CatmullRomCurve::geometryTypeId() const -> GeometryType
{
  return TYPE_CATMULLROMCURVE;
}

auto
CatmullRomCurve::geometryType() const -> std::string
{
  return "CatmullRomCurve";
}

auto
CatmullRomCurve::clone() const -> CatmullRomCurve *
{
  return new CatmullRomCurve(*this);
}

auto
CatmullRomCurve::isEmpty() const -> bool
{
  return _points.empty();
}

auto
CatmullRomCurve::is3D() const -> bool
{
  return !_points.empty() && _points[0].is3D();
}

auto
CatmullRomCurve::isMeasured() const -> bool
{
  return !_points.empty() && _points[0].isMeasured();
}

auto
CatmullRomCurve::point(size_t index) const -> const Point &
{
  BOOST_ASSERT(index < _points.size());
  return _points[index];
}

auto
CatmullRomCurve::point(size_t index) -> Point &
{
  BOOST_ASSERT(index < _points.size());
  return _points[index];
}

void
CatmullRomCurve::setPoint(size_t index, const Point &point)
{
  BOOST_ASSERT(index < _points.size());
  _points[index] = point;
  updateParameterValues();
}

void
CatmullRomCurve::addPoint(const Point &point)
{
  _points.push_back(point);
  updateParameterValues();
}

void
CatmullRomCurve::insertPoint(size_t index, const Point &point)
{
  BOOST_ASSERT(index <= _points.size());
  _points.insert(_points.begin() + static_cast<ptrdiff_t>(index), point);
  updateParameterValues();
}

void
CatmullRomCurve::removePoint(size_t index)
{
  BOOST_ASSERT(index < _points.size());
  _points.erase(_points.begin() + static_cast<ptrdiff_t>(index));
  updateParameterValues();
}

void
CatmullRomCurve::clear()
{
  _points.clear();
  _parameterValues.clear();
}

void
CatmullRomCurve::setTension(double tension)
{
  _tension = std::max(0.0, std::min(1.0, tension));
}

void
CatmullRomCurve::setBoundaryType(BoundaryType type)
{
  _boundaryType = type;
}

void
CatmullRomCurve::setParametrization(ParametrizationType type)
{
  _parametrization = type;
  updateParameterValues();
}

void
CatmullRomCurve::updateParameterValues() const
{
  _parameterValues.clear();

  if (_points.size() < 2) {
    return;
  }

  _parameterValues.emplace_back(0.0);

  if (_parametrization == UNIFORM) {
    // Uniform spacing
    for (size_t index = 1; index < _points.size(); ++index) {
      _parameterValues.emplace_back(static_cast<Kernel::FT>(index));
    }
  } else {
    // Chord length or centripetal
    Kernel::FT totalLength = 0.0;

    for (size_t index = 1; index < _points.size(); ++index) {
      auto deltaX = _points[index].x() - _points[index - 1].x();
      auto deltaY = _points[index].y() - _points[index - 1].y();
      auto deltaZ =
          is3D() ? (_points[index].z() - _points[index - 1].z()) : 0.0;

      auto segmentLength = std::sqrt(CGAL::to_double(
          (deltaX * deltaX) + (deltaY * deltaY) + (deltaZ * deltaZ)));

      if (_parametrization == CENTRIPETAL) {
        segmentLength = std::sqrt(
            CGAL::to_double(segmentLength)); // Square root for centripetal
      }

      totalLength += segmentLength;
      _parameterValues.push_back(totalLength);
    }

    // Normalize to [0, 1]
    if (totalLength > 0) {
      for (size_t index = 1; index < _parameterValues.size(); ++index) {
        _parameterValues[index] /= totalLength;
      }
    }
  }
}

auto
CatmullRomCurve::findSegment(double parameter) const -> size_t
{
  if (_points.size() < 2) {
    return 0;
  }

  if (_parameterValues.empty()) {
    updateParameterValues();
  }

  // Binary search for the segment
  size_t left  = 0;
  size_t right = _parameterValues.size() - 1;

  while (right - left > 1) {
    size_t mid = (left + right) / 2;
    if (parameter < _parameterValues[mid]) {
      right = mid;
    } else {
      left = mid;
    }
  }

  return left;
}

auto
CatmullRomCurve::catmullRomInterpolate(const Point &point0, const Point &point1,
                                       const Point &point2, const Point &point3,
                                       double parameter) const -> Point
{
  // Catmull-Rom matrix coefficients
  auto paramSquared = parameter * parameter;
  auto paramCubed   = paramSquared * parameter;

  // Apply tension parameter
  auto tau = _tension;

  // Catmull-Rom basis functions
  auto basis0 =
      ((-tau * paramCubed) + (2.0 * tau * paramSquared)) - (tau * parameter);
  auto basis1 = ((2.0 - tau) * paramCubed) + ((tau - 3.0) * paramSquared) + 1.0;
  auto basis2 = ((tau - 2.0) * paramCubed) +
                ((3.0 - 2.0 * tau) * paramSquared) + (tau * parameter);
  auto basis3 = (tau * paramCubed) - (tau * paramSquared);

  auto coordX = (basis0 * point0.x()) + (basis1 * point1.x()) +
                (basis2 * point2.x()) + (basis3 * point3.x());
  auto coordY = (basis0 * point0.y()) + (basis1 * point1.y()) +
                (basis2 * point2.y()) + (basis3 * point3.y());

  if (is3D()) {
    auto coordZ = (basis0 * point0.z()) + (basis1 * point1.z()) +
                  (basis2 * point2.z()) + (basis3 * point3.z());
    if (isMeasured()) {
      // Linear interpolation for M values
      auto measure =
          ((1.0 - parameter) * point1.m()) + (parameter * point2.m());
      return {coordX, coordY, coordZ, measure};
    }
    return {coordX, coordY, coordZ};
  }
  if (isMeasured()) {
    auto measure = ((1.0 - parameter) * point1.m()) + (parameter * point2.m());
    auto resultPoint = Point(coordX, coordY);
    resultPoint.setM(measure);
    return resultPoint;
  }

  return {coordX, coordY};
}

// Helper function to get boundary point for first segment
auto
CatmullRomCurve::getFirstSegmentBoundaryPoint() const -> Point
{
  switch (_boundaryType) {
  case BOUNDARY_ZERO:
    return _points[0]; // Duplicate first point

  case BOUNDARY_CLAMPED:
    // Reflect based on first tangent
    if (_points.size() > 1) {
      auto deltaX = _points[1].x() - _points[0].x();
      auto deltaY = _points[1].y() - _points[0].y();

      auto result = Point(_points[0].x() - deltaX, _points[0].y() - deltaY);
      if (is3D()) {
        auto deltaZ = _points[1].z() - _points[0].z();
        result      = Point(result.x(), result.y(), _points[0].z() - deltaZ);
      }
      return result;
    }
    return _points[0];

  case BOUNDARY_CYCLIC:
    return _points[_points.size() - 2];

  case BOUNDARY_NATURAL:
    // Natural boundary - extrapolate
    if (_points.size() > 2) {
      auto deltaX = _points[1].x() - _points[2].x();
      auto deltaY = _points[1].y() - _points[2].y();

      auto result = Point(_points[0].x() + deltaX, _points[0].y() + deltaY);
      if (is3D()) {
        auto deltaZ = _points[1].z() - _points[2].z();
        result      = Point(result.x(), result.y(), _points[0].z() + deltaZ);
      }
      return result;
    }
    return _points[0];
  }
  return _points[0]; // fallback
}

// Helper function to get boundary point for last segment
auto
CatmullRomCurve::getLastSegmentBoundaryPoint() const -> Point
{
  size_t pointCount = _points.size();

  switch (_boundaryType) {
  case BOUNDARY_ZERO:
    return _points[pointCount - 1]; // Duplicate last point

  case BOUNDARY_CLAMPED:
    // Reflect based on last tangent
    if (pointCount > 1) {
      auto deltaX = _points[pointCount - 1].x() - _points[pointCount - 2].x();
      auto deltaY = _points[pointCount - 1].y() - _points[pointCount - 2].y();

      auto result = Point(_points[pointCount - 1].x() + deltaX,
                          _points[pointCount - 1].y() + deltaY);
      if (is3D()) {
        auto deltaZ = _points[pointCount - 1].z() - _points[pointCount - 2].z();
        result =
            Point(result.x(), result.y(), _points[pointCount - 1].z() + deltaZ);
      }
      return result;
    }
    return _points[pointCount - 1];

  case BOUNDARY_CYCLIC:
    return _points[1];

  case BOUNDARY_NATURAL:
    // Natural boundary - extrapolate
    if (pointCount > 2) {
      auto deltaX = _points[pointCount - 2].x() - _points[pointCount - 3].x();
      auto deltaY = _points[pointCount - 2].y() - _points[pointCount - 3].y();

      auto result = Point(_points[pointCount - 1].x() + deltaX,
                          _points[pointCount - 1].y() + deltaY);
      if (is3D()) {
        auto deltaZ = _points[pointCount - 2].z() - _points[pointCount - 3].z();
        result =
            Point(result.x(), result.y(), _points[pointCount - 1].z() + deltaZ);
      }
      return result;
    }
    return _points[pointCount - 1];
  }
  return _points[pointCount - 1]; // fallback
}

void
CatmullRomCurve::getBoundaryPoints(Point &point0, Point &point3,
                                   size_t segment) const
{
  // Get phantom points for boundaries
  if (segment == 0) {
    point0 = getFirstSegmentBoundaryPoint();
  }

  if (segment == _points.size() - 2) {
    point3 = getLastSegmentBoundaryPoint();
  }
}

auto
CatmullRomCurve::evaluate(double parameter) const -> Point
{
  if (_points.empty()) {
    throw InappropriateGeometryException(
        "Cannot evaluate empty Catmull-Rom curve");
  }

  if (_points.size() == 1) {
    return _points[0];
  }

  if (_points.size() == 2) {
    // Linear interpolation for two points
    auto coordX =
        ((1.0 - parameter) * _points[0].x()) + (parameter * _points[1].x());
    auto coordY =
        ((1.0 - parameter) * _points[0].y()) + (parameter * _points[1].y());

    if (is3D()) {
      auto coordZ =
          ((1.0 - parameter) * _points[0].z()) + (parameter * _points[1].z());
      return {coordX, coordY, coordZ};
    }
    return {coordX, coordY};
  }

  auto bounds = parameterBounds();
  if (parameter < bounds.first || parameter > bounds.second) {
    throw InappropriateGeometryException("Parameter t is outside valid range");
  }

  // Handle endpoints exactly
  if (parameter <= 0.0) {
    return _points[0];
  }
  if (parameter >= 1.0) {
    return _points.back();
  }

  // Find the segment
  size_t segment = findSegment(parameter);
  if (segment >= _points.size() - 1) {
    segment = _points.size() - 2;
  }

  // Get the four control points
  Point point0;
  Point point1;
  Point point2;
  Point point3;

  point1 = _points[segment];
  point2 = _points[segment + 1];

  if (segment == 0) {
    getBoundaryPoints(point0, point3, segment);
    point3 = (segment + 2 < _points.size()) ? _points[segment + 2] : point3;
  } else if (segment == _points.size() - 2) {
    point0 = _points[segment - 1];
    getBoundaryPoints(point0, point3, segment);
  } else {
    point0 = _points[segment - 1];
    point3 = _points[segment + 2];
  }

  // Calculate local parameter
  auto param0 = _parameterValues[segment];
  auto param1 = _parameterValues[segment + 1];
  auto localT = (parameter - param0) / (param1 - param0);

  return catmullRomInterpolate(point0, point1, point2, point3,
                               CGAL::to_double(localT));
}

// Helper function to calculate first derivative for linear case
auto
CatmullRomCurve::calculateLinearDerivative() const -> Point
{
  auto deltaX = _points[1].x() - _points[0].x();
  auto deltaY = _points[1].y() - _points[0].y();

  if (is3D()) {
    auto deltaZ = _points[1].z() - _points[0].z();
    return {deltaX, deltaY, deltaZ};
  }
  return {deltaX, deltaY};
}

// Helper function to get control points for derivative calculation
void
CatmullRomCurve::getDerivativeControlPoints(size_t segment, Point &point0,
                                            Point &point1, Point &point2,
                                            Point &point3) const
{
  point1 = _points[segment];
  point2 = _points[segment + 1];

  if (segment == 0) {
    getBoundaryPoints(point0, point3, segment);
    point3 = (segment + 2 < _points.size()) ? _points[segment + 2] : point3;
  } else if (segment == _points.size() - 2) {
    point0 = _points[segment - 1];
    getBoundaryPoints(point0, point3, segment);
  } else {
    point0 = _points[segment - 1];
    point3 = _points[segment + 2];
  }
}

// Helper function to calculate Catmull-Rom derivative basis functions
auto
CatmullRomCurve::calculateDerivativeBasisFunctions(double localParam) const
    -> std::array<double, 4>
{
  auto paramSquared = localParam * localParam;
  auto tau          = _tension;

  return {((-3.0 * tau * paramSquared) + (4.0 * tau * localParam)) - tau,
          (3.0 * (2.0 - tau) * paramSquared) + (2.0 * (tau - 3.0) * localParam),
          (3.0 * (tau - 2.0) * paramSquared) +
              (2.0 * (3.0 - 2.0 * tau) * localParam) + tau,
          (3.0 * tau * paramSquared) - (2.0 * tau * localParam)};
}

auto
CatmullRomCurve::derivative(double parameter, unsigned int order) const -> Point
{
  if (order == 0) {
    return evaluate(parameter);
  }

  if (_points.size() < 2) {
    return is3D() ? Point(0, 0, 0) : Point(0, 0);
  }

  if (order == 1) {
    // Analytical first derivative
    if (_points.size() == 2) {
      return calculateLinearDerivative();
    }

    // Find segment
    size_t segment = findSegment(parameter);
    if (segment >= _points.size() - 1) {
      segment = _points.size() - 2;
    }

    // Get the four control points
    Point point0;
    Point point1;
    Point point2;
    Point point3;
    getDerivativeControlPoints(segment, point0, point1, point2, point3);

    // Calculate local parameter
    auto param0     = _parameterValues[segment];
    auto param1     = _parameterValues[segment + 1];
    auto localParam = (parameter - param0) / (param1 - param0);

    // Get derivative basis functions
    auto basisFunctions =
        calculateDerivativeBasisFunctions(CGAL::to_double(localParam));

    auto deltaX =
        (basisFunctions[0] * point0.x()) + (basisFunctions[1] * point1.x()) +
        (basisFunctions[2] * point2.x()) + (basisFunctions[3] * point3.x());
    auto deltaY =
        (basisFunctions[0] * point0.y()) + (basisFunctions[1] * point1.y()) +
        (basisFunctions[2] * point2.y()) + (basisFunctions[3] * point3.y());

    // Scale by inverse of segment parameter length
    auto scale = 1.0 / (param1 - param0);
    deltaX *= scale;
    deltaY *= scale;

    if (is3D()) {
      auto deltaZ =
          (basisFunctions[0] * point0.z()) + (basisFunctions[1] * point1.z()) +
          (basisFunctions[2] * point2.z()) + (basisFunctions[3] * point3.z());
      deltaZ *= scale;
      return {deltaX, deltaY, deltaZ};
    }

    return {deltaX, deltaY};
  }

  // Higher order derivatives using finite differences
  const auto stepSize = 1e-8;
  auto       bounds   = parameterBounds();

  Point derivative1 = (parameter - stepSize >= bounds.first)
                          ? derivative(parameter - stepSize, order - 1)
                          : derivative(parameter, order - 1);

  Point derivative2 = (parameter + stepSize <= bounds.second)
                          ? derivative(parameter + stepSize, order - 1)
                          : derivative(parameter, order - 1);

  auto deltaX = (derivative2.x() - derivative1.x()) / (2 * stepSize);
  auto deltaY = (derivative2.y() - derivative1.y()) / (2 * stepSize);

  if (is3D()) {
    auto deltaZ = (derivative2.z() - derivative1.z()) / (2 * stepSize);
    return {deltaX, deltaY, deltaZ};
  }

  return {deltaX, deltaY};
}

auto
CatmullRomCurve::toLineString(unsigned int numSegments) const
    -> std::unique_ptr<LineString>
{
  if (_points.empty()) {
    return std::make_unique<LineString>();
  }

  if (_points.size() == 1) {
    auto lineString = std::make_unique<LineString>();
    lineString->addPoint(_points[0]);
    return lineString;
  }

  numSegments = std::max(numSegments,
                         static_cast<unsigned int>(
                             _points.size() * 2)); // Ensure adequate sampling

  auto lineString = std::make_unique<LineString>();

  for (unsigned int index = 0; index <= numSegments; ++index) {
    auto parameter =
        static_cast<double>(index) / static_cast<double>(numSegments);
    lineString->addPoint(evaluate(parameter));
  }

  return lineString;
}

auto
CatmullRomCurve::parameterBounds() const -> std::pair<double, double>
{
  return std::make_pair(0.0, 1.0);
}

auto
CatmullRomCurve::isValid() const -> bool
{
  if (_points.empty()) {
    return true;
  }

  // Check dimensional consistency
  bool first3D = _points[0].is3D();
  bool firstM  = _points[0].isMeasured();

  for (size_t index = 1; index < _points.size(); ++index) {
    if (_points[index].is3D() != first3D ||
        _points[index].isMeasured() != firstM) {
      return false;
    }
  }

  // Check tension is valid - simplified using DeMorgan's theorem
  return (_tension >= 0.0 && _tension <= 1.0);
}

auto
CatmullRomCurve::toHermiteCurve() const -> std::unique_ptr<HermiteCurve>
{
  auto hermite = std::make_unique<HermiteCurve>();

  if (_points.empty()) {
    return hermite;
  }

  for (size_t index = 0; index < _points.size(); ++index) {
    HermiteCurve::ControlPoint controlPoint(_points[index]);

    // Calculate tangents
    Point tangent = getTangent(index);

    // Set handles based on tangent
    auto scale = 0.33; // Handle length factor

    controlPoint.inHandle = Point(_points[index].x() - (tangent.x() * scale),
                                  _points[index].y() - (tangent.y() * scale));

    controlPoint.outHandle = Point(_points[index].x() + (tangent.x() * scale),
                                   _points[index].y() + (tangent.y() * scale));

    if (is3D()) {
      controlPoint.inHandle =
          Point(controlPoint.inHandle.x(), controlPoint.inHandle.y(),
                _points[index].z() - (tangent.z() * scale));
      controlPoint.outHandle =
          Point(controlPoint.outHandle.x(), controlPoint.outHandle.y(),
                _points[index].z() + (tangent.z() * scale));
    }

    hermite->addControlPoint(controlPoint);
  }

  return hermite;
}

auto
CatmullRomCurve::getTangent(size_t index) const -> Point
{
  if (_points.size() < 2) {
    return is3D() ? Point(0, 0, 0) : Point(0, 0);
  }

  Point tangent;

  if (index == 0) {
    // First point
    tangent = Point((_points[1].x() - _points[0].x()) * _tension,
                    (_points[1].y() - _points[0].y()) * _tension);

    if (is3D()) {
      tangent = Point(tangent.x(), tangent.y(),
                      (_points[1].z() - _points[0].z()) * _tension);
    }
  } else if (index == _points.size() - 1) {
    // Last point
    size_t pointCount = _points.size();
    tangent           = Point(
        (_points[pointCount - 1].x() - _points[pointCount - 2].x()) * _tension,
        (_points[pointCount - 1].y() - _points[pointCount - 2].y()) * _tension);

    if (is3D()) {
      tangent =
          Point(tangent.x(), tangent.y(),
                (_points[pointCount - 1].z() - _points[pointCount - 2].z()) *
                    _tension);
    }
  } else {
    // Middle points
    tangent = Point(
        (_points[index + 1].x() - _points[index - 1].x()) * 0.5 * _tension,
        (_points[index + 1].y() - _points[index - 1].y()) * 0.5 * _tension);

    if (is3D()) {
      tangent = Point(tangent.x(), tangent.y(),
                      (_points[index + 1].z() - _points[index - 1].z()) * 0.5 *
                          _tension);
    }
  }

  return tangent;
}

auto
CatmullRomCurve::toBSplineCurve() const -> std::unique_ptr<BSplineCurve>
{
  // Convert to cubic B-spline
  // This requires solving for B-spline control points that interpolate our
  // points

  if (_points.size() < 2) {
    return std::make_unique<BSplineCurve>(_points, 3);
  }

  // For now, return a simple approximation
  // A proper implementation would solve the interpolation system
  std::vector<Point> controlPoints;

  // Add extra control points for cubic B-spline
  controlPoints.push_back(_points[0]);

  for (const auto &point : _points) {
    controlPoints.push_back(point);
  }

  controlPoints.push_back(_points.back());

  return std::make_unique<BSplineCurve>(controlPoints, 3);
}

} // namespace SFCGAL
