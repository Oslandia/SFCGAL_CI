#include "SFCGAL/BezierCurve.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryVisitor.h"
#include "SFCGAL/config.h"
#include <CGAL/number_utils.h>
#include <algorithm>
#include <cmath>

namespace SFCGAL {

BezierCurve::BezierCurve(const std::vector<Point> &controlPoints)
    : _controlPoints(controlPoints)
{
}

auto
BezierCurve::geometryTypeId() const -> GeometryType
{
  return TYPE_BEZIERCURVE;
}

auto
BezierCurve::geometryType() const -> std::string
{
  return "BezierCurve";
}

auto
BezierCurve::clone() const -> BezierCurve *
{
  return new BezierCurve(*this);
}

void
BezierCurve::accept(GeometryVisitor &visitor)
{
  visitor.visit(*this);
}

void
BezierCurve::accept(ConstGeometryVisitor &visitor) const
{
  visitor.visit(*this);
}

auto
BezierCurve::isEmpty() const -> bool
{
  return _controlPoints.empty();
}

auto
BezierCurve::is3D() const -> bool
{
  return !_controlPoints.empty() && _controlPoints[0].is3D();
}

auto
BezierCurve::isMeasured() const -> bool
{
  return !_controlPoints.empty() && _controlPoints[0].isMeasured();
}

auto
BezierCurve::degree() const -> unsigned int
{
  return _controlPoints.empty()
             ? 0
             : static_cast<unsigned int>(_controlPoints.size() - 1);
}

auto
BezierCurve::controlPointAt(size_t index) const -> const Point &
{
  BOOST_ASSERT(index < _controlPoints.size());
  return _controlPoints[index];
}

auto
BezierCurve::controlPointAt(size_t index) -> Point &
{
  BOOST_ASSERT(index < _controlPoints.size());
  return _controlPoints[index];
}

void
BezierCurve::setControlPoint(size_t index, const Point &point)
{
  BOOST_ASSERT(index < _controlPoints.size());
  _controlPoints[index] = point;
}

auto
BezierCurve::controlPoints() const -> std::vector<Point>
{
  return _controlPoints;
}

auto
BezierCurve::lerp(const Point &point1, const Point &point2, double parameter)
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

auto
BezierCurve::evaluate(double parameter) const -> Point
{
  if (_controlPoints.empty()) {
    BOOST_THROW_EXCEPTION(Exception("Cannot evaluate empty Bezier curve"));
  }

  if (parameter < 0.0 || parameter > 1.0) {
    BOOST_THROW_EXCEPTION(Exception("Parameter t must be in [0, 1]"));
  }

  // De Casteljau's algorithm
  std::vector<Point> points = _controlPoints;

  for (size_t level = 0; level < degree(); ++level) {
    for (size_t index = 0; index < points.size() - 1; ++index) {
      points[index] = lerp(points[index], points[index + 1], parameter);
    }
    points.pop_back();
  }

  return points[0];
}

auto
BezierCurve::derivative(double parameter, unsigned int order) const -> Point
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

  if (degree() < order || _controlPoints.size() < 2) {
    return {0.0, 0.0, 0.0, 0.0, coordType};
  }

  std::vector<Point> derivativePoints;
  auto               factor = static_cast<double>(degree());

  for (size_t index = 0; index < _controlPoints.size() - 1; ++index) {
    double deltaX = CGAL::to_double(_controlPoints[index + 1].x() -
                                    _controlPoints[index].x()) *
                    factor;
    double deltaY = CGAL::to_double(_controlPoints[index + 1].y() -
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

  if (derivativePoints.empty()) {
    return {0.0, 0.0, 0.0, 0.0, coordType};
  }

  BezierCurve derivativeCurve(derivativePoints);

  if (order == 1) {
    return derivativeCurve.evaluate(parameter);
  }

  return derivativeCurve.derivative(parameter, order - 1);
}

auto
BezierCurve::toLineString(unsigned int numSegments) const
    -> std::unique_ptr<LineString>
{
  if (_controlPoints.empty()) {
    return std::make_unique<LineString>();
  }

  if (numSegments < 2) {
    BOOST_THROW_EXCEPTION(Exception("Number of segments must be at least 2"));
  }

  auto lineString = std::make_unique<LineString>();

  for (unsigned int index = 0; index <= numSegments; ++index) {
    double parameter = static_cast<double>(index) / numSegments;
    lineString->addPoint(evaluate(parameter));
  }

  return lineString;
}

auto
BezierCurve::parameterBounds() const -> std::pair<double, double>
{
  return std::make_pair(0.0, 1.0);
}

auto
BezierCurve::isValid() const -> bool
{
  if (_controlPoints.empty()) {
    return true;
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

  return true;
}

void
BezierCurve::addControlPoint(const Point &point)
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
}

void
BezierCurve::removeControlPoint(size_t index)
{
  BOOST_ASSERT(index < _controlPoints.size());
  _controlPoints.erase(_controlPoints.begin() + static_cast<ptrdiff_t>(index));
}

void
BezierCurve::clear()
{
  _controlPoints.clear();
}

auto
BezierCurve::evaluateWithSteps(double parameter) const
    -> std::vector<std::vector<Point>>
{
  std::vector<std::vector<Point>> allSteps;

  if (_controlPoints.empty()) {
    return allSteps;
  }

  std::vector<Point> points = _controlPoints;
  allSteps.push_back(points);

  for (size_t level = 0; level < degree(); ++level) {
    for (size_t index = 0; index < points.size() - 1; ++index) {
      points[index] = lerp(points[index], points[index + 1], parameter);
    }
    points.pop_back();
    if (!points.empty()) {
      allSteps.push_back(points);
    }
  }

  return allSteps;
}

} // namespace SFCGAL
