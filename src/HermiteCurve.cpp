#include "SFCGAL/HermiteCurve.h"
#include "SFCGAL/BezierCurve.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/config.h"
#include <CGAL/number_utils.h>
#include <algorithm>
#include <cmath>

namespace SFCGAL {

HermiteCurve::HermiteCurve() = default;

HermiteCurve::HermiteCurve(const std::vector<ControlPoint> &controlPoints)
    : _controlPoints(controlPoints)
{
}

HermiteCurve::HermiteCurve(const HermiteCurve &other) = default;

auto
HermiteCurve::operator=(const HermiteCurve &other) -> HermiteCurve &
{
  if (this != &other) {
    Curve::operator=(other);
    _controlPoints = other._controlPoints;
  }
  return *this;
}

auto
HermiteCurve::geometryTypeId() const -> GeometryType
{
  return TYPE_HERMITECURVE;
}

auto
HermiteCurve::geometryType() const -> std::string
{
  return "HermiteCurve";
}

auto
HermiteCurve::clone() const -> HermiteCurve *
{
  return new HermiteCurve(*this);
}

auto
HermiteCurve::isEmpty() const -> bool
{
  return _controlPoints.empty();
}

auto
HermiteCurve::is3D() const -> bool
{
  return !_controlPoints.empty() && _controlPoints[0].position.is3D();
}

auto
HermiteCurve::isMeasured() const -> bool
{
  return !_controlPoints.empty() && _controlPoints[0].position.isMeasured();
}

auto
HermiteCurve::getControlPoint(size_t index) const
    -> const HermiteCurve::ControlPoint &
{
  BOOST_ASSERT(index < _controlPoints.size());
  return _controlPoints[index];
}

auto
HermiteCurve::getControlPoint(size_t index) -> HermiteCurve::ControlPoint &
{
  BOOST_ASSERT(index < _controlPoints.size());
  return _controlPoints[index];
}

void
HermiteCurve::setControlPoint(size_t index, const ControlPoint &controlPoint)
{
  BOOST_ASSERT(index < _controlPoints.size());
  _controlPoints[index] = controlPoint;
}

void
HermiteCurve::addControlPoint(const ControlPoint &controlPoint)
{
  _controlPoints.push_back(controlPoint);
}

void
HermiteCurve::setInHandle(size_t index, const Point &handle)
{
  BOOST_ASSERT(index < _controlPoints.size());
  _controlPoints[index].inHandle = handle;
}

void
HermiteCurve::setOutHandle(size_t index, const Point &handle)
{
  BOOST_ASSERT(index < _controlPoints.size());
  _controlPoints[index].outHandle = handle;
}

auto
HermiteCurve::getInHandle(size_t index) const -> const Point &
{
  BOOST_ASSERT(index < _controlPoints.size());
  return _controlPoints[index].inHandle;
}

auto
HermiteCurve::getOutHandle(size_t index) const -> const Point &
{
  BOOST_ASSERT(index < _controlPoints.size());
  return _controlPoints[index].outHandle;
}

// Helper function to calculate smooth handles for the first control point
void
HermiteCurve::calculateFirstPointHandles(size_t index)
{
  if (_controlPoints.size() <= 1) {
    return;
  }

  Point &position = _controlPoints[index].position;
  Point &next     = _controlPoints[1].position;

  double deltaX = CGAL::to_double(next.x() - position.x());
  double deltaY = CGAL::to_double(next.y() - position.y());
  double deltaZ = is3D() ? CGAL::to_double(next.z() - position.z()) : 0.0;

  // Place handles at 1/3 distance
  _controlPoints[index].outHandle = {
      CGAL::to_double(position.x()) + (deltaX / 3.0),
      CGAL::to_double(position.y()) + (deltaY / 3.0),
      is3D() ? CGAL::to_double(position.z()) + (deltaZ / 3.0)
             : CGAL::to_double(position.z())};
  _controlPoints[index].inHandle = {
      CGAL::to_double(position.x()) - (deltaX / 6.0),
      CGAL::to_double(position.y()) - (deltaY / 6.0),
      is3D() ? CGAL::to_double(position.z()) - (deltaZ / 6.0)
             : CGAL::to_double(position.z())};
}

// Helper function to calculate smooth handles for the last control point
void
HermiteCurve::calculateLastPointHandles(size_t index)
{
  Point &position = _controlPoints[index].position;
  Point &previous = _controlPoints[index - 1].position;

  double deltaX = CGAL::to_double(position.x() - previous.x());
  double deltaY = CGAL::to_double(position.y() - previous.y());
  double deltaZ = is3D() ? CGAL::to_double(position.z() - previous.z()) : 0.0;

  _controlPoints[index].inHandle = {
      CGAL::to_double(position.x()) - (deltaX / 3.0),
      CGAL::to_double(position.y()) - (deltaY / 3.0),
      is3D() ? CGAL::to_double(position.z()) - (deltaZ / 3.0)
             : CGAL::to_double(position.z())};
  _controlPoints[index].outHandle = {
      CGAL::to_double(position.x()) + (deltaX / 6.0),
      CGAL::to_double(position.y()) + (deltaY / 6.0),
      is3D() ? CGAL::to_double(position.z()) + (deltaZ / 6.0)
             : CGAL::to_double(position.z())};
}

// Helper function to calculate distances for middle point handle calculations
auto
HermiteCurve::calculateNeighborDistances(size_t index) const
    -> std::pair<double, double>
{
  const Point &position = _controlPoints[index].position;
  const Point &previous = _controlPoints[index - 1].position;
  const Point &next     = _controlPoints[index + 1].position;

  double previousDistance = std::sqrt(
      std::pow(CGAL::to_double(position.x() - previous.x()), 2) +
      std::pow(CGAL::to_double(position.y() - previous.y()), 2) +
      (is3D() ? std::pow(CGAL::to_double(position.z() - previous.z()), 2)
              : 0.0));

  double nextDistance = std::sqrt(
      std::pow(CGAL::to_double(next.x() - position.x()), 2) +
      std::pow(CGAL::to_double(next.y() - position.y()), 2) +
      (is3D() ? std::pow(CGAL::to_double(next.z() - position.z()), 2) : 0.0));

  return {previousDistance, nextDistance};
}

// Helper function to calculate smooth handles for middle control points
void
HermiteCurve::calculateMiddlePointHandles(size_t index)
{
  Point &position = _controlPoints[index].position;
  Point &previous = _controlPoints[index - 1].position;
  Point &next     = _controlPoints[index + 1].position;

  double deltaX = CGAL::to_double(next.x() - previous.x());
  double deltaY = CGAL::to_double(next.y() - previous.y());
  double deltaZ = is3D() ? CGAL::to_double(next.z() - previous.z()) : 0.0;

  auto [previousDistance, nextDistance] = calculateNeighborDistances(index);
  double factor                         = 0.25; // Handle length factor

  _controlPoints[index].inHandle = {
      CGAL::to_double(position.x()) - ((deltaX * factor * previousDistance) /
                                       (previousDistance + nextDistance)),
      CGAL::to_double(position.y()) - ((deltaY * factor * previousDistance) /
                                       (previousDistance + nextDistance)),
      is3D() ? CGAL::to_double(position.z()) -
                   ((deltaZ * factor * previousDistance) /
                    (previousDistance + nextDistance))
             : CGAL::to_double(position.z())};

  _controlPoints[index].outHandle = {
      CGAL::to_double(position.x()) + ((deltaX * factor * nextDistance) /
                                       (previousDistance + nextDistance)),
      CGAL::to_double(position.y()) + ((deltaY * factor * nextDistance) /
                                       (previousDistance + nextDistance)),
      is3D()
          ? CGAL::to_double(position.z()) + ((deltaZ * factor * nextDistance) /
                                             (previousDistance + nextDistance))
          : CGAL::to_double(position.z())};
}

void
HermiteCurve::autoSmoothHandles(size_t index)
{
  if (_controlPoints.size() < 2 || index >= _controlPoints.size()) {
    return;
  }

  if (index == 0) {
    calculateFirstPointHandles(index);
  } else if (index == _controlPoints.size() - 1) {
    calculateLastPointHandles(index);
  } else {
    calculateMiddlePointHandles(index);
  }
}

void
HermiteCurve::autoSmoothAllHandles()
{
  for (size_t pointIndex = 0; pointIndex < _controlPoints.size();
       ++pointIndex) {
    autoSmoothHandles(pointIndex);
  }
}

auto
HermiteCurve::hermiteInterpolate(const ControlPoint &point1,
                                 const ControlPoint &point2,
                                 double              parameter) const -> Point
{
  // Hermite basis functions
  double paramSquared = parameter * parameter;
  double paramCubed   = paramSquared * parameter;

  double hermiteBasis1 =
      (2 * paramCubed) - (3 * paramSquared) + 1; // Position at point1
  double hermiteBasis2 =
      (-2 * paramCubed) + (3 * paramSquared); // Position at point2
  double hermiteBasis3 =
      paramCubed - (2 * paramSquared) + parameter;  // Tangent at point1
  double hermiteBasis4 = paramCubed - paramSquared; // Tangent at point2

  // Calculate tangent vectors from handles
  double tangent1X =
      CGAL::to_double(point1.outHandle.x() - point1.position.x()) * 3.0;
  double tangent1Y =
      CGAL::to_double(point1.outHandle.y() - point1.position.y()) * 3.0;
  double tangent2X =
      CGAL::to_double(point2.position.x() - point2.inHandle.x()) * 3.0;
  double tangent2Y =
      CGAL::to_double(point2.position.y() - point2.inHandle.y()) * 3.0;

  double coordX = (hermiteBasis1 * CGAL::to_double(point1.position.x())) +
                  (hermiteBasis2 * CGAL::to_double(point2.position.x())) +
                  (hermiteBasis3 * tangent1X) + (hermiteBasis4 * tangent2X);
  double coordY = (hermiteBasis1 * CGAL::to_double(point1.position.y())) +
                  (hermiteBasis2 * CGAL::to_double(point2.position.y())) +
                  (hermiteBasis3 * tangent1Y) + (hermiteBasis4 * tangent2Y);

  if (is3D()) {
    double tangent1Z =
        CGAL::to_double(point1.outHandle.z() - point1.position.z()) * 3.0;
    double tangent2Z =
        CGAL::to_double(point2.position.z() - point2.inHandle.z()) * 3.0;
    double coordZ = (hermiteBasis1 * CGAL::to_double(point1.position.z())) +
                    (hermiteBasis2 * CGAL::to_double(point2.position.z())) +
                    (hermiteBasis3 * tangent1Z) + (hermiteBasis4 * tangent2Z);

    if (isMeasured()) {
      double measure = (hermiteBasis1 * point1.position.m()) +
                       (hermiteBasis2 * point2.position.m());
      return {coordX, coordY, coordZ, measure};
    }
    return {coordX, coordY, coordZ};
  }

  if (isMeasured()) {
    double measure = (hermiteBasis1 * point1.position.m()) +
                     (hermiteBasis2 * point2.position.m());
    return {Kernel::FT(coordX), Kernel::FT(coordY), Kernel::FT(NaN()), measure};
  }

  return {coordX, coordY};
}

auto
HermiteCurve::hermiteTangent(const ControlPoint &point1,
                             const ControlPoint &point2, double parameter) const
    -> Point
{
  // Derivatives of Hermite basis functions
  double paramSquared = parameter * parameter;

  double derivativeBasis1 = (6 * paramSquared) - (6 * parameter);
  double derivativeBasis2 = (-6 * paramSquared) + (6 * parameter);
  double derivativeBasis3 = (3 * paramSquared) - (4 * parameter) + 1;
  double derivativeBasis4 = (3 * paramSquared) - (2 * parameter);

  // Calculate tangent vectors from handles
  double tangent1X =
      CGAL::to_double(point1.outHandle.x() - point1.position.x()) * 3.0;
  double tangent1Y =
      CGAL::to_double(point1.outHandle.y() - point1.position.y()) * 3.0;
  double tangent2X =
      CGAL::to_double(point2.position.x() - point2.inHandle.x()) * 3.0;
  double tangent2Y =
      CGAL::to_double(point2.position.y() - point2.inHandle.y()) * 3.0;

  double deltaX = (derivativeBasis1 * CGAL::to_double(point1.position.x())) +
                  (derivativeBasis2 * CGAL::to_double(point2.position.x())) +
                  (derivativeBasis3 * tangent1X) +
                  (derivativeBasis4 * tangent2X);
  double deltaY = (derivativeBasis1 * CGAL::to_double(point1.position.y())) +
                  (derivativeBasis2 * CGAL::to_double(point2.position.y())) +
                  (derivativeBasis3 * tangent1Y) +
                  (derivativeBasis4 * tangent2Y);

  if (is3D()) {
    double tangent1Z =
        CGAL::to_double(point1.outHandle.z() - point1.position.z()) * 3.0;
    double tangent2Z =
        CGAL::to_double(point2.position.z() - point2.inHandle.z()) * 3.0;
    double deltaZ = (derivativeBasis1 * CGAL::to_double(point1.position.z())) +
                    (derivativeBasis2 * CGAL::to_double(point2.position.z())) +
                    (derivativeBasis3 * tangent1Z) +
                    (derivativeBasis4 * tangent2Z);
    return {deltaX, deltaY, deltaZ};
  }

  return {deltaX, deltaY};
}

auto
HermiteCurve::evaluate(double parameter) const -> Point
{
  if (_controlPoints.empty()) {
    BOOST_THROW_EXCEPTION(Exception("Cannot evaluate empty Hermite curve"));
  }

  if (_controlPoints.size() == 1) {
    return _controlPoints[0].position;
  }

  auto bounds = parameterBounds();
  if (parameter < bounds.first || parameter > bounds.second) {
    BOOST_THROW_EXCEPTION(Exception("Parameter t is outside valid range"));
  }

  // Find the segment
  size_t numSegments   = _controlPoints.size() - 1;
  double segmentLength = 1.0 / static_cast<double>(numSegments);
  auto   segment       = static_cast<size_t>(parameter / segmentLength);

  if (segment >= numSegments) {
    segment = numSegments - 1;
  }

  // Local parameter within segment
  double localParameter =
      (parameter - (static_cast<double>(segment) * segmentLength)) /
      segmentLength;

  return hermiteInterpolate(_controlPoints[segment],
                            _controlPoints[segment + 1], localParameter);
}

auto
HermiteCurve::derivative(double parameter, unsigned int order) const -> Point
{
  if (order == 0) {
    return evaluate(parameter);
  }

  if (_controlPoints.size() < 2) {
    return is3D() ? Point{0, 0, 0} : Point{0, 0};
  }

  if (order == 1) {
    auto bounds = parameterBounds();
    if (parameter < bounds.first || parameter > bounds.second) {
      BOOST_THROW_EXCEPTION(Exception("Parameter t is outside valid range"));
    }

    // Find the segment
    size_t numSegments   = _controlPoints.size() - 1;
    double segmentLength = 1.0 / static_cast<double>(numSegments);
    auto   segment       = static_cast<size_t>(parameter / segmentLength);

    if (segment >= numSegments) {
      segment = numSegments - 1;
    }

    // Local parameter within segment
    double localParameter =
        (parameter - (static_cast<double>(segment) * segmentLength)) /
        segmentLength;

    Point tangent = hermiteTangent(_controlPoints[segment],
                                   _controlPoints[segment + 1], localParameter);

    // Scale by inverse segment length
    auto scale = static_cast<double>(numSegments);

    if (is3D()) {
      return {CGAL::to_double(tangent.x()) * scale,
              CGAL::to_double(tangent.y()) * scale,
              CGAL::to_double(tangent.z()) * scale};
    }
    return {CGAL::to_double(tangent.x()) * scale,
            CGAL::to_double(tangent.y()) * scale};
  }

  // Higher order derivatives using finite differences
  const double stepSize         = 1e-8;
  Point        firstDerivative  = derivative(parameter + stepSize, order - 1);
  Point        secondDerivative = derivative(parameter - stepSize, order - 1);

  double deltaX = CGAL::to_double(firstDerivative.x() - secondDerivative.x()) /
                  (2 * stepSize);
  double deltaY = CGAL::to_double(firstDerivative.y() - secondDerivative.y()) /
                  (2 * stepSize);

  if (is3D()) {
    double deltaZ =
        CGAL::to_double(firstDerivative.z() - secondDerivative.z()) /
        (2 * stepSize);
    return {deltaX, deltaY, deltaZ};
  }

  return {deltaX, deltaY};
}

auto
HermiteCurve::toLineString(unsigned int numSegments) const
    -> std::unique_ptr<LineString>
{
  if (_controlPoints.empty()) {
    return std::make_unique<LineString>();
  }

  if (numSegments < 2) {
    BOOST_THROW_EXCEPTION(Exception("Number of segments must be at least 2"));
  }

  auto lineString = std::make_unique<LineString>();

  for (unsigned int segmentIndex = 0; segmentIndex <= numSegments;
       ++segmentIndex) {
    double parameter =
        static_cast<double>(segmentIndex) / static_cast<double>(numSegments);
    lineString->addPoint(evaluate(parameter));
  }

  return lineString;
}

auto
HermiteCurve::parameterBounds() const -> std::pair<double, double>
{
  return std::make_pair(0.0, 1.0);
}

auto
HermiteCurve::isValid() const -> bool
{
  if (_controlPoints.empty()) {
    return true;
  }

  // Check dimensional consistency
  bool first3D = _controlPoints[0].position.is3D();
  bool firstM  = _controlPoints[0].position.isMeasured();

  return std::all_of(_controlPoints.begin(), _controlPoints.end(),
                     [first3D, firstM](const auto &controlPoint) {
                       // Check position dimensions
                       if (controlPoint.position.is3D() != first3D ||
                           controlPoint.position.isMeasured() != firstM) {
                         return false;
                       }

                       // Check handle dimensions
                       return controlPoint.inHandle.is3D() == first3D &&
                              controlPoint.outHandle.is3D() == first3D;
                     });
}

auto
HermiteCurve::toBezierCurve() const -> std::unique_ptr<BezierCurve>
{
  if (_controlPoints.empty()) {
    return std::make_unique<BezierCurve>();
  }

  if (_controlPoints.size() == 1) {
    std::vector<Point> points = {_controlPoints[0].position};
    return std::make_unique<BezierCurve>(points);
  }

  // Convert each segment to cubic Bezier
  std::vector<Point> bezierPoints;

  for (size_t segmentIndex = 0; segmentIndex < _controlPoints.size() - 1;
       ++segmentIndex) {
    const ControlPoint &controlPoint1 = _controlPoints[segmentIndex];
    const ControlPoint &controlPoint2 = _controlPoints[segmentIndex + 1];

    if (segmentIndex == 0) {
      bezierPoints.push_back(controlPoint1.position);
    }

    bezierPoints.push_back(controlPoint1.outHandle);
    bezierPoints.push_back(controlPoint2.inHandle);
    bezierPoints.push_back(controlPoint2.position);
  }

  return std::make_unique<BezierCurve>(bezierPoints);
}

auto
HermiteCurve::fromBezierCurve(const BezierCurve &bezier)
    -> std::unique_ptr<HermiteCurve>
{
  auto hermite = std::make_unique<HermiteCurve>();

  size_t numControlPoints = bezier.numControlPoints();
  if (numControlPoints == 0) {
    return hermite;
  }

  if (numControlPoints == 1) {
    ControlPoint controlPoint(bezier.controlPointAt(0));
    hermite->addControlPoint(controlPoint);
    return hermite;
  }

  // Assume cubic Bezier segments
  if ((numControlPoints - 1) % 3 != 0) {
    BOOST_THROW_EXCEPTION(Exception("Bezier curve must have 3n+1 control "
                                    "points for conversion to Hermite"));
  }

  size_t numSegments = (numControlPoints - 1) / 3;

  for (size_t segmentIndex = 0; segmentIndex <= numSegments; ++segmentIndex) {
    size_t       pointIndex = segmentIndex * 3;
    ControlPoint controlPoint;

    controlPoint.position = bezier.controlPointAt(pointIndex);

    if (segmentIndex > 0) {
      controlPoint.inHandle = bezier.controlPointAt(pointIndex - 1);
    } else {
      controlPoint.inHandle = controlPoint.position;
    }

    if (segmentIndex < numSegments) {
      controlPoint.outHandle = bezier.controlPointAt(pointIndex + 1);
    } else {
      controlPoint.outHandle = controlPoint.position;
    }

    hermite->addControlPoint(controlPoint);
  }

  return hermite;
}

} // namespace SFCGAL
