// Copyright (c) 2025-2025, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/Segment.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/algorithm/distance.h"
#include "SFCGAL/algorithm/distance3d.h"
#include <CGAL/intersections.h>
#include <algorithm>

namespace SFCGAL {

Segment::Segment() : _source(), _target() {}

// Private helper method for validation
static auto
validatePointDimensions(const Point &point1, const Point &point2,
                        const char *point1Name, const char *point2Name) -> void
{
  // Check spatial dimension consistency (2D/3D)
  if (point1.is3D() != point2.is3D()) {
    std::string message =
        "Segment endpoints must have the same spatial dimension. ";
    message += point1.is3D() ? std::string(point1Name) + " is 3D, "
                             : std::string(point1Name) + " is 2D, ";
    message += point2.is3D() ? "but " + std::string(point2Name) + " is 3D."
                             : "but " + std::string(point2Name) + " is 2D.";
    throw Exception(message);
  }

  // Check measurement dimension consistency (M value)
  if (point1.isMeasured() != point2.isMeasured()) {
    std::string message =
        "Segment endpoints must both be measured or both unmeasured. ";
    message += point1.isMeasured()
                   ? std::string(point1Name) + " is measured, "
                   : std::string(point1Name) + " is not measured, ";
    message += point2.isMeasured()
                   ? "but " + std::string(point2Name) + " is measured."
                   : "but " + std::string(point2Name) + " is not measured.";
    throw Exception(message);
  }
}

Segment::Segment(const Point &p1, const Point &p2)
{
  // Check for emptiness first
  if (p1.isEmpty() || p2.isEmpty()) {
    _source = Point(); // Empty point
    _target = Point(); // Empty point
    return;
  }

  // Use the helper method for dimension validation
  validatePointDimensions(p1, p2, "First point", "second point");

  _source = p1;
  _target = p2;
}

auto
Segment::setSource(const Point &p) -> void
{
  // If segment is empty, refuse to set just one point
  if (isEmpty()) {
    throw Exception("Cannot set source point on an empty segment. Use "
                    "setPoints to define both endpoints.");
  }

  validatePointDimensions(p, _target, "New point", "target");
  _source = p;
}

auto
Segment::setTarget(const Point &p) -> void
{
  // If segment is empty, refuse to set just one point
  if (isEmpty()) {
    throw Exception("Cannot set target point on an empty segment. Use "
                    "setPoints to define both endpoints.");
  }

  validatePointDimensions(p, _source, "New point", "source");
  _target = p;
}

auto
Segment::setPoints(const Point &source, const Point &target) -> void
{
  // Use the constructor logic by creating a temporary segment
  Segment temp(source, target);

  // If we got here, no exception was thrown, so the points are compatible
  _source = source;
  _target = target;
}

auto
Segment::clear() -> void
{
  _source = Point();
  _target = Point();
}

auto
Segment::toSegment_2() const -> Kernel::Segment_2
{
  return Kernel::Segment_2(_source.toPoint_2(), _target.toPoint_2());
}

auto
Segment::toSegment_3() const -> Kernel::Segment_3
{
  return Kernel::Segment_3(_source.toPoint_3(), _target.toPoint_3());
}

auto
Segment::isEmpty() const -> bool
{
  return _source.isEmpty() || _target.isEmpty();
}

auto
Segment::hasSameDimension() const -> bool
{
  return _source.is3D() == _target.is3D() &&
         _source.isMeasured() == _target.isMeasured();
}

auto
Segment::isDegenerate() const -> bool
{
  if (isEmpty()) {
    return true;
  }

  return applyByDimension([this]() { return toSegment_2().is_degenerate(); },
                          [this]() { return toSegment_3().is_degenerate(); });
}

auto
Segment::squaredLength() const -> Kernel::FT
{
  if (isEmpty()) {
    return 0;
  }

  return applyByDimension([this]() { return toSegment_2().squared_length(); },
                          [this]() { return toSegment_3().squared_length(); });
}

auto
Segment::length() const -> Kernel::FT
{
  if (isEmpty()) {
    return 0;
  }

  return CGAL::sqrt(CGAL::to_double(squaredLength()));
}

auto
Segment::is3D() const -> bool
{
  // Segment is 3D if any endpoint is 3D
  // Unlike to OGC, Segment, here, can not be EMPTY Z
  return !isEmpty() && _source.is3D();
}

auto
Segment::isMeasured() const -> bool
{
  // Segment is measured if any endpoint is measured
  // Unlike to OGC, Segment, here, can not be EMPTY M
  return !isEmpty() && _source.isMeasured();
}

auto
Segment::interpolate(double t) const -> Point
{
  if (isEmpty()) {
    return Point();
  }

  // Clamp parameter to [0,1]
  t = std::clamp(t, 0.0, 1.0);

  // Create base point with interpolated coordinates
  Point result = applyByDimension(
      [this, t]() {
        const auto p1 = _source.toPoint_2();
        const auto p2 = _target.toPoint_2();
        return Point(p1 + t * (p2 - p1));
      },
      [this, t]() {
        const auto p1 = _source.toPoint_3();
        const auto p2 = _target.toPoint_3();
        return Point(p1 + t * (p2 - p1));
      });

  // Handle M value interpolation if segment is measured
  if (isMeasured()) {
    double m1 = _source.m();
    double m2 = _target.m();
    result.setM(m1 + t * (m2 - m1));
  }

  return result;
}

auto
Segment::midpoint() const -> Point
{
  // The midpoint is simply the interpolation at t=0.5
  return interpolate(0.5);
}

auto
Segment::reverse() -> void
{
  if (!isEmpty()) {
    std::swap(_source, _target);
  }
}

} // namespace SFCGAL
