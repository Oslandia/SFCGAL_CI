// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/LineString.h"
#include "SFCGAL/GeometryVisitor.h"

namespace SFCGAL {

LineString::LineString() = default;

LineString::LineString(const std::vector<Point> &points)
{
  for (const auto &point : points) {
    _points.push_back(point.clone());
  }
}

LineString::LineString(const Point &startPoint, const Point &endPoint)

{
  _points.push_back(startPoint.clone());
  _points.push_back(endPoint.clone());
}

LineString::LineString(const LineString &other) : GeometryImpl(other)
{
  _points.reserve(other._points.size());
  for (const auto &point : other._points) {
    _points.push_back(point->clone());
  }
}

auto
LineString::operator=(LineString other) -> LineString &
{
  swap(other);
  return *this;
}

LineString::~LineString() = default;

auto
LineString::geometryTypeId() const -> GeometryType
{
  return TYPE_LINESTRING;
}

auto
LineString::geometryType() const -> std::string
{
  return "LineString";
}

auto
LineString::dimension() const -> int
{
  return 1;
}

auto
LineString::coordinateDimension() const -> int
{
  return isEmpty() ? 0 : _points[0]->coordinateDimension();
}

auto
LineString::isEmpty() const -> bool
{
  return _points.empty();
}

auto
LineString::is3D() const -> bool
{
  return !isEmpty() && startPoint().is3D();
}

auto
LineString::isMeasured() const -> bool
{
  return !isEmpty() && startPoint().isMeasured();
}

auto
LineString::dropZ() -> bool
{
  if (!is3D()) {
    return false;
  }

  for (auto &_point : _points) {
    _point->dropZ();
  }

  return true;
}

auto
LineString::dropM() -> bool
{
  if (!isMeasured()) {
    return false;
  }

  for (auto &_point : _points) {
    _point->dropM();
  }

  return true;
}

auto
LineString::swapXY() -> void
{
  for (auto &_point : _points) {
    _point->swapXY();
  }
}

void
LineString::clear()
{
  _points.clear();
}

void
LineString::reverse()
{
  std::reverse(_points.begin(), _points.end());
}

auto
LineString::numSegments() const -> size_t
{
  if (_points.empty()) {
    return 0;
  }
  return _points.size() - 1;
}

auto
LineString::isClosed() const -> bool
{
  return (!isEmpty()) && (startPoint() == endPoint());
}

void
LineString::reserve(const size_t &n)
{
  _points.reserve(n);
}

void
LineString::accept(GeometryVisitor &visitor)
{
  visitor.visit(*this);
}

void
LineString::accept(ConstGeometryVisitor &visitor) const
{
  visitor.visit(*this);
}

auto
LineString::toPolygon_2(bool fixOrientation) const -> CGAL::Polygon_2<Kernel>
{
  if (isEmpty()) {
    return {};
  }

  Point_2_const_iterator pointEndIterator = points_2_end();
  // skip the last point
  pointEndIterator--;

  // skip double points
  // TODO: what to do with cycles ?
  std::list<Kernel::Point_2> points;
  Kernel::Point_2            lastPoint;

  for (Point_2_const_iterator pointIterator = points_2_begin();
       pointIterator != pointEndIterator; ++pointIterator) {
    if (pointIterator == points_2_begin()) {
      lastPoint = *pointIterator;
      points.push_back(*pointIterator);
      continue;
    }

    if (lastPoint != *pointIterator) {
      points.push_back(*pointIterator);
    }

    lastPoint = *pointIterator;
  }

  CGAL::Polygon_2<Kernel> result(points.begin(), points.end());

  if (fixOrientation && result.orientation() == CGAL::CLOCKWISE) {
    result.reverse_orientation();
  }

  return result;
}

} // namespace SFCGAL
