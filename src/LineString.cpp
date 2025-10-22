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
    _points.push_back(std::unique_ptr<Point>(point.clone()));
  }
}

LineString::LineString(const Point &startPoint, const Point &endPoint)

{
  _points.push_back(std::unique_ptr<Point>(startPoint.clone()));
  _points.push_back(std::unique_ptr<Point>(endPoint.clone()));
}

LineString::LineString(const LineString &other) : Geometry(other)
{
  _points.reserve(other._points.size());
  for (const auto &point : other._points) {
    _points.push_back(std::unique_ptr<Point>(point->clone()));
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
LineString::clone() const -> LineString *
{
  return new LineString(*this);
}

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
  return visitor.visit(*this);
}

void
LineString::accept(ConstGeometryVisitor &visitor) const
{
  return visitor.visit(*this);
}

auto
LineString::toPolygon_2(bool fixOrientation) const -> CGAL::Polygon_2<Kernel>
{
  if (isEmpty()) {
    return {};
  }

  Point_2_const_iterator pend = points_2_end();
  // skip the last point
  pend--;

  // skip double points
  // TODO: what to do with cycles ?
  std::list<Kernel::Point_2> points;
  Kernel::Point_2            lastP;

  for (Point_2_const_iterator pit = points_2_begin(); pit != pend; ++pit) {
    if (pit == points_2_begin()) {
      lastP = *pit;
      points.push_back(*pit);
      continue;
    }

    if (lastP != *pit) {
      points.push_back(*pit);
    }

    lastP = *pit;
  }

  CGAL::Polygon_2<Kernel> result(points.begin(), points.end());

  if (fixOrientation && result.orientation() == CGAL::CLOCKWISE) {
    result.reverse_orientation();
  }

  return result;
}

} // namespace SFCGAL
