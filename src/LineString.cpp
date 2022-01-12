// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/GeometryVisitor.h>
#include <SFCGAL/LineString.h>

namespace SFCGAL {

///
///
///
LineString::LineString() : Geometry(), _points() {}

///
///
///
LineString::LineString(const std::vector<Point> &points) : Geometry(), _points()
{
  for (size_t i = 0; i < points.size(); i++) {
    _points.push_back(points[i].clone());
  }
}

///
///
///
LineString::LineString(const Point &startPoint, const Point &endPoint)
    : Geometry(), _points()
{
  _points.push_back(startPoint.clone());
  _points.push_back(endPoint.clone());
}

///
///
///
LineString::LineString(const LineString &other) : Geometry(other)
{
  for (size_t i = 0; i < other.numPoints(); i++) {
    _points.push_back(other.pointN(i).clone());
  }
}

///
///
///
LineString &
LineString::operator=(LineString other)
{
  swap(other);
  return *this;
}

///
///
///
LineString::~LineString() {}

///
///
///
LineString *
LineString::clone() const
{
  return new LineString(*this);
}

///
///
///
GeometryType
LineString::geometryTypeId() const
{
  return TYPE_LINESTRING;
}

///
///
///
std::string
LineString::geometryType() const
{
  return "LineString";
}

///
///
///
int
LineString::dimension() const
{
  return 1;
}

///
///
int
LineString::coordinateDimension() const
{
  return isEmpty() ? 0 : _points[0].coordinateDimension();
}

///
///
///
bool
LineString::isEmpty() const
{
  return _points.empty();
}

///
///
///
bool
LineString::is3D() const
{
  return !isEmpty() && startPoint().is3D();
}

///
///
///
bool
LineString::isMeasured() const
{
  return !isEmpty() && startPoint().isMeasured();
}

///
///
///
void
LineString::clear()
{
  _points.clear();
}

///
///
///
void
LineString::reverse()
{
  std::reverse(_points.begin(), _points.end());
}

///
///
///
size_t
LineString::numSegments() const
{
  if (_points.empty()) {
    return 0;
  } else {
    return _points.size() - 1;
  }
}

///
///
///
bool
LineString::isClosed() const
{
  return (!isEmpty()) && (startPoint() == endPoint());
}

///
///
///
void
LineString::reserve(const size_t &n)
{
  _points.reserve(n);
}

///
///
///
void
LineString::accept(GeometryVisitor &visitor)
{
  return visitor.visit(*this);
}

///
///
///
void
LineString::accept(ConstGeometryVisitor &visitor) const
{
  return visitor.visit(*this);
}

///
///
///
CGAL::Polygon_2<Kernel>
LineString::toPolygon_2(bool fixOrientation) const
{
  if (isEmpty()) {
    return CGAL::Polygon_2<Kernel>();
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
