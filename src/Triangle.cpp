// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/Triangle.h"
#include "SFCGAL/GeometryVisitor.h"

#include "SFCGAL/Polygon.h"

namespace SFCGAL {

Triangle::Triangle()
{
  _vertices[0] = Point();
  _vertices[1] = Point();
  _vertices[2] = Point();
}

Triangle::Triangle(const Kernel::Triangle_2 &triangle)
{
  for (int i = 0; i < 3; i++) {
    _vertices[i] = triangle.vertex(i);
  }
}

Triangle::Triangle(const Kernel::Triangle_3 &triangle)
{
  for (int i = 0; i < 3; i++) {
    _vertices[i] = triangle.vertex(i);
  }
}

Triangle::Triangle(const Point &p, const Point &q, const Point &r)
{
  _vertices[0] = p;
  _vertices[1] = q;
  _vertices[2] = r;
}

Triangle::Triangle(const Triangle &other) : Surface(other)
{
  _vertices[0] = other._vertices[0];
  _vertices[1] = other._vertices[1];
  _vertices[2] = other._vertices[2];
}

auto
Triangle::operator=(const Triangle &other) -> Triangle &
{
  _vertices[0] = other._vertices[0];
  _vertices[1] = other._vertices[1];
  _vertices[2] = other._vertices[2];
  return *this;
}

Triangle::~Triangle() = default;

auto
Triangle::clone() const -> Triangle *
{
  return new Triangle(*this);
}

auto
Triangle::geometryType() const -> std::string
{
  return "Triangle";
}

auto
Triangle::geometryTypeId() const -> GeometryType
{
  return TYPE_TRIANGLE;
}

auto
Triangle::coordinateDimension() const -> int
{
  return _vertices[0].coordinateDimension();
}

auto
Triangle::isEmpty() const -> bool
{
  return _vertices[0].isEmpty();
}

auto
Triangle::is3D() const -> bool
{
  return _vertices[0].is3D();
}

auto
Triangle::isMeasured() const -> bool
{
  return _vertices[0].isMeasured();
}

auto
Triangle::dropZ() -> bool
{
  if (!is3D()) {
    return false;
  }

  for (auto &_vertex : _vertices) {
    _vertex.dropZ();
  }

  return true;
}

auto
Triangle::dropM() -> bool
{
  if (!isMeasured()) {
    return false;
  }

  for (auto &_vertex : _vertices) {
    _vertex.dropM();
  }

  return true;
}

auto
Triangle::swapXY() -> void
{
  for (auto &_vertex : _vertices) {
    _vertex.swapXY();
  }
}

void
Triangle::reverse()
{
  // note : first point kept to simplify testing
  std::swap(_vertices[1], _vertices[2]);
}

auto
Triangle::toPolygon() const -> Polygon
{
  if (isEmpty()) {
    return {};
  }

  std::vector<Point> points;

  for (size_t i = 0; i < 4; i++) {
    points.push_back(vertex(i));
  }

  return Polygon(LineString(points));
}

void
Triangle::accept(GeometryVisitor &visitor)
{
  return visitor.visit(*this);
}

void
Triangle::accept(ConstGeometryVisitor &visitor) const
{
  return visitor.visit(*this);
}

} // namespace SFCGAL
