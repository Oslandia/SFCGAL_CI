// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/Point.h"
#include "SFCGAL/GeometryVisitor.h"

#include "SFCGAL/Exception.h"

using namespace SFCGAL::detail;

namespace SFCGAL {

///
///
///
Point::Point() : _m(NaN()) {}

///
///
///
Point::Point(const Coordinate &coordinate) : _coordinate(coordinate), _m(NaN())
{
}

///
///
///
Point::Point(const Kernel::FT &x, const Kernel::FT &y)
    : _coordinate(x, y), _m(NaN())
{
}

///
///
///
Point::Point(const Kernel::FT &x, const Kernel::FT &y, const Kernel::FT &z,
             const double &m)
    : _coordinate(x, y, z), _m(m)
{
}

///
///
///
Point::Point(const double &x, const double &y) : _coordinate(x, y), _m(NaN()) {}

///
///
///
Point::Point(const double &x, const double &y, const double &z)
    : _coordinate(x, y, z), _m(NaN())
{
}

///
///
///
Point::Point(const double &x, const double &y, const double &z, const double &m)
    : _coordinate(x, y, z), _m(m)
{
}

///
///
///
Point::Point(const Kernel::Point_2 &other) : _coordinate(other), _m(NaN()) {}

///
///
///
Point::Point(const Kernel::Point_3 &other) : _coordinate(other), _m(NaN()) {}

///
///
///
Point::Point(const Point &other)

    = default;

///
///
///
auto
Point::operator=(const Point &other) -> Point &
{
  _coordinate = other._coordinate;
  _m          = other._m;
  return *this;
}

///
///
///
Point::~Point() = default;

///
///
///
auto
Point::clone() const -> Point *
{
  return new Point(*this);
}

///
///
///
auto
Point::geometryType() const -> std::string
{
  return "Point";
}

///
///
///
auto
Point::geometryTypeId() const -> GeometryType
{
  return TYPE_POINT;
}

///
///
///
auto
Point::dimension() const -> int
{
  return 0;
}

///
///
///
auto
Point::coordinateDimension() const -> int
{
  return _coordinate.coordinateDimension() + (isMeasured() ? 1 : 0);
}

///
///
///
auto
Point::isEmpty() const -> bool
{
  return _coordinate.isEmpty();
}

///
///
///
auto
Point::is3D() const -> bool
{
  return _coordinate.is3D();
}

///
///
///
auto
Point::isMeasured() const -> bool
{
  return !std::isnan(_m);
}

///
///
///
void
Point::accept(GeometryVisitor &visitor)
{
  return visitor.visit(*this);
}

///
///
///
void
Point::accept(ConstGeometryVisitor &visitor) const
{
  return visitor.visit(*this);
}

///
///
///
auto
Point::operator<(const Point &other) const -> bool
{
  return _coordinate < other._coordinate;
}

///
///
///
auto
Point::operator==(const Point &other) const -> bool
{
  return _coordinate == other._coordinate;
}

///
///
///
auto
Point::operator!=(const Point &other) const -> bool
{
  return _coordinate != other._coordinate;
}

auto
Point::almostEqual(const Point &other, const double tolerance) const -> bool
{
  return _coordinate.almostEqual(other._coordinate, tolerance);
}

/// @{
/// @privatesection
///
/// Private structures used to implement partial function specialization
template <int D>
struct do_toPoint_d {
  static auto
  toPoint(const Point *p) -> CGAL::Point_2<Kernel>
  {
    return p->toPoint_2();
  }
};

template <>
struct do_toPoint_d<3> {
  static auto
  toPoint(const Point *p) -> CGAL::Point_3<Kernel>
  {
    return p->toPoint_3();
  }
};

/// @} end of private section

template <int Dim>
auto
Point::toPoint_d() const -> typename TypeForDimension<Dim>::Point
{
  return do_toPoint_d<Dim>::toPoint(this);
}
// template instanciations
template CGAL::Point_2<Kernel>
Point::toPoint_d<2>() const;
template CGAL::Point_3<Kernel>
Point::toPoint_d<3>() const;

} // namespace SFCGAL
