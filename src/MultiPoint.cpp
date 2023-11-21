// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/GeometryVisitor.h>
#include <SFCGAL/MultiPoint.h>

namespace SFCGAL {

///
///
///
MultiPoint::MultiPoint() = default;

///
///
///
MultiPoint::MultiPoint(MultiPoint const &other) = default;

///
///
///
auto
MultiPoint::operator=(MultiPoint other) -> MultiPoint &
{
  swap(other);
  return *this;
}

///
///
///
MultiPoint::~MultiPoint() = default;

///
///
///
auto
MultiPoint::clone() const -> MultiPoint *
{
  return new MultiPoint(*this);
}

///
///
///
auto
MultiPoint::geometryType() const -> std::string
{
  return "MultiPoint";
}

///
///
///
auto
MultiPoint::geometryTypeId() const -> GeometryType
{
  return TYPE_MULTIPOINT;
}

///
///
///
auto
MultiPoint::isAllowed(Geometry const &g) -> bool
{
  return g.geometryTypeId() == TYPE_POINT;
}

///
///
///
void
MultiPoint::accept(GeometryVisitor &visitor)
{
  return visitor.visit(*this);
}

///
///
///
void
MultiPoint::accept(ConstGeometryVisitor &visitor) const
{
  return visitor.visit(*this);
}

} // namespace SFCGAL
