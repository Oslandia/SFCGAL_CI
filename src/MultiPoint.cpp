// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/GeometryVisitor.h>
#include <SFCGAL/MultiPoint.h>

namespace SFCGAL {

///
///
///
MultiPoint::MultiPoint() : GeometryCollection() {}

///
///
///
MultiPoint::MultiPoint(MultiPoint const &other) : GeometryCollection(other) {}

///
///
///
MultiPoint &
MultiPoint::operator=(MultiPoint other)
{
  swap(other);
  return *this;
}

///
///
///
MultiPoint::~MultiPoint() {}

///
///
///
MultiPoint *
MultiPoint::clone() const
{
  return new MultiPoint(*this);
}

///
///
///
std::string
MultiPoint::geometryType() const
{
  return "MultiPoint";
}

///
///
///
GeometryType
MultiPoint::geometryTypeId() const
{
  return TYPE_MULTIPOINT;
}

///
///
///
bool
MultiPoint::isAllowed(Geometry const &g)
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
