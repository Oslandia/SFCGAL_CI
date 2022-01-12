// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/GeometryVisitor.h>
#include <SFCGAL/MultiPolygon.h>

namespace SFCGAL {

///
///
///
MultiPolygon::MultiPolygon() : GeometryCollection() {}

///
///
///
MultiPolygon::MultiPolygon(MultiPolygon const &other)
    : GeometryCollection(other)
{
}

///
///
///
MultiPolygon &
MultiPolygon::operator=(MultiPolygon other)
{
  swap(other);
  return *this;
}

///
///
///
MultiPolygon::~MultiPolygon() {}

///
///
///
MultiPolygon *
MultiPolygon::clone() const
{
  return new MultiPolygon(*this);
}

///
///
///
std::string
MultiPolygon::geometryType() const
{
  return "MultiPolygon";
}

///
///
///
GeometryType
MultiPolygon::geometryTypeId() const
{
  return TYPE_MULTIPOLYGON;
}

///
///
///
bool
MultiPolygon::isAllowed(Geometry const &g)
{
  return g.geometryTypeId() == TYPE_POLYGON;
}

///
///
///
void
MultiPolygon::accept(GeometryVisitor &visitor)
{
  return visitor.visit(*this);
}

///
///
///
void
MultiPolygon::accept(ConstGeometryVisitor &visitor) const
{
  return visitor.visit(*this);
}

} // namespace SFCGAL
