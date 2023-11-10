// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/GeometryVisitor.h>
#include <SFCGAL/MultiPolygon.h>

namespace SFCGAL {

///
///
///
MultiPolygon::MultiPolygon()  = default;

///
///
///
MultiPolygon::MultiPolygon(MultiPolygon const &other)

    = default;

///
///
///
auto
MultiPolygon::operator=(MultiPolygon other) -> MultiPolygon &
{
  swap(other);
  return *this;
}

///
///
///
MultiPolygon::~MultiPolygon() = default;

///
///
///
auto
MultiPolygon::clone() const -> MultiPolygon *
{
  return new MultiPolygon(*this);
}

///
///
///
auto
MultiPolygon::geometryType() const -> std::string
{
  return "MultiPolygon";
}

///
///
///
auto
MultiPolygon::geometryTypeId() const -> GeometryType
{
  return TYPE_MULTIPOLYGON;
}

///
///
///
auto
MultiPolygon::isAllowed(Geometry const &g) -> bool
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
