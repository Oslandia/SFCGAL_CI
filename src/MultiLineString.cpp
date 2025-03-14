// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/GeometryVisitor.h"

namespace SFCGAL {

MultiLineString::MultiLineString() = default;

MultiLineString::MultiLineString(MultiLineString const &other)

    = default;

auto
MultiLineString::operator=(MultiLineString other) -> MultiLineString &
{
  swap(other);
  return *this;
}

MultiLineString::~MultiLineString() = default;

auto
MultiLineString::clone() const -> MultiLineString *
{
  return new MultiLineString(*this);
}

auto
MultiLineString::geometryType() const -> std::string
{
  return "MultiLineString";
}

auto
MultiLineString::geometryTypeId() const -> GeometryType
{
  return TYPE_MULTILINESTRING;
}

auto
MultiLineString::isAllowed(Geometry const &g) -> bool
{
  return g.geometryTypeId() == TYPE_LINESTRING;
}

void
MultiLineString::accept(GeometryVisitor &visitor)
{
  return visitor.visit(*this);
}

void
MultiLineString::accept(ConstGeometryVisitor &visitor) const
{
  return visitor.visit(*this);
}

} // namespace SFCGAL
