// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/GeometryVisitor.h>
#include <SFCGAL/MultiSolid.h>

namespace SFCGAL {

///
///
///
MultiSolid::MultiSolid()  {}

///
///
///
MultiSolid::MultiSolid(MultiSolid const &other) = default;

///
///
///
auto
MultiSolid::operator=(MultiSolid other) -> MultiSolid &
{
  swap(other);
  return *this;
}

///
///
///
MultiSolid::~MultiSolid() = default;

///
///
///
auto
MultiSolid::clone() const -> MultiSolid *
{
  return new MultiSolid(*this);
}

///
///
///
auto
MultiSolid::geometryType() const -> std::string
{
  return "MultiSolid";
}

///
///
///
auto
MultiSolid::geometryTypeId() const -> GeometryType
{
  return TYPE_MULTISOLID;
}

///
///
///
auto
MultiSolid::isAllowed(Geometry const &g) -> bool
{
  return g.geometryTypeId() == TYPE_SOLID;
}

///
///
///
void
MultiSolid::accept(GeometryVisitor &visitor)
{
  return visitor.visit(*this);
}

///
///
///
void
MultiSolid::accept(ConstGeometryVisitor &visitor) const
{
  return visitor.visit(*this);
}

} // namespace SFCGAL
