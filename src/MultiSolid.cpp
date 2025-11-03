// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/GeometryVisitor.h"

namespace SFCGAL {

MultiSolid::MultiSolid() = default;

MultiSolid::MultiSolid(MultiSolid const &other) = default;

auto
MultiSolid::operator=(MultiSolid other) -> MultiSolid &
{
  swap(other);
  return *this;
}

MultiSolid::~MultiSolid() = default;

auto
MultiSolid::geometryType() const -> std::string
{
  return "MultiSolid";
}

auto
MultiSolid::geometryTypeId() const -> GeometryType
{
  return TYPE_MULTISOLID;
}

auto
MultiSolid::isAllowed(Geometry const &geometry) -> bool
{
  return geometry.geometryTypeId() == TYPE_SOLID;
}

auto
MultiSolid::accept(GeometryVisitor &visitor) -> void
{
  return visitor.visit(*this);
}

auto
MultiSolid::accept(ConstGeometryVisitor &visitor) const -> void
{
  return visitor.visit(*this);
}

} // namespace SFCGAL
