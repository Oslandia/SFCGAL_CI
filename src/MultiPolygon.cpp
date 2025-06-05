// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/GeometryVisitor.h"

namespace SFCGAL {

MultiPolygon::MultiPolygon() = default;

MultiPolygon::MultiPolygon(MultiPolygon const &other)

    = default;

auto
MultiPolygon::operator=(MultiPolygon other) -> MultiPolygon &
{
  swap(other);
  return *this;
}

MultiPolygon::~MultiPolygon() = default;

auto
MultiPolygon::clone() const -> MultiPolygon *
{
  return new MultiPolygon(*this);
}

auto
MultiPolygon::geometryType() const -> std::string
{
  return "MultiPolygon";
}

auto
MultiPolygon::geometryTypeId() const -> GeometryType
{
  return TYPE_MULTIPOLYGON;
}

auto
MultiPolygon::isAllowed(Geometry const &g) -> bool
{
  return g.geometryTypeId() == TYPE_POLYGON;
}

void
MultiPolygon::accept(GeometryVisitor &visitor)
{
  return visitor.visit(*this);
}

void
MultiPolygon::accept(ConstGeometryVisitor &visitor) const
{
  return visitor.visit(*this);
}

} // namespace SFCGAL
