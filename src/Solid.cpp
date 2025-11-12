// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/Solid.h"
#include "SFCGAL/GeometryVisitor.h"
#include "SFCGAL/PolyhedralSurface.h"
#include <memory>

namespace SFCGAL {

Solid::Solid() { _shells.push_back(std::make_unique<PolyhedralSurface>()); }

Solid::Solid(const PolyhedralSurface &exteriorShell)
{
  _shells.push_back(exteriorShell.clone());
}

Solid::Solid(PolyhedralSurface *exteriorShell)
{
  _shells.emplace_back(exteriorShell);
}

Solid::Solid(const std::vector<PolyhedralSurface> &shells)
{
  if (shells.empty()) {
    _shells.emplace_back(std::make_unique<PolyhedralSurface>());
  } else {
    _shells.reserve(shells.size());
    for (const auto &shell : shells) {
      _shells.emplace_back(shell.clone());
    }
  }
}

Solid::Solid(const Solid &other) : GeometryImpl(other)
{
  _shells.reserve(other._shells.size());
  for (const auto &shell : other._shells) {
    _shells.emplace_back(shell->clone());
  }
}

auto
Solid::operator=(Solid other) -> Solid &
{
  swap(other);
  return *this;
}

Solid::~Solid() = default;

auto
Solid::geometryType() const -> std::string
{
  return "Solid";
}

auto
Solid::geometryTypeId() const -> GeometryType
{
  return TYPE_SOLID;
}

auto
Solid::dimension() const -> int
{
  return 3;
}

auto
Solid::coordinateDimension() const -> int
{
  return exteriorShell().coordinateDimension();
}

auto
Solid::isEmpty() const -> bool
{
  return exteriorShell().isEmpty();
}

auto
Solid::is3D() const -> bool
{
  return exteriorShell().is3D();
}

auto
Solid::isMeasured() const -> bool
{
  return exteriorShell().isMeasured();
}

auto
Solid::dropZ() -> bool
{
  if (!is3D()) {
    return false;
  }

  for (auto &shell : _shells) {
    shell->dropZ();
  }

  return true;
}

auto
Solid::dropM() -> bool
{
  if (!isMeasured()) {
    return false;
  }

  for (auto &shell : _shells) {
    shell->dropM();
  }

  return true;
}

auto
Solid::swapXY() -> void
{
  for (auto &shell : _shells) {
    shell->swapXY();
  }
}

void
Solid::accept(GeometryVisitor &visitor)
{
  visitor.visit(*this);
}

void
Solid::accept(ConstGeometryVisitor &visitor) const
{
  visitor.visit(*this);
}

auto
Solid::toSurfaceMesh() const -> Surface_mesh_3
{
  return exteriorShell().toSurfaceMesh();
}

} // namespace SFCGAL
