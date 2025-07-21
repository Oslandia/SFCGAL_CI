// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/primitive3d/Cube.h"

namespace SFCGAL {

Cube::Cube(const Kernel::FT &size) : m_box(Box(size, size, size)) {}

auto
Cube::operator=(Cube other) -> Cube &
{
  std::swap(m_box, other.m_box);
  return *this;
}

void
Cube::setSize(const Kernel::FT &size)
{
  m_box.setYExtent(size);
  m_box.setXExtent(size);
  m_box.setZExtent(size);
}

auto
Cube::generatePolyhedralSurface() -> PolyhedralSurface
{
  return m_box.generatePolyhedralSurface();
}

auto
Cube::area() const -> double
{
  return m_box.area();
}

auto
Cube::volume() const -> double
{
  return m_box.volume();
}

auto
Cube::size() const -> const Kernel::FT &
{
  return m_box.xExtent();
}

} // namespace SFCGAL
