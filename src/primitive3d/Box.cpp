// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <utility>

#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/primitive3d/Box.h"

namespace SFCGAL {

Box::Box(Kernel::FT x_extent, Kernel::FT y_extent, Kernel::FT z_extent)
    : m_x_extent(std::move(x_extent)), m_y_extent(std::move(y_extent)),
      m_z_extent(std::move(z_extent))
{
}

auto
Box::operator=(Box other) -> Box &
{
  std::swap(m_x_extent, other.m_x_extent);
  std::swap(m_y_extent, other.m_y_extent);
  std::swap(m_z_extent, other.m_z_extent);
  std::swap(m_polyhedral_surface, other.m_polyhedral_surface);
  return *this;
}

void
Box::setXExtent(const Kernel::FT &x_extent)
{
  m_x_extent = x_extent;
  invalidateCache();
}

void
Box::setYExtent(const Kernel::FT &y_extent)
{
  m_y_extent = y_extent;
  invalidateCache();
}

void
Box::setZExtent(const Kernel::FT &z_extent)
{
  m_z_extent = z_extent;
  invalidateCache();
}

auto
Box::xExtent() const -> const Kernel::FT &
{
  return m_x_extent;
}

auto
Box::yExtent() const -> const Kernel::FT &
{
  return m_y_extent;
}

auto
Box::zExtent() const -> const Kernel::FT &
{
  return m_z_extent;
}

void
Box::invalidateCache()
{
  m_polyhedral_surface.reset();
}

auto
Box::generatePolyhedralSurface() -> PolyhedralSurface
{
  if (m_polyhedral_surface) {
    return *m_polyhedral_surface;
  }

  m_polyhedral_surface = PolyhedralSurface();

  // vertices
  std::vector<Point> points{{0, 0, 0},
                            {m_x_extent, 0, 0},
                            {0, m_y_extent, 0},
                            {m_x_extent, m_y_extent, 0},
                            {0, 0, m_z_extent},
                            {m_x_extent, 0, m_z_extent},
                            {0, m_y_extent, m_z_extent},
                            {m_x_extent, m_y_extent, m_z_extent}};

  // faces
  m_polyhedral_surface->addPatch(
      LineString({points[0], points[2], points[3], points[1], points[0]}));
  m_polyhedral_surface->addPatch(
      LineString({points[4], points[5], points[7], points[6], points[4]}));
  m_polyhedral_surface->addPatch(
      LineString({points[0], points[1], points[5], points[4], points[0]}));
  m_polyhedral_surface->addPatch(
      LineString({points[2], points[6], points[7], points[3], points[2]}));
  m_polyhedral_surface->addPatch(
      LineString({points[1], points[3], points[7], points[5], points[1]}));
  m_polyhedral_surface->addPatch(
      LineString({points[0], points[4], points[6], points[2], points[0]}));

  // create polyhedral surface
  return *m_polyhedral_surface;
}

auto
Box::area() const -> double
{
  return CGAL::to_double(2 *
                         (m_x_extent * m_y_extent + m_y_extent * m_z_extent +
                          m_z_extent * m_x_extent));
}

auto
Box::volume() const -> double
{
  return CGAL::to_double(m_x_extent * m_y_extent * m_z_extent);
}

} // namespace SFCGAL
