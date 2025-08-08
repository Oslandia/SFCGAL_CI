// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <utility>

#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/primitive3d/Torus.h"

// TODO: how to handle invalid radiuses?
//         - negative
//         - tube_radius >= main_radius

namespace SFCGAL {

Torus::Torus(Kernel::FT main_radius, Kernel::FT tube_radius,
             int main_num_radial, int tube_num_radial)
    : m_main_radius(std::move(main_radius)),
      m_tube_radius(std::move(tube_radius)), m_main_num_radial(main_num_radial),
      m_tube_num_radial(tube_num_radial)
{
}

auto
Torus::operator=(Torus other) -> Torus &
{
  std::swap(m_main_radius, other.m_main_radius);
  std::swap(m_tube_radius, other.m_tube_radius);
  std::swap(m_polyhedral_surface, other.m_polyhedral_surface);
  std::swap(m_main_num_radial, other.m_main_num_radial);
  std::swap(m_tube_num_radial, other.m_tube_num_radial);
  return *this;
}

void
Torus::setMainRadius(const Kernel::FT &main_radius)
{
  m_main_radius = main_radius;
  invalidateCache();
}

void
Torus::setTubeRadius(const Kernel::FT &tube_radius)
{
  m_tube_radius = tube_radius;
  invalidateCache();
}

void
Torus::setMainNumRadial(const int &main_num_radial)
{
  m_main_num_radial = main_num_radial;
  invalidateCache();
}

void
Torus::setTubeNumRadial(const int &tube_num_radial)
{
  m_tube_num_radial = tube_num_radial;
  invalidateCache();
}

void
Torus::invalidateCache()
{
  m_polyhedral_surface.reset();
}

auto
Torus::generatePolyhedralSurface() const -> PolyhedralSurface
{
  if (m_polyhedral_surface) {
    return *m_polyhedral_surface;
  }

  // FIXME: time and memory consuming starting from 999 m_main_num_radial

  m_polyhedral_surface = PolyhedralSurface();
  for (int i = 0; i < m_main_num_radial; ++i) {
    for (int j = 0; j < m_tube_num_radial; ++j) {
      std::vector<Point> points;
      for (int k = 0; k < 4; ++k) {
        double main_angle = 2.0 * CGAL_PI *
                            ((i + (k & 1)) % m_main_num_radial) /
                            m_main_num_radial;
        double tube_angle = 2.0 * CGAL_PI *
                            ((j + (k >> 1)) % m_tube_num_radial) /
                            m_tube_num_radial;

        double x = (CGAL::to_double(m_main_radius) +
                    CGAL::to_double(m_tube_radius) * std::cos(tube_angle)) *
                   std::cos(main_angle);
        double y = (CGAL::to_double(m_main_radius) +
                    CGAL::to_double(m_tube_radius) * std::cos(tube_angle)) *
                   std::sin(main_angle);
        double z = CGAL::to_double(m_tube_radius) * std::sin(tube_angle);
        points.emplace_back(x, y, z);
      }
      m_polyhedral_surface->addPatch(
          LineString({points[0], points[1], points[3], points[2], points[0]}));
    }
  }

  return *m_polyhedral_surface;
}

auto
Torus::area() const -> double
{
  return 4.0 * std::pow(CGAL_PI, 2) *
         CGAL::to_double(m_main_radius * m_tube_radius);
}

auto
Torus::volume() const -> double
{
  return 2.0 * std::pow(CGAL_PI, 2) *
         CGAL::to_double(m_main_radius * m_tube_radius * m_tube_radius);
}

auto
Torus::mainRadius() const -> const Kernel::FT &
{
  return m_main_radius;
}
auto

Torus::tubeRadius() const -> const Kernel::FT &
{
  return m_tube_radius;
}

auto
Torus::mainNumRadial() const -> int
{
  return m_main_num_radial;
}

auto
Torus::tubeNumRadial() const -> int
{
  return m_tube_num_radial;
}

} // namespace SFCGAL
