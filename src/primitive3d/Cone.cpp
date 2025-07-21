// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <utility>

#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/primitive3d/Cone.h"

// TODO: how to handle invalid (negative) radiuses?

namespace SFCGAL {

Cone::Cone(Kernel::FT bottom_radius, Kernel::FT top_radius, Kernel::FT height,
           int num_radial)
    : m_bottom_radius(std::move(bottom_radius)),
      m_top_radius(std::move(top_radius)), m_height(std::move(height)),
      m_num_radial(num_radial)
{
}

auto
Cone::operator=(Cone other) -> Cone &
{
  std::swap(m_height, other.m_height);
  std::swap(m_top_radius, other.m_top_radius);
  std::swap(m_bottom_radius, other.m_bottom_radius);
  std::swap(m_polyhedral_surface, other.m_polyhedral_surface);
  std::swap(m_num_radial, other.m_num_radial);
  return *this;
}

void
Cone::setHeight(const Kernel::FT &height)
{
  m_height = height;
  invalidateCache();
}

void
Cone::setBottomRadius(const Kernel::FT &bottom_radius)
{
  m_bottom_radius = bottom_radius;
  invalidateCache();
}

void
Cone::setTopRadius(const Kernel::FT &top_radius)
{
  m_top_radius = top_radius;
  invalidateCache();
}

void
Cone::setNumRadial(const int &num_radial)
{
  m_num_radial = num_radial;
  invalidateCache();
}

void
Cone::invalidateCache()
{
  m_polyhedral_surface.reset();
}

auto
Cone::generatePolyhedralSurface() -> PolyhedralSurface
{
  if (m_polyhedral_surface) {
    return *m_polyhedral_surface;
  }

  // faces
  m_polyhedral_surface = PolyhedralSurface();

  // bottom face
  std::vector<Point> bottom_points;
  if (m_bottom_radius > 0) {
    for (int i = 0; i < m_num_radial; ++i) {
      double     angle = -2.0 * CGAL_PI * i / m_num_radial;
      Kernel::FT x     = m_bottom_radius * std::cos(angle);
      Kernel::FT y     = m_bottom_radius * std::sin(angle);
      bottom_points.emplace_back(x, y, 0);
    }
    bottom_points.push_back(bottom_points[0]);
    m_polyhedral_surface->addPatch(LineString(bottom_points));
  }

  // top face
  std::vector<Point> top_points;
  if (m_top_radius > 0) {
    for (int i = 0; i < m_num_radial; ++i) {
      double     angle = 2.0 * CGAL_PI * i / m_num_radial;
      Kernel::FT x     = m_top_radius * std::cos(angle);
      Kernel::FT y     = m_top_radius * std::sin(angle);
      top_points.emplace_back(x, y, m_height);
    }
    top_points.push_back(top_points[0]);
    m_polyhedral_surface->addPatch(LineString(top_points));
  }

  // side faces if top radius is a single point (true cone)
  if (m_top_radius == 0 && m_bottom_radius > 0) {
    for (int i = 0; i < m_num_radial; ++i) {
      m_polyhedral_surface->addPatch(
          LineString({bottom_points[i], Point(0, 0, m_height),
                      bottom_points[i + 1], bottom_points[i]}));
    }
  }

  // side faces if bottom radius is a single point (true cone upside-down)
  else if (m_bottom_radius == 0 && m_top_radius > 0) {
    for (int i = 0; i < m_num_radial; ++i) {
      m_polyhedral_surface->addPatch(LineString({
          top_points[i],
          Point(0, 0, 0),
          top_points[i + 1],
          top_points[i],
      }));
    }
  }

  // side faces if we have a frustum (part of a cone)
  else if (m_bottom_radius > 0 && m_top_radius > 0) {
    for (int i = 0; i < m_num_radial; ++i) {
      m_polyhedral_surface->addPatch(
          LineString({bottom_points[i], top_points[m_num_radial - i],
                      top_points[m_num_radial - 1 - i], bottom_points[i + 1],
                      bottom_points[i]}));
    }
  }

  return *m_polyhedral_surface;
}

auto
Cone::area(bool withDiscretization) const -> double
{
  double diff_radius =
      std::abs(CGAL::to_double(m_bottom_radius - m_top_radius));
  double slant_height = std::sqrt(std::pow(diff_radius, 2) +
                                  std::pow(CGAL::to_double(m_height), 2));
  if (withDiscretization) {
    double trapezoid_height =
        std::sqrt(std::pow(slant_height, 2) -
                  std::pow(diff_radius * std::sin(CGAL_PI / m_num_radial), 2));
    return (std::sin(CGAL_PI / m_num_radial) *
                CGAL::to_double((m_bottom_radius + m_top_radius)) *
                trapezoid_height +
            std::sin(2.0 * CGAL_PI / m_num_radial) *
                CGAL::to_double((m_bottom_radius * m_bottom_radius +
                                 m_top_radius * m_top_radius)) /
                2.0) *
           m_num_radial;
  }
  return CGAL::to_double(CGAL_PI *
                         ((m_bottom_radius + m_top_radius) * slant_height +
                          m_bottom_radius * m_bottom_radius +
                          m_top_radius * m_top_radius));
}

auto
Cone::volume(bool withDiscretization) const -> double
{
  if (withDiscretization) {
    double bottom_surface = std::pow(CGAL::to_double(m_bottom_radius), 2) *
                            m_num_radial / 2.0 *
                            std::sin(2.0 * CGAL_PI / m_num_radial);
    double top_surface = std::pow(CGAL::to_double(m_top_radius), 2) *
                         m_num_radial / 2.0 *
                         std::sin(2.0 * CGAL_PI / m_num_radial);
    return CGAL::to_double(m_height) / 3.0 *
           (bottom_surface + std::sqrt(bottom_surface * top_surface) +
            top_surface);
  }

  return CGAL::to_double(m_height * CGAL_PI / 3.0 *
                         (m_bottom_radius * m_bottom_radius +
                          m_bottom_radius * m_top_radius +
                          m_top_radius * m_top_radius));
}

auto
Cone::height() const -> const Kernel::FT &
{
  return m_height;
}

auto
Cone::bottomRadius() const -> const Kernel::FT &
{
  return m_bottom_radius;
}

auto
Cone::topRadius() const -> const Kernel::FT &
{
  return m_top_radius;
}

auto
Cone::numRadial() const -> int
{
  return m_num_radial;
}
} // namespace SFCGAL
