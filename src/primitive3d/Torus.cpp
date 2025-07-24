// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <utility>

#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/primitive3d/Torus.h"

namespace SFCGAL {

Torus::Torus(const Kernel::FT &main_radius, const Kernel::FT &tube_radius,
             unsigned int main_num_radial, unsigned int tube_num_radial)
{
  if (main_radius <= 0. || tube_radius <= 0.) {
    BOOST_THROW_EXCEPTION(Exception("Torus parameters cannot be negative."));
  } else if (tube_radius >= main_radius) {
    BOOST_THROW_EXCEPTION(
        Exception("Tube radius cannot be greater than main radius."));
  }

  m_parameters["main_radius"]     = main_radius;
  m_parameters["tube_radius"]     = tube_radius;
  m_parameters["main_num_radial"] = main_num_radial;
  m_parameters["tube_num_radial"] = tube_num_radial;
}

auto
Torus::operator=(Torus &other) -> Torus &
{
  Primitive::operator=(other);
  return *this;
}

auto
Torus::primitiveType() const -> std::string
{
  return "Torus";
}

auto
Torus::primitiveTypeId() const -> PrimitiveType
{
  return PrimitiveType::TYPE_TORUS;
}

void
Torus::setMainRadius(const Kernel::FT &main_radius)
{
  if (main_radius <= 0.) {
    BOOST_THROW_EXCEPTION(Exception("Main radius cannot be negative."));
  } else if (tubeRadius() >= main_radius) {
    BOOST_THROW_EXCEPTION(
        Exception("Tube radius cannot be greater than main radius."));
  }

  m_parameters.at("main_radius") = main_radius;
  invalidateCache();
}

void
Torus::setTubeRadius(const Kernel::FT &tube_radius)
{
  if (tube_radius <= 0.) {
    BOOST_THROW_EXCEPTION(Exception("Tube radius cannot be negative."));
  } else if (tube_radius >= mainRadius()) {
    BOOST_THROW_EXCEPTION(
        Exception("Tube radius cannot be greater than main radius."));
  }

  m_parameters.at("tube_radius") = tube_radius;
  invalidateCache();
}

void
Torus::setMainNumRadial(const unsigned int &main_num_radial)
{
  m_parameters.at("main_num_radial") = main_num_radial;
  invalidateCache();
}

void
Torus::setTubeNumRadial(const unsigned int &tube_num_radial)
{
  m_parameters.at("tube_num_radial") = tube_num_radial;
  invalidateCache();
}

auto
Torus::generatePolyhedralSurface() const -> PolyhedralSurface
{
  if (m_polyhedral_surface) {
    return *m_polyhedral_surface;
  }

  m_polyhedral_surface = PolyhedralSurface();
  for (unsigned int i = 0; i < mainNumRadial(); ++i) {
    for (unsigned int j = 0; j < tubeNumRadial(); ++j) {
      std::vector<Point> points;
      for (unsigned int k = 0; k < 4; ++k) {
        double main_angle =
            2.0 * CGAL_PI * ((i + (k & 1)) % mainNumRadial()) / mainNumRadial();
        double tube_angle = 2.0 * CGAL_PI * ((j + (k >> 1)) % tubeNumRadial()) /
                            tubeNumRadial();

        double x = (CGAL::to_double(mainRadius()) +
                    CGAL::to_double(tubeRadius()) * std::cos(tube_angle)) *
                   std::cos(main_angle);
        double y = (CGAL::to_double(mainRadius()) +
                    CGAL::to_double(tubeRadius()) * std::cos(tube_angle)) *
                   std::sin(main_angle);
        double z = CGAL::to_double(tubeRadius()) * std::sin(tube_angle);
        points.emplace_back(x, y, z);
      }
      m_polyhedral_surface->addPatch(
          LineString({points[0], points[1], points[3], points[2], points[0]}));
    }
  }

  return *m_polyhedral_surface;
}

auto
Torus::area3D() const -> double
{
  return 4.0 * std::pow(CGAL_PI, 2) *
         CGAL::to_double(mainRadius() * tubeRadius());
}

auto
Torus::volume() const -> double
{
  return 2.0 * std::pow(CGAL_PI, 2) *
         CGAL::to_double(mainRadius() * tubeRadius() * tubeRadius());
}

auto
Torus::mainRadius() const -> const Kernel::FT &
{
  return std::get<Kernel::FT>(m_parameters.at("main_radius"));
}
auto

Torus::tubeRadius() const -> const Kernel::FT &
{
  return std::get<Kernel::FT>(m_parameters.at("tube_radius"));
}

auto
Torus::mainNumRadial() const -> unsigned int
{
  return std::get<unsigned int>(m_parameters.at("main_num_radial"));
}

auto
Torus::tubeNumRadial() const -> unsigned int
{
  return std::get<unsigned int>(m_parameters.at("tube_num_radial"));
}

} // namespace SFCGAL
