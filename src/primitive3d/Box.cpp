// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <utility>

#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/primitive3d/Box.h"

namespace SFCGAL {

Box::Box(const Kernel::FT &x_extent, const Kernel::FT &y_extent,
         const Kernel::FT &z_extent)
{
  m_parameters["x_extent"] = x_extent;
  m_parameters["y_extent"] = y_extent;
  m_parameters["z_extent"] = z_extent;

  Box::validateParameters(m_parameters);
}

auto
Box::operator=(Box &other) -> Box &
{
  if (this != &other) {
    Primitive::operator=(other);
  }
  return *this;
}

auto
Box::primitiveType() const -> std::string
{
  return "Box";
}

auto
Box::primitiveTypeId() const -> PrimitiveType
{
  return PrimitiveType::TYPE_BOX;
}

void
Box::setXExtent(const Kernel::FT &x_extent)
{
  validateAndSetParameter("x_extent", x_extent);
}

void
Box::setYExtent(const Kernel::FT &y_extent)
{
  validateAndSetParameter("y_extent", y_extent);
}

void
Box::setZExtent(const Kernel::FT &z_extent)
{
  validateAndSetParameter("z_extent", z_extent);
}

void
Box::validateParameters(
    std::unordered_map<std::string, PrimitiveParameter> const &tempParameters)
    const
{
  const double xExtent =
      CGAL::to_double(std::get<Kernel::FT>(tempParameters.at("x_extent")));
  const double yExtent =
      CGAL::to_double(std::get<Kernel::FT>(tempParameters.at("y_extent")));
  const double zExtent =
      CGAL::to_double(std::get<Kernel::FT>(tempParameters.at("z_extent")));

  if (xExtent < 0. || yExtent < 0. || zExtent < 0) {
    BOOST_THROW_EXCEPTION(Exception("Box parameters cannot be negative."));
  }
}

auto
Box::xExtent() const -> const Kernel::FT &
{
  return std::get<Kernel::FT>(m_parameters.at("x_extent"));
}

auto
Box::yExtent() const -> const Kernel::FT &
{
  return std::get<Kernel::FT>(m_parameters.at("y_extent"));
}

auto
Box::zExtent() const -> const Kernel::FT &
{
  return std::get<Kernel::FT>(m_parameters.at("z_extent"));
}

auto
Box::generatePolyhedralSurface() const -> PolyhedralSurface
{
  if (m_polyhedral_surface) {
    return *m_polyhedral_surface;
  }

  m_polyhedral_surface = PolyhedralSurface();

  // vertices
  std::vector<Point> points{{0, 0, 0},
                            {xExtent(), 0, 0},
                            {0, yExtent(), 0},
                            {xExtent(), yExtent(), 0},
                            {0, 0, zExtent()},
                            {xExtent(), 0, zExtent()},
                            {0, yExtent(), zExtent()},
                            {xExtent(), yExtent(), zExtent()}};

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
Box::area3D() const -> double
{
  return CGAL::to_double(2 * (xExtent() * yExtent() + yExtent() * zExtent() +
                              zExtent() * xExtent()));
}

auto
Box::volume() const -> double
{
  return CGAL::to_double(xExtent() * yExtent() * zExtent());
}

} // namespace SFCGAL
