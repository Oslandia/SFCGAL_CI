// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <utility>

#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/primitive3d/Cone.h"
#include "SFCGAL/primitive3d/Primitive.h"

namespace SFCGAL {

Cone::Cone(const Kernel::FT &bottom_radius, const Kernel::FT &top_radius,
           const Kernel::FT &height, unsigned int num_radial)
{
  m_parameters["bottom_radius"] = bottom_radius;
  m_parameters["top_radius"]    = top_radius;
  m_parameters["height"]        = height;
  m_parameters["num_radial"]    = num_radial;

  Cone::validateParameters(m_parameters);
}

auto
Cone::primitiveType() const -> std::string
{
  return "Cone";
}

auto
Cone::primitiveTypeId() const -> PrimitiveType
{
  return PrimitiveType::TYPE_CONE;
}

void
Cone::setHeight(const Kernel::FT &height)
{
  validateAndSetParameter("height", height);
}

void
Cone::setBottomRadius(const Kernel::FT &bottom_radius)
{
  validateAndSetParameter("bottom_radius", bottom_radius);
}

void
Cone::setTopRadius(const Kernel::FT &top_radius)
{
  validateAndSetParameter("top_radius", top_radius);
}

void
Cone::setNumRadial(const unsigned int &num_radial)
{
  validateAndSetParameter("num_radial", num_radial);
}

void
Cone::validateParameters(
    std::unordered_map<std::string, PrimitiveParameter> const &tempParameters)
    const
{
  const double height =
      CGAL::to_double(std::get<Kernel::FT>(tempParameters.at("height")));
  const double bottomRadius =
      CGAL::to_double(std::get<Kernel::FT>(tempParameters.at("bottom_radius")));
  const double topRadius =
      CGAL::to_double(std::get<Kernel::FT>(tempParameters.at("top_radius")));

  if (height <= 0. || bottomRadius < 0. || topRadius < 0) {
    BOOST_THROW_EXCEPTION(Exception("Cone parameters cannot be negative."));
  }
}

auto
Cone::generatePolyhedralSurface() const -> PolyhedralSurface
{
  if (m_polyhedral_surface) {
    return *m_polyhedral_surface;
  }

  // faces
  m_polyhedral_surface = PolyhedralSurface();

  // bottom face
  std::vector<Point> bottom_points;
  if (bottomRadius() > 0) {
    for (unsigned int i = 0; i < numRadial(); ++i) {
      double     angle = -2.0 * CGAL_PI * i / numRadial();
      Kernel::FT x     = bottomRadius() * std::cos(angle);
      Kernel::FT y     = bottomRadius() * std::sin(angle);
      bottom_points.emplace_back(x, y, 0);
    }
    bottom_points.push_back(bottom_points[0]);
    m_polyhedral_surface->addPatch(LineString(bottom_points));
  }

  // top face
  std::vector<Point> top_points;
  if (topRadius() > 0) {
    for (unsigned int i = 0; i < numRadial(); ++i) {
      double     angle = 2.0 * CGAL_PI * i / numRadial();
      Kernel::FT x     = topRadius() * std::cos(angle);
      Kernel::FT y     = topRadius() * std::sin(angle);
      top_points.emplace_back(x, y, height());
    }
    top_points.push_back(top_points[0]);
    m_polyhedral_surface->addPatch(LineString(top_points));
  }

  // side faces if top radius is a single point (true cone)
  if (topRadius() == 0 && bottomRadius() > 0) {
    for (unsigned int i = 0; i < numRadial(); ++i) {
      m_polyhedral_surface->addPatch(
          LineString({bottom_points[i], Point(0, 0, height()),
                      bottom_points[i + 1], bottom_points[i]}));
    }
  }

  // side faces if bottom radius is a single point (true cone upside-down)
  else if (bottomRadius() == 0 && topRadius() > 0) {
    for (unsigned int i = 0; i < numRadial(); ++i) {
      m_polyhedral_surface->addPatch(LineString({
          top_points[i],
          Point(0, 0, 0),
          top_points[i + 1],
          top_points[i],
      }));
    }
  }

  // side faces if we have a frustum (part of a cone)
  else if (bottomRadius() > 0 && topRadius() > 0) {
    for (unsigned int i = 0; i < numRadial(); ++i) {
      m_polyhedral_surface->addPatch(
          LineString({bottom_points[i], top_points[numRadial() - i],
                      top_points[numRadial() - 1 - i], bottom_points[i + 1],
                      bottom_points[i]}));
    }
  }

  return *m_polyhedral_surface;
}

auto
Cone::area3D(bool withDiscretization) const -> double
{
  double diff_radius  = std::abs(CGAL::to_double(bottomRadius() - topRadius()));
  double slant_height = std::sqrt(std::pow(diff_radius, 2) +
                                  std::pow(CGAL::to_double(height()), 2));
  if (withDiscretization) {
    double trapezoid_height =
        std::sqrt(std::pow(slant_height, 2) -
                  std::pow(diff_radius * std::sin(CGAL_PI / numRadial()), 2));
    return (std::sin(CGAL_PI / numRadial()) *
                CGAL::to_double((bottomRadius() + topRadius())) *
                trapezoid_height +
            std::sin(2.0 * CGAL_PI / numRadial()) *
                CGAL::to_double((bottomRadius() * bottomRadius() +
                                 topRadius() * topRadius())) /
                2.0) *
           numRadial();
  }
  return CGAL::to_double(
      CGAL_PI * ((bottomRadius() + topRadius()) * slant_height +
                 bottomRadius() * bottomRadius() + topRadius() * topRadius()));
}

auto
Cone::volume(bool withDiscretization) const -> double
{
  if (withDiscretization) {
    double bottom_surface = std::pow(CGAL::to_double(bottomRadius()), 2) *
                            numRadial() / 2.0 *
                            std::sin(2.0 * CGAL_PI / numRadial());
    double top_surface = std::pow(CGAL::to_double(topRadius()), 2) *
                         numRadial() / 2.0 *
                         std::sin(2.0 * CGAL_PI / numRadial());
    return CGAL::to_double(height()) / 3.0 *
           (bottom_surface + std::sqrt(bottom_surface * top_surface) +
            top_surface);
  }

  return CGAL::to_double(height() * CGAL_PI / 3.0 *
                         (bottomRadius() * bottomRadius() +
                          bottomRadius() * topRadius() +
                          topRadius() * topRadius()));
}

auto
Cone::height() const -> const Kernel::FT &
{
  return std::get<Kernel::FT>(m_parameters.at("height"));
}

auto
Cone::bottomRadius() const -> const Kernel::FT &
{
  return std::get<Kernel::FT>(m_parameters.at("bottom_radius"));
}

auto
Cone::topRadius() const -> const Kernel::FT &
{
  return std::get<Kernel::FT>(m_parameters.at("top_radius"));
}

auto
Cone::numRadial() const -> unsigned int
{
  return std::get<unsigned int>(m_parameters.at("num_radial"));
}
} // namespace SFCGAL
