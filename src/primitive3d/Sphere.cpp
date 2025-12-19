// Copyright (c) 2024-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/primitive3d/Sphere.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/primitive3d/Primitive.h"

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Vector_3.h>
#include <cmath>
#include <utility>

namespace SFCGAL {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

/**
 * @brief Helper class for building the sphere polyhedron
 */
template <class HDS>
class Sphere_builder : public CGAL::Modifier_base<HDS> {
public:
  Sphere_builder(double radius, unsigned int num_vertical,
                 unsigned int num_horizontal, Point_3 center,
                 const Kernel::Vector_3 &direction)
      : radius(radius), num_vertical(num_vertical),
        num_horizontal(num_horizontal), center(std::move(center)),
        direction(normalizeVector(direction))
  {
  }

  void
  operator()(HDS &hds) override
  {
    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

    int num_vertices = ((num_vertical - 1) * num_horizontal) + 2;
    int num_faces    = num_vertical * num_horizontal * 2;

    B.begin_surface(num_vertices, num_faces);

    // Add vertices
    addVertices(B);

    // Add faces
    addTopFaces(B);
    addMiddleFaces(B);
    addBottomFaces(B);

    B.end_surface();
  }

private:
  // Function to get an orthogonal vector in the XY plane
  auto
  get_orthogonal_vector(const Kernel::Vector_3 &vec) -> Kernel::Vector_3
  {
    if (vec.x() != 0 || vec.y() != 0) {
      return Kernel::Vector_3(-vec.y(), vec.x(), 0);
    }
    return Kernel::Vector_3(0, -vec.z(), vec.y());
  }

  void
  addVertices(CGAL::Polyhedron_incremental_builder_3<HDS> &B)
  {
    Kernel::Vector_3 v1 = normalizeVector(get_orthogonal_vector(direction));
    Kernel::Vector_3 v2 = normalizeVector(CGAL::cross_product(direction, v1));

    // Add top vertex
    B.add_vertex(Point_3(center + direction * radius));

    // Add middle vertices
    for (unsigned int i = 1; i < num_vertical; ++i) {
      double phi = M_PI * double(i) / double(num_vertical);
      double z   = radius * std::cos(phi);
      double r   = radius * std::sin(phi);
      for (unsigned int j = 0; j < num_horizontal; ++j) {
        double           theta = 2 * M_PI * double(j) / double(num_horizontal);
        Kernel::Vector_3 point_vec =
            r * (std::cos(theta) * v1 + std::sin(theta) * v2) + z * direction;
        B.add_vertex(Point_3(center + point_vec));
      }
    }

    // Add bottom vertex
    B.add_vertex(Point_3(center - direction * radius));
  }

  void
  addTopFaces(CGAL::Polyhedron_incremental_builder_3<HDS> &B)
  {
    for (unsigned int j = 0; j < num_horizontal; ++j) {
      B.begin_facet();
      B.add_vertex_to_facet(0);
      B.add_vertex_to_facet(1 + ((j + 1) % num_horizontal));
      B.add_vertex_to_facet(1 + j);
      B.end_facet();
    }
  }

  void
  addMiddleFaces(CGAL::Polyhedron_incremental_builder_3<HDS> &B)
  {
    for (unsigned int i = 1; i < num_vertical - 1; ++i) {
      for (unsigned int j = 0; j < num_horizontal; ++j) {
        int current = 1 + ((i - 1) * num_horizontal) + j;
        int next = 1 + ((i - 1) * num_horizontal) + ((j + 1) % num_horizontal);
        int below_current = 1 + (i * num_horizontal) + j;
        int below_next = 1 + (i * num_horizontal) + ((j + 1) % num_horizontal);

        B.begin_facet();
        B.add_vertex_to_facet(current);
        B.add_vertex_to_facet(next);
        B.add_vertex_to_facet(below_next);
        B.end_facet();

        B.begin_facet();
        B.add_vertex_to_facet(current);
        B.add_vertex_to_facet(below_next);
        B.add_vertex_to_facet(below_current);
        B.end_facet();
      }
    }
  }

  void
  addBottomFaces(CGAL::Polyhedron_incremental_builder_3<HDS> &B)
  {
    unsigned int last_vertex = ((num_vertical - 1) * num_horizontal) + 1;
    unsigned int last_row    = 1 + ((num_vertical - 2) * num_horizontal);
    for (unsigned int j = 0; j < num_horizontal; ++j) {
      B.begin_facet();
      B.add_vertex_to_facet(last_vertex);
      B.add_vertex_to_facet(last_row + j);
      B.add_vertex_to_facet(last_row + ((j + 1) % num_horizontal));
      B.end_facet();
    }
  }

  double           radius;
  unsigned int     num_vertical;
  unsigned int     num_horizontal;
  Point_3          center;
  Kernel::Vector_3 direction;
};

/// @} end of private section

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

Sphere::Sphere(const Kernel::FT &radius, const Kernel::Point_3 &center,
               unsigned int num_vertical, unsigned int num_horizontal,
               const Kernel::Vector_3 &direction)
{
  m_parameters["radius"]         = Kernel::FT(radius);
  m_parameters["num_vertical"]   = num_vertical;
  m_parameters["num_horizontal"] = num_horizontal;
  m_parameters["center"]         = center;
  m_parameters["direction"]      = normalizeVector(direction);

  Sphere::validateParameters(m_parameters);
}

auto
Sphere::primitiveType() const -> std::string
{
  return "Sphere";
}

auto
Sphere::primitiveTypeId() const -> PrimitiveType
{
  return PrimitiveType::TYPE_SPHERE;
}

void
Sphere::invalidateCache()
{
  Primitive::invalidateCache();
  m_polyhedron.reset();
  m_points.reset();
}

void
Sphere::validateParameters(
    std::unordered_map<std::string, PrimitiveParameter> const &tempParameters)
    const
{
  const double radius =
      CGAL::to_double(std::get<Kernel::FT>(tempParameters.at("radius")));

  if (radius <= 0.) {
    BOOST_THROW_EXCEPTION(Exception("Sphere radius cannot be negative."));
  }
}

// Generate the polyhedron representation of the sphere
auto
Sphere::generateSpherePolyhedron() const -> Polyhedron_3
{
  Polyhedron_3                             P;
  Sphere_builder<Polyhedron_3::HalfedgeDS> builder(
      CGAL::to_double(radius()), numVertical(), numHorizontal(), center(),
      direction());
  P.delegate(builder);
  return P;
}

auto
Sphere::generatePolyhedron() const -> Polyhedron_3
{
  if (!m_polyhedron) {
    m_polyhedron = generateSpherePolyhedron();
  }
  return *m_polyhedron;
}

auto
Sphere::generatePoints() const -> std::vector<Point_3>
{
  if (!m_points) {
    m_points = generateSpherePoints();
  }
  return *m_points;
}

auto
Sphere::generatePolyhedralSurface() const -> PolyhedralSurface
{
  if (m_polyhedral_surface) {
    return *m_polyhedral_surface;
  }

  m_polyhedral_surface = PolyhedralSurface(generatePolyhedron());
  return *m_polyhedral_surface;
}

// Generate points on the sphere's surface
auto
Sphere::generateSpherePoints() const -> std::vector<Point_3>
{
  std::vector<Point_3> points;
  points.reserve(static_cast<unsigned long>(numVertical()) * numHorizontal());

  Kernel::Vector_3 vec1 = normalizeVector(get_orthogonal_vector(direction()));
  Kernel::Vector_3 vec2 =
      normalizeVector(CGAL::cross_product(direction(), vec1));

  Kernel::FT d_lat = CGAL_PI / static_cast<double>(numVertical() - 1);
  Kernel::FT d_lon = 2 * CGAL_PI / static_cast<double>(numHorizontal());

  for (unsigned int i = 0; i < numVertical(); ++i) {
    Kernel::FT lat        = CGAL_PI / 2 - static_cast<int>(i) * d_lat;
    Kernel::FT z          = radius() * std::sin(CGAL::to_double(lat));
    Kernel::FT lat_radius = radius() * std::cos(CGAL::to_double(lat));
    for (unsigned int j = 0; j < numHorizontal(); ++j) {
      Kernel::FT       lon = static_cast<int>(j) * d_lon;
      Kernel::Vector_3 point_vec =
          lat_radius * (std::cos(CGAL::to_double(lon)) * vec1 +
                        std::sin(CGAL::to_double(lon)) * vec2) +
          z * direction();
      points.emplace_back(center() + point_vec);
    }
  }
  return points;
}

auto
Sphere::volume(bool /*withDiscretization*/) const -> double
{
  return CGAL::to_double((4.0 / 3.0) * radius() * radius() * radius() *
                         CGAL_PI);
}

auto
Sphere::area3D(bool /*withDiscretization*/) const -> double
{
  return CGAL::to_double(4 * radius() * radius() * CGAL_PI);
}

} // namespace SFCGAL
