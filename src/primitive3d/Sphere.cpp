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
#include <map>
#include <stdexcept>
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
  Sphere_builder(double radius, unsigned int num_subdivisions, Point_3 center,
                 const Kernel::Vector_3 &direction)
      : radius(radius), num_subdivisions(num_subdivisions),
        center(std::move(center)), direction(normalizeVector(direction))
  {
  }

  void
  operator()(HDS &hds) override
  {
    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

    // Create icosahedron vertices
    const double phi      = (1.0 + std::sqrt(5.0)) / 2.0; // golden ratio
    const double inv_norm = 1.0 / std::sqrt(1.0 + (phi * phi));

    // 12 base icosahedron vertices
    std::vector<Point_3> base_vertices = {
        Point_3(-1 * inv_norm, phi * inv_norm, 0),
        Point_3(1 * inv_norm, phi * inv_norm, 0),
        Point_3(-1 * inv_norm, -phi * inv_norm, 0),
        Point_3(1 * inv_norm, -phi * inv_norm, 0),
        Point_3(0, -1 * inv_norm, phi * inv_norm),
        Point_3(0, 1 * inv_norm, phi * inv_norm),
        Point_3(0, -1 * inv_norm, -phi * inv_norm),
        Point_3(0, 1 * inv_norm, -phi * inv_norm),
        Point_3(phi * inv_norm, 0, -1 * inv_norm),
        Point_3(phi * inv_norm, 0, 1 * inv_norm),
        Point_3(-phi * inv_norm, 0, -1 * inv_norm),
        Point_3(-phi * inv_norm, 0, 1 * inv_norm)};

    // 20 base icosahedron faces
    std::vector<std::array<int, 3>> base_faces = {
        {0, 11, 5}, {0, 5, 1},  {0, 1, 7},   {0, 7, 10}, {0, 10, 11},
        {1, 5, 9},  {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
        {3, 9, 4},  {3, 4, 2},  {3, 2, 6},   {3, 6, 8},  {3, 8, 9},
        {4, 9, 5},  {2, 4, 11}, {6, 2, 10},  {8, 6, 7},  {9, 8, 1}};

    // Subdivide the icosahedron
    std::vector<Point_3>            vertices = base_vertices;
    std::vector<std::array<int, 3>> faces    = base_faces;

    for (unsigned int level = 0; level < num_subdivisions; ++level) {
      subdivideIcosahedron(vertices, faces);
    }

    // Project all vertices to sphere surface, scale, orient and translate
    for (auto &vertex : vertices) {
      // Normalize to unit sphere
      Kernel::Vector_3 vec(vertex.x(), vertex.y(), vertex.z());
      double length = std::sqrt(CGAL::to_double(vec.squared_length()));
      if (length > 1e-10) {
        vec = vec / length;
      }

      // Apply radius scaling
      vec = vec * radius;

      // Apply orientation (simple case - keep Z-up for now)
      // In future, could apply rotation matrix based on direction vector

      // Translate to center
      vertex = Point_3(center.x() + vec.x(), center.y() + vec.y(),
                       center.z() + vec.z());
    }

    B.begin_surface(vertices.size(), faces.size());

    try {
      // Add vertices
      for (const auto &vertex : vertices) {
        B.add_vertex(vertex);
      }

      // Add faces
      for (const auto &face : faces) {
        B.begin_facet();
        B.add_vertex_to_facet(face[0]);
        B.add_vertex_to_facet(face[1]);
        B.add_vertex_to_facet(face[2]);
        B.end_facet();
      }

      B.end_surface();

      if (B.error()) {
        throw std::runtime_error("CGAL polyhedron builder reported an error "
                                 "during icosahedron construction");
      }
    } catch (const std::exception &e) {
      if (!B.error()) {
        B.rollback();
      }
      throw std::runtime_error(
          std::string("Failed to build subdivided icosahedron: ") + e.what());
    }
  }

private:
  double           radius;
  unsigned int     num_subdivisions;
  Point_3          center;
  Kernel::Vector_3 direction;

  void
  subdivideIcosahedron(std::vector<Point_3>            &vertices,
                       std::vector<std::array<int, 3>> &faces)
  {
    std::vector<std::array<int, 3>>    new_faces;
    std::map<std::pair<int, int>, int> edge_to_vertex;

    // Helper function to get or create a midpoint vertex
    auto getMidpoint = [&](int vertexA, int vertexB) -> int {
      auto key = std::make_pair(std::min(vertexA, vertexB),
                                std::max(vertexA, vertexB));
      auto it  = edge_to_vertex.find(key);
      if (it != edge_to_vertex.end()) {
        return it->second;
      }

      // Create new midpoint vertex
      Point_3 point1 = vertices[vertexA];
      Point_3 point2 = vertices[vertexB];
      Point_3 midpoint((point1.x() + point2.x()) / 2.0,
                       (point1.y() + point2.y()) / 2.0,
                       (point1.z() + point2.z()) / 2.0);

      int new_index = static_cast<int>(vertices.size());
      vertices.push_back(midpoint);
      edge_to_vertex[key] = new_index;
      return new_index;
    };

    // Subdivide each triangle into 4 smaller triangles
    for (const auto &face : faces) {
      int vertex0 = face[0];
      int vertex1 = face[1];
      int vertex2 = face[2];

      // Get midpoints of each edge
      int m01 = getMidpoint(vertex0, vertex1);
      int m12 = getMidpoint(vertex1, vertex2);
      int m20 = getMidpoint(vertex2, vertex0);

      // Create 4 new triangles
      new_faces.push_back({vertex0, m01, m20});
      new_faces.push_back({vertex1, m12, m01});
      new_faces.push_back({vertex2, m20, m12});
      new_faces.push_back({m01, m12, m20});
    }

    faces = new_faces;
  }
};

/// @} end of private section

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

Sphere::Sphere(const Kernel::FT &radius, const Kernel::Point_3 &center,
               unsigned int num_subdivisions, const Kernel::Vector_3 &direction)
{
  m_parameters["radius"]           = Kernel::FT(radius);
  m_parameters["num_subdivisions"] = num_subdivisions;
  m_parameters["center"]           = center;
  m_parameters["direction"]        = normalizeVector(direction);

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
    BOOST_THROW_EXCEPTION(Exception("Sphere radius must be positive."));
  }

  const unsigned int num_subdivisions =
      std::get<unsigned int>(tempParameters.at("num_subdivisions"));

  if (num_subdivisions > 6) {
    BOOST_THROW_EXCEPTION(Exception(
        "Sphere subdivisions should not exceed 6 (too many vertices)."));
  }
}

// Generate the polyhedron representation of the sphere
auto
Sphere::generateSpherePolyhedron() const -> Polyhedron_3
{
  Polyhedron_3                             P;
  Sphere_builder<Polyhedron_3::HalfedgeDS> builder(
      CGAL::to_double(radius()), numSubdivisions(), center(), direction());
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

  try {
    Polyhedron_3 polyhedron = generatePolyhedron();

    // Validate the polyhedron before converting to PolyhedralSurface
    if (polyhedron.empty()) {
      throw std::runtime_error("Generated sphere polyhedron is empty");
    }

    if (!polyhedron.is_valid()) {
      throw std::runtime_error("Generated sphere polyhedron is invalid");
    }

    if (!polyhedron.is_closed()) {
      throw std::runtime_error("Generated sphere polyhedron is not closed");
    }

    m_polyhedral_surface = PolyhedralSurface(polyhedron);

    return *m_polyhedral_surface;
  } catch (const std::exception &e) {
    throw std::runtime_error(std::string("Sphere surface generation failed: ") +
                             e.what());
  }
}

// Generate points on the sphere's surface from icosahedron vertices
auto
Sphere::generateSpherePoints() const -> std::vector<Point_3>
{
  // Generate the same vertices as the polyhedron
  Polyhedron_3         poly = generatePolyhedron();
  std::vector<Point_3> points;
  points.reserve(poly.size_of_vertices());

  for (auto vit = poly.vertices_begin(); vit != poly.vertices_end(); ++vit) {
    points.emplace_back(vit->point());
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
