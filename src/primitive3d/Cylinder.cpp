// Copyright (c) 2024-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <utility>

#include "SFCGAL/numeric.h"
#include "SFCGAL/primitive3d/Cylinder.h"

namespace SFCGAL {

Cylinder::Cylinder(const Point_3 &base_center, const Vector_3 &axis,
                   const Kernel::FT &radius, const Kernel::FT &height,
                   unsigned int num_radial)
{
  m_parameters["base_center"] = base_center;
  m_parameters["axis"]        = axis;
  m_parameters["radius"]      = radius;
  m_parameters["height"]      = height;
  m_parameters["num_radial"]  = num_radial;

  Cylinder::validateParameters(m_parameters);
}

auto
Cylinder::primitiveType() const -> std::string
{
  return "Cylinder";
}

auto
Cylinder::primitiveTypeId() const -> PrimitiveType
{
  return PrimitiveType::TYPE_CYLINDER;
}

void
Cylinder::setBaseCenter(const Point_3 &base_center)
{
  validateAndSetParameter("base_center", base_center);
}

void
Cylinder::setAxis(const Vector_3 &axis)
{
  validateAndSetParameter("axis", axis);
}

void
Cylinder::setRadius(const Kernel::FT &radius)
{
  validateAndSetParameter("radius", radius);
}

void
Cylinder::setHeight(const Kernel::FT &height)
{
  validateAndSetParameter("height", height);
}

void
Cylinder::setNumRadial(unsigned int num)
{
  validateAndSetParameter("num_radial", num);
}

void
Cylinder::validateParameters(
    std::unordered_map<std::string, PrimitiveParameter> const &tempParameters)
    const
{
  const double radius =
      CGAL::to_double(std::get<Kernel::FT>(tempParameters.at("radius")));
  const double height =
      CGAL::to_double(std::get<Kernel::FT>(tempParameters.at("height")));
  const unsigned int num_radial =
      std::get<unsigned int>(tempParameters.at("num_radial"));

  if (radius <= 0.) {
    BOOST_THROW_EXCEPTION(Exception("Cylinder radius cannot be negative."));
  } else if (height <= 0.) {
    BOOST_THROW_EXCEPTION(Exception("Cylinder height cannot be negative."));
  } else if (num_radial < 3) {
    BOOST_THROW_EXCEPTION(
        Exception("Cylinder requires at least 3 radial segments."));
  }
}

void
Cylinder::invalidateCache()
{
  Primitive::invalidateCache();
  m_polyhedron.reset();
  m_surface_mesh.reset();
  m_polyhedral_surface.reset();
}

auto
Cylinder::normalize(const Vector_3 &vector) -> Vector_3
{
  double length = std::sqrt(CGAL::to_double(vector.squared_length()));
  if (length < EPSILON) {
    return vector;
  }
  return vector / length;
}

auto
Cylinder::generatePolyhedron() const -> Polyhedron_3
{
  if (m_polyhedron) {
    return *m_polyhedron;
  }

  Polyhedron_3   poly;
  Surface_mesh_3 sm = generateSurfaceMesh();
  CGAL::copy_face_graph(sm, poly);
  m_polyhedron = poly;
  return poly;
}

auto
Cylinder::generateSurfaceMesh() const -> Surface_mesh_3
{
  if (m_surface_mesh) {
    return *m_surface_mesh;
  }

  Surface_mesh_3 mesh;

  Vector_3 normalized_axis = normalize(axis());

  // Find a perpendicular vector
  Vector_3 perpendicular;
  if (CGAL::abs(normalized_axis.z()) < 0.9) {
    // If axis is not too close to Z, use Z cross axis
    perpendicular =
        normalize(CGAL::cross_product(normalized_axis, Vector_3(0, 0, 1)));
  } else {
    // If axis is close to Z, use Y cross axis
    perpendicular =
        normalize(CGAL::cross_product(normalized_axis, Vector_3(0, 1, 0)));
  }
  Vector_3 perpendicular2 =
      normalize(CGAL::cross_product(normalized_axis, perpendicular));

  std::vector<Surface_mesh_3::Vertex_index> base_vertices;
  std::vector<Surface_mesh_3::Vertex_index> top_vertices;

  // Create vertices for the base and top
  for (unsigned int i = 0; i < numRadial(); ++i) {
    double   angle  = 2.0 * M_PI * i / numRadial();
    Vector_3 offset = radius() * (std::cos(angle) * perpendicular +
                                  std::sin(angle) * perpendicular2);
    base_vertices.push_back(mesh.add_vertex(baseCenter() + offset));
    top_vertices.push_back(
        mesh.add_vertex(baseCenter() + offset + height() * normalized_axis));
  }

  // Add side faces
  for (unsigned int i = 0; i < numRadial(); ++i) {
    unsigned int next = (i + 1) % numRadial();
    mesh.add_face(base_vertices[i], top_vertices[i], top_vertices[next]);
    mesh.add_face(base_vertices[i], top_vertices[next], base_vertices[next]);
  }

  // Add base and top faces
  Surface_mesh_3::Vertex_index base_center = mesh.add_vertex(baseCenter());
  Surface_mesh_3::Vertex_index top_center =
      mesh.add_vertex(baseCenter() + height() * normalized_axis);

  for (unsigned int i = 0; i < numRadial(); ++i) {
    unsigned int next = (i + 1) % numRadial();
    mesh.add_face(base_center, base_vertices[i], base_vertices[next]);
    mesh.add_face(top_center, top_vertices[next], top_vertices[i]);
  }
  m_surface_mesh = mesh;
  return mesh;
}

auto
Cylinder::generatePolyhedralSurface() const -> PolyhedralSurface
{
  if (m_polyhedral_surface) {
    return *m_polyhedral_surface;
  }

  m_polyhedral_surface = PolyhedralSurface(generateSurfaceMesh());
  return *m_polyhedral_surface;
}

auto
Cylinder::volume(bool /*withDiscretization*/) const -> double
{
  return CGAL::to_double(radius() * radius() * height() * CGAL_PI);
}

auto
Cylinder::area3D(bool /*withDiscretization*/) const -> double
{
  return CGAL::to_double(2 * radius() * radius() * CGAL_PI +
                         2 * radius() * height() * CGAL_PI);
}

} // namespace SFCGAL
