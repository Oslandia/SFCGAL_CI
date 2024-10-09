#include <utility>

#include "SFCGAL/Cylinder.h"

namespace SFCGAL {

Cylinder::Cylinder(const Point_3 &base_center, const Vector_3 &axis,
                   const Kernel::FT &radius, const Kernel::FT &height,
                   int num_radial)
    : m_base_center(base_center), m_axis(axis), m_radius(radius),
      m_height(height), m_num_radial(num_radial)
{
}

auto
Cylinder::operator=(Cylinder other) -> Cylinder &
{
  std::swap(m_base_center, other.m_base_center);
  std::swap(m_axis, other.m_axis);
  std::swap(m_radius, other.m_radius);
  std::swap(m_height, other.m_height);
  std::swap(m_num_radial, other.m_num_radial);
  std::swap(m_polyhedron, other.m_polyhedron);
  std::swap(m_surface_mesh, other.m_surface_mesh);
  return *this;
}

void
Cylinder::setBaseCenter(const Point_3 &base_center)
{
  m_base_center = base_center;
  invalidateCache();
}

void
Cylinder::setAxis(const Vector_3 &axis)
{
  m_axis = axis;
  invalidateCache();
}

void
Cylinder::setRadius(const Kernel::FT &radius)
{
  m_radius = radius;
  invalidateCache();
}

void
Cylinder::setHeight(const Kernel::FT &height)
{
  m_height = height;
  invalidateCache();
}

void
Cylinder::setNumRadial(int num)
{
  m_num_radial = num;
  invalidateCache();
}

void
Cylinder::invalidateCache()
{
  m_polyhedron.reset();
  m_surface_mesh.reset();
}

auto
Cylinder::normalize(const Vector_3 &v) -> Cylinder::Vector_3
{
  double length = std::sqrt(CGAL::to_double(v.squared_length()));
  if (length < 1e-8) {
    return v;
  }
  return v / length;
}

auto
Cylinder::generatePolyhedron() -> Cylinder::Polyhedron_3
{
  if (m_polyhedron) {
    return *m_polyhedron;
  }

  Polyhedron_3 poly;
  Surface_mesh sm = generateSurfaceMesh();
  CGAL::copy_face_graph(sm, poly);
  m_polyhedron = poly;
  return poly;
}

auto
Cylinder::generateSurfaceMesh() -> Cylinder::Surface_mesh
{
  if (m_surface_mesh) {
    return *m_surface_mesh;
  }

  Surface_mesh mesh;

  Vector_3 normalized_axis = normalize(m_axis);
  Vector_3 perpendicular =
      normalize(CGAL::cross_product(normalized_axis, Vector_3(0, 0, 1)));
  if (perpendicular.squared_length() < 1e-8) {
    perpendicular =
        normalize(CGAL::cross_product(normalized_axis, Vector_3(0, 1, 0)));
  }
  Vector_3 perpendicular2 = CGAL::cross_product(normalized_axis, perpendicular);

  std::vector<Surface_mesh::Vertex_index> base_vertices;
  std::vector<Surface_mesh::Vertex_index> top_vertices;

  // Create vertices for the base and top
  for (int i = 0; i < m_num_radial; ++i) {
    double   angle  = 2.0 * M_PI * i / m_num_radial;
    Vector_3 offset = m_radius * (std::cos(angle) * perpendicular +
                                  std::sin(angle) * perpendicular2);
    base_vertices.push_back(mesh.add_vertex(m_base_center + offset));
    top_vertices.push_back(
        mesh.add_vertex(m_base_center + offset + m_height * normalized_axis));
  }

  // Add side faces
  for (int i = 0; i < m_num_radial; ++i) {
    int next = (i + 1) % m_num_radial;
    mesh.add_face(base_vertices[i], top_vertices[i], top_vertices[next]);
    mesh.add_face(base_vertices[i], top_vertices[next], base_vertices[next]);
  }

  // Add base and top faces
  Surface_mesh::Vertex_index base_center = mesh.add_vertex(m_base_center);
  Surface_mesh::Vertex_index top_center =
      mesh.add_vertex(m_base_center + m_height * normalized_axis);

  for (int i = 0; i < m_num_radial; ++i) {
    int next = (i + 1) % m_num_radial;
    mesh.add_face(base_center, base_vertices[i], base_vertices[next]);
    mesh.add_face(top_center, top_vertices[next], top_vertices[i]);
  }
  m_surface_mesh = mesh;
  return mesh;
}

} // namespace SFCGAL
