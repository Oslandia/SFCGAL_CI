#include "SFCGAL/Sphere.h"
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Vector_3.h>
#include <cmath>
#include <utility>

namespace SFCGAL {

/**
 * @brief Helper class for building the sphere polyhedron
 */
template <class HDS>
class Sphere_builder : public CGAL::Modifier_base<HDS> {
public:
  Sphere_builder(double radius, int num_vertical, int num_horizontal,
                 Point_3 center, const Kernel::Vector_3 &direction)
      : radius(radius), num_vertical(num_vertical),
        num_horizontal(num_horizontal), center(std::move(center)),
        direction(normalizeVector(direction))
  {
  }

  void
  operator()(HDS &hds) override
  {
    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

    int num_vertices = (num_vertical - 1) * num_horizontal + 2;
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
    for (int i = 1; i < num_vertical; ++i) {
      double phi = M_PI * double(i) / double(num_vertical);
      double z   = radius * std::cos(phi);
      double r   = radius * std::sin(phi);
      for (int j = 0; j < num_horizontal; ++j) {
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
    for (int j = 0; j < num_horizontal; ++j) {
      B.begin_facet();
      B.add_vertex_to_facet(0);
      B.add_vertex_to_facet(1 + (j + 1) % num_horizontal);
      B.add_vertex_to_facet(1 + j);
      B.end_facet();
    }
  }

  void
  addMiddleFaces(CGAL::Polyhedron_incremental_builder_3<HDS> &B)
  {
    for (int i = 1; i < num_vertical - 1; ++i) {
      for (int j = 0; j < num_horizontal; ++j) {
        int current = 1 + (i - 1) * num_horizontal + j;
        int next    = 1 + (i - 1) * num_horizontal + (j + 1) % num_horizontal;
        int below_current = 1 + i * num_horizontal + j;
        int below_next    = 1 + i * num_horizontal + (j + 1) % num_horizontal;

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
    int last_vertex = (num_vertical - 1) * num_horizontal + 1;
    int last_row    = 1 + (num_vertical - 2) * num_horizontal;
    for (int j = 0; j < num_horizontal; ++j) {
      B.begin_facet();
      B.add_vertex_to_facet(last_vertex);
      B.add_vertex_to_facet(last_row + j);
      B.add_vertex_to_facet(last_row + (j + 1) % num_horizontal);
      B.end_facet();
    }
  }

  double           radius;
  int              num_vertical;
  int              num_horizontal;
  Point_3          center;
  Kernel::Vector_3 direction;
};

Sphere::Sphere(const Kernel::FT &radius, const Kernel::Point_3 &center,
               int num_vertical, int num_horizontal,
               const Kernel::Vector_3 &direction)
    : m_radius(std::move(radius)), m_center(std::move(center)),
      m_num_vertical(num_vertical), m_num_horizontal(num_horizontal),
      m_direction(normalizeVector(direction))
{
}

auto
Sphere::operator=(Sphere other) -> Sphere &
{
  std::swap(m_radius, other.m_radius);
  std::swap(m_center, other.m_center);
  std::swap(m_num_vertical, other.m_num_vertical);
  std::swap(m_num_horizontal, other.m_num_horizontal);
  std::swap(m_direction, other.m_direction);
  std::swap(m_polyhedron, other.m_polyhedron);
  std::swap(m_points, other.m_points);
  return *this;
}

void
Sphere::invalidateCache()
{
  m_polyhedron.reset();
  m_points.reset();
}

// Generate the polyhedron representation of the sphere
auto
Sphere::generateSpherePolyhedron() -> Polyhedron_3
{
  Polyhedron_3                             P;
  Sphere_builder<Polyhedron_3::HalfedgeDS> builder(
      CGAL::to_double(m_radius), m_num_vertical, m_num_horizontal, m_center,
      m_direction);
  P.delegate(builder);
  return P;
}

auto
Sphere::generatePolyhedron() -> Polyhedron_3
{
  if (!m_polyhedron) {
    m_polyhedron = generateSpherePolyhedron();
  }
  return *m_polyhedron;
}

auto
Sphere::generatePoints() -> std::vector<Point_3>
{
  if (!m_points) {
    m_points = generateSpherePoints();
  }
  return *m_points;
}

// Generate points on the sphere's surface
auto
Sphere::generateSpherePoints() -> std::vector<Point_3>
{
  std::vector<Point_3> points;
  points.reserve(static_cast<size_t>(m_num_vertical) *
                 static_cast<size_t>(m_num_horizontal));

  Kernel::Vector_3 v1 = normalizeVector(get_orthogonal_vector(m_direction));
  Kernel::Vector_3 v2 = normalizeVector(CGAL::cross_product(m_direction, v1));

  Kernel::FT d_lat = CGAL_PI / (m_num_vertical - 1);
  Kernel::FT d_lon = 2 * CGAL_PI / m_num_horizontal;

  for (int i = 0; i < m_num_vertical; ++i) {
    Kernel::FT lat = CGAL_PI / 2 - i * d_lat;
    Kernel::FT z   = m_radius * std::sin(CGAL::to_double(lat));
    Kernel::FT r   = m_radius * std::cos(CGAL::to_double(lat));
    for (int j = 0; j < m_num_horizontal; ++j) {
      Kernel::FT       lon       = j * d_lon;
      Kernel::Vector_3 point_vec = r * (std::cos(CGAL::to_double(lon)) * v1 +
                                        std::sin(CGAL::to_double(lon)) * v2) +
                                   z * m_direction;
      points.emplace_back(m_center + point_vec);
    }
  }
  return points;
}

} // namespace SFCGAL
