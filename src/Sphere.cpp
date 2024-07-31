#include "SFCGAL/Sphere.h"
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <cmath>

namespace SFCGAL {
using Point_3      = Kernel::Point_3;
using Polyhedron_3 = CGAL::Polyhedron_3<Kernel>;
template <class HDS>
class Sphere_builder : public CGAL::Modifier_base<HDS> {
public:
  Sphere_builder(double radius, int num_vertical, int num_horizontal,
                 Point_3 center = Point_3(0, 0, 0))
      : radius(radius), num_vertical(num_vertical),
        num_horizontal(num_horizontal), center(center)
  {
  }

  void
  operator()(HDS &hds)
  {
    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

    int num_vertices = (num_vertical - 1) * num_horizontal + 2;
    int num_faces    = num_vertical * num_horizontal * 2;

    B.begin_surface(num_vertices, num_faces);

    // Ajouter le sommet supérieur
    B.add_vertex(Point_3(center.x(), center.y(), center.z() + radius));

    // Ajouter les sommets intermédiaires
    for (int i = 1; i < num_vertical; ++i) {
      double phi = M_PI * double(i) / double(num_vertical);
      double z   = radius * std::cos(phi);
      double r   = radius * std::sin(phi);
      for (int j = 0; j < num_horizontal; ++j) {
        double theta = 2 * M_PI * double(j) / double(num_horizontal);
        double x     = r * std::cos(theta);
        double y     = r * std::sin(theta);
        B.add_vertex(Point_3(center.x() + x, center.y() + y, center.z() + z));
      }
    }

    // Ajouter le sommet inférieur
    B.add_vertex(Point_3(center.x(), center.y(), center.z() - radius));

    // Ajouter les faces supérieures
    for (int j = 0; j < num_horizontal; ++j) {
      B.begin_facet();
      B.add_vertex_to_facet(0);
      B.add_vertex_to_facet(1 + (j + 1) % num_horizontal);
      B.add_vertex_to_facet(1 + j);
      B.end_facet();
    }

    // Ajouter les faces intermédiaires
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

    // Ajouter les faces inférieures
    int last_vertex = num_vertices - 1;
    int last_row    = 1 + (num_vertical - 2) * num_horizontal;
    for (int j = 0; j < num_horizontal; ++j) {
      B.begin_facet();
      B.add_vertex_to_facet(last_vertex);
      B.add_vertex_to_facet(last_row + j);
      B.add_vertex_to_facet(last_row + (j + 1) % num_horizontal);
      B.end_facet();
    }

    B.end_surface();
  }

private:
  double  radius;
  int     num_vertical;
  int     num_horizontal;
  Point_3 center;
};

Sphere::Sphere(const Kernel::FT &radius, const Point_3 &center,
               int num_vertical, int num_horizontal)
    : m_radius(radius), m_center(center), m_num_vertical(num_vertical),
      m_num_horizontal(num_horizontal)
{
}

Sphere &
Sphere::operator=(Sphere other)
{
  std::swap(m_radius, other.m_radius);
  std::swap(m_center, other.m_center);
  std::swap(m_num_vertical, other.m_num_vertical);
  std::swap(m_num_horizontal, other.m_num_horizontal);
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

Polyhedron_3
Sphere::generateSpherePolyhedron()
{
  Polyhedron_3                             P;
  Sphere_builder<Polyhedron_3::HalfedgeDS> builder(
      CGAL::to_double(m_radius), m_num_vertical, m_num_horizontal, m_center);
  P.delegate(builder);
  return P;
}

Polyhedron_3
Sphere::generatePolyhedron()
{
  if (!m_polyhedron) {
    m_polyhedron = generateSpherePolyhedron();
  }
  return *m_polyhedron;
}

std::vector<Point_3>
Sphere::generatePoints()
{
  if (!m_points) {
    m_points = generateSpherePoints();
  }
  return *m_points;
}

std::vector<Point_3>
Sphere::generateSpherePoints()
{
  std::vector<Point_3> points;
  points.reserve(m_num_vertical * m_num_horizontal);

  Kernel::FT d_lat = CGAL_PI / (m_num_vertical - 1);
  Kernel::FT d_lon = 2 * CGAL_PI / m_num_horizontal;

  for (int i = 0; i < m_num_vertical; ++i) {
    Kernel::FT lat = CGAL_PI / 2 - i * d_lat;
    Kernel::FT z   = m_radius * std::sin(CGAL::to_double(lat));
    Kernel::FT r   = m_radius * std::cos(CGAL::to_double(lat));
    for (int j = 0; j < m_num_horizontal; ++j) {
      Kernel::FT lon = j * d_lon;
      points.emplace_back(m_center.x() + r * std::cos(CGAL::to_double(lon)),
                          m_center.y() + r * std::sin(CGAL::to_double(lon)),
                          m_center.z() + z);
    }
  }
  return points;
}

} // namespace SFCGAL
