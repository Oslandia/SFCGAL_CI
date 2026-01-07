// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/buffer3D.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/algorithm/minkowskiSum3D.h"
#include "SFCGAL/algorithm/union.h"
#include "SFCGAL/detail/GeometrySet.h"
#include "SFCGAL/numeric.h"
#include "SFCGAL/primitive3d/Cylinder.h"
#include "SFCGAL/primitive3d/Sphere.h"
#include "SFCGAL/triangulate/triangulatePolygon.h"
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/minkowski_sum_3.h>

namespace PMP = CGAL::Polygon_mesh_processing;

namespace SFCGAL::algorithm {

Buffer3D::Buffer3D(const Geometry &inputGeometry, double radius, int segments)
    : _radius(radius), _segments(segments)
{
  if (inputGeometry.is<Point>()) {
    _inputPoints.push_back(inputGeometry.as<Point>().toPoint_3());
  } else if (inputGeometry.is<LineString>()) {
    const auto &ls = inputGeometry.as<LineString>();
    for (size_t i = 0; i < ls.numPoints(); ++i) {
      _inputPoints.push_back(ls.pointN(i).toPoint_3());
    }
  } else {
    throw std::invalid_argument("Input geometry must be a Point or LineString");
  }
}

auto
Buffer3D::compute(BufferType type) const -> std::unique_ptr<PolyhedralSurface>
{
  if (_inputPoints.size() == 1) {
    return computePointBuffer();
  }

  switch (type) {
  case ROUND:
    return computeRoundBuffer();
  case CYLSPHERE:
    return computeCylSphereBuffer();
  case FLAT:
    return computeFlatBuffer();
  default:
    throw std::invalid_argument("Invalid buffer type");
  }
}

auto
Buffer3D::computePointBuffer() const -> std::unique_ptr<PolyhedralSurface>
{
  Kernel::Point_3 center(_inputPoints[0].x(), _inputPoints[0].y(),
                         _inputPoints[0].z());
  // Convert segments to subdivision level for icosahedron
  unsigned int subdivision_level =
      std::max(1U, static_cast<unsigned int>(_segments / 16));
  Sphere sphere(_radius, center, subdivision_level);
  return std::make_unique<PolyhedralSurface>(
      sphere.generatePolyhedralSurface());
}

auto
Buffer3D::computeRoundBuffer() const -> std::unique_ptr<PolyhedralSurface>
{
  using point_iterator = Point_3 *;
  using point_range    = std::pair<point_iterator, point_iterator>;
  using polyline       = std::list<point_range>;
  using Nef_polyhedron = CGAL::Nef_polyhedron_3<Kernel>;

  // Create sphere
  Point_3 center(0, 0, 0);
  // Convert segments to subdivision level for icosahedron
  unsigned int subdivision_level =
      std::max(1U, static_cast<unsigned int>(_segments / 16));
  SFCGAL::Sphere sphere(_radius, center, subdivision_level);

  // Generate polyhedron from sphere
  CGAL::Polyhedron_3<Kernel> spherePolyhedron = sphere.generatePolyhedron();

  // Convert Polyhedron to a Nef_polyhedron
  Nef_polyhedron N0(spherePolyhedron);

  // Create a polyline from _inputPoints
  polyline poly;
  if (!_inputPoints.empty()) {
    // Create a non-const copy of the points
    std::vector<Point_3> points_copy(_inputPoints.begin(), _inputPoints.end());
    poly.emplace_back(&points_copy.front(), &points_copy.back() + 1);

    // Create a Nef_polyhedron from the polyline
    Nef_polyhedron N1(poly.begin(), poly.end(),
                      Nef_polyhedron::Polylines_tag());

    // Perform Minkowski sum
    Nef_polyhedron result = CGAL::minkowski_sum_3(N0, N1);

    // Convert result to SFCGAL::PolyhedralSurface
    SFCGAL::detail::MarkedPolyhedron out;
    result.convert_to_polyhedron(out);
    return std::make_unique<PolyhedralSurface>(out);
  } // If _inputPoints is empty, return an empty PolyhedralSurface
  return std::make_unique<PolyhedralSurface>();
}

auto
Buffer3D::computeCylSphereBuffer() const -> std::unique_ptr<PolyhedralSurface>
{
  using Nef_polyhedron = CGAL::Nef_polyhedron_3<Kernel>;
  Nef_polyhedron result;

  // Add a sphere at the first point of the line
  if (!_inputPoints.empty()) {
    Kernel::Point_3 start_sphere_center(
        _inputPoints[0].x(), _inputPoints[0].y(), _inputPoints[0].z());
    Kernel::Vector_3 start_direction;
    if (_inputPoints.size() > 1) {
      start_direction = _inputPoints[1] - _inputPoints[0];
    } else {
      start_direction =
          Kernel::Vector_3(0, 0, 1); // Default direction if only one point
    }

    unsigned int subdivision_level =
        std::max(1U, static_cast<unsigned int>(_segments / 16));
    Sphere start_sphere(_radius, start_sphere_center, subdivision_level,
                        start_direction);
    CGAL::Polyhedron_3<Kernel> start_sphere_poly =
        start_sphere.generatePolyhedron();
    Nef_polyhedron start_sphere_nef(start_sphere_poly);

    result = start_sphere_nef;
  }

  // Create a cylinder and spheres for each segment of the line
  for (size_t i = 0; i < _inputPoints.size() - 1; ++i) {
    // Create a cylinder between each successive point
    Kernel::Vector_3 axis(_inputPoints[i + 1].x() - _inputPoints[i].x(),
                          _inputPoints[i + 1].y() - _inputPoints[i].y(),
                          _inputPoints[i + 1].z() - _inputPoints[i].z());
    Kernel::FT      height = CGAL::sqrt(CGAL::to_double(axis.squared_length()));
    Kernel::Point_3 base(_inputPoints[i].x(), _inputPoints[i].y(),
                         _inputPoints[i].z());
    Cylinder        cyl(base, axis, _radius, height, _segments);
    CGAL::Polyhedron_3<Kernel> cyl_poly = cyl.generatePolyhedron();
    Nef_polyhedron             cyl_nef(cyl_poly);

    result = result.join(cyl_nef);

    // Add a sphere at the junctions (rounded corners)
    if (i < _inputPoints.size() - 1) {
      Kernel::Point_3  sphereCenter(_inputPoints[i + 1].x(),
                                    _inputPoints[i + 1].y(),
                                    _inputPoints[i + 1].z());
      Kernel::Vector_3 sphere_direction;

      if (i < _inputPoints.size() - 2) {
        // For intermediate points, use the bisector of the two adjacent
        // segments
        Kernel::Vector_3 prev_dir = _inputPoints[i + 1] - _inputPoints[i];
        Kernel::Vector_3 next_dir = _inputPoints[i + 2] - _inputPoints[i + 1];
        prev_dir =
            prev_dir / CGAL::sqrt(CGAL::to_double(prev_dir.squared_length()));
        next_dir =
            next_dir / CGAL::sqrt(CGAL::to_double(next_dir.squared_length()));
        sphere_direction = prev_dir + next_dir;
      } else {
        // For the last point, use the direction of the last segment
        sphere_direction = _inputPoints[i + 1] - _inputPoints[i];
      }
      sphere_direction =
          sphere_direction /
          CGAL::sqrt(CGAL::to_double(sphere_direction.squared_length()));

      unsigned int subdivision_level =
          std::max(1U, static_cast<unsigned int>(_segments / 16));
      Sphere sphere(_radius, sphereCenter, subdivision_level, sphere_direction);
      CGAL::Polyhedron_3<Kernel> sphere_poly = sphere.generatePolyhedron();
      Nef_polyhedron             sphere_nef(sphere_poly);
      result = result.join(sphere_nef);
    }
  }

  // Convert the Nef_polyhedron to Polyhedron_3
  CGAL::Polyhedron_3<Kernel> merged_mesh;
  result.convert_to_polyhedron(merged_mesh);

  // Clean up the geometry
  PMP::remove_connected_components_of_negligible_size(merged_mesh);

  // Convert the merged mesh to PolyhedralSurface and return
  auto resultSurface = std::make_unique<PolyhedralSurface>();
  resultSurface->addPatchs(merged_mesh);

  return resultSurface;
}

auto
Buffer3D::computeFlatBuffer() const -> std::unique_ptr<PolyhedralSurface>
{
  std::vector<Kernel::Point_3> line_points;
  line_points.reserve(_inputPoints.size());
  for (const auto &p : _inputPoints) {
    line_points.emplace_back(p.x(), p.y(), p.z());
  }

  Surface_mesh_3                                         buffer;
  std::vector<std::vector<Surface_mesh_3::Vertex_index>> rings;

  std::vector<Kernel::Plane_3> bisector_planes;
  for (size_t i = 1; i < line_points.size() - 1; ++i) {
    bisector_planes.push_back(compute_bisector_plane(
        line_points[i - 1], line_points[i], line_points[i + 1]));
  }

  for (size_t i = 0; i < line_points.size() - 1; ++i) {
    Kernel::Vector_3 axis =
        normalizeVector(line_points[i + 1] - line_points[i]);
    Kernel::Point_3 extended_start =
        extend_point(line_points[i], -axis, _radius);
    Kernel::Point_3 extended_end =
        extend_point(line_points[i + 1], axis, _radius);

    std::vector<Kernel::Point_3> start_circle =
        create_circle_points(extended_start, axis, _radius, _segments);
    std::vector<Kernel::Point_3> end_circle =
        create_circle_points(extended_end, axis, _radius, _segments);

    std::vector<Surface_mesh_3::Vertex_index> start_ring;
    std::vector<Surface_mesh_3::Vertex_index> end_ring;

    for (int j = 0; j < _segments; ++j) {
      Kernel::Point_3 start_point = start_circle[j];
      Kernel::Point_3 end_point   = end_circle[j];

      if (i > 0) {
        start_point = intersect_segment_plane(start_point, end_point,
                                              bisector_planes[i - 1]);
      }
      if (i < line_points.size() - 2) {
        end_point =
            intersect_segment_plane(start_point, end_point, bisector_planes[i]);
      }

      Surface_mesh_3::Vertex_index v1 = buffer.add_vertex(start_point);
      Surface_mesh_3::Vertex_index v2 = buffer.add_vertex(end_point);
      start_ring.push_back(v1);
      end_ring.push_back(v2);

      if (j > 0) {
        buffer.add_face(start_ring[j - 1], start_ring[j], end_ring[j],
                        end_ring[j - 1]);
      }
    }
    buffer.add_face(start_ring.back(), start_ring.front(), end_ring.front(),
                    end_ring.back());

    rings.push_back(start_ring);
    if (i == line_points.size() - 2) {
      rings.push_back(end_ring);
    }
  }

  // Add caps
  Kernel::Point_3 start_center = extend_point(
      line_points.front(),
      normalizeVector(line_points[1] - line_points.front()), -_radius);
  Kernel::Point_3 end_center = extend_point(
      line_points.back(),
      normalizeVector(line_points.back() - line_points[line_points.size() - 2]),
      _radius);
  Surface_mesh_3::Vertex_index start_center_index =
      buffer.add_vertex(start_center);
  Surface_mesh_3::Vertex_index end_center_index = buffer.add_vertex(end_center);

  for (int i = 0; i < _segments; ++i) {
    buffer.add_face(start_center_index, rings.front()[(i + 1) % _segments],
                    rings.front()[i]);
    buffer.add_face(end_center_index, rings.back()[i],
                    rings.back()[(i + 1) % _segments]);
  }

  return std::make_unique<PolyhedralSurface>(buffer);
}

auto
Buffer3D::extend_point(const CGAL::Point_3<Kernel>  &point,
                       const CGAL::Vector_3<Kernel> &direction,
                       double distance) const -> CGAL::Point_3<Kernel>
{
  return point + direction * distance;
}

auto
Buffer3D::create_circle_points(const CGAL::Point_3<Kernel>  &center,
                               const CGAL::Vector_3<Kernel> &axis,
                               double radius, int segments) const
    -> std::vector<CGAL::Point_3<Kernel>>
{
  std::vector<CGAL::Point_3<Kernel>> points;
  CGAL::Vector_3<Kernel>             perpendicular =
      CGAL::cross_product(axis, CGAL::Vector_3<Kernel>(0, 0, 1));
  if (perpendicular == CGAL::NULL_VECTOR) {
    perpendicular = CGAL::cross_product(axis, Kernel::Vector_3(0, 1, 0));
  }
  perpendicular = normalizeVector(perpendicular);
  Kernel::Vector_3 perpendicular2 =
      normalizeVector(CGAL::cross_product(axis, perpendicular));

  for (int i = 0; i < segments; ++i) {
    double           angle  = 2.0 * M_PI * i / segments;
    Kernel::Vector_3 offset = radius * (std::cos(angle) * perpendicular +
                                        std::sin(angle) * perpendicular2);
    points.push_back(center + offset);
  }
  return points;
}

auto
Buffer3D::compute_bisector_plane(const CGAL::Point_3<Kernel> &p1,
                                 const CGAL::Point_3<Kernel> &p2,
                                 const CGAL::Point_3<Kernel> &p3) const
    -> CGAL::Plane_3<Kernel>
{
  CGAL::Vector_3<Kernel> v1       = normalizeVector(p2 - p1);
  CGAL::Vector_3<Kernel> v2       = normalizeVector(p3 - p2);
  CGAL::Vector_3<Kernel> bisector = v1 + v2;
  return {p2, bisector};
}

auto
Buffer3D::intersect_segment_plane(const CGAL::Point_3<Kernel> &p1,
                                  const CGAL::Point_3<Kernel> &p2,
                                  const CGAL::Plane_3<Kernel> &plane) const
    -> CGAL::Point_3<Kernel>
{
  CGAL::Vector_3<Kernel> v = p2 - p1;
  CGAL::Epeck::FT        t =
      -plane.a() * p1.x() - plane.b() * p1.y() - plane.c() * p1.z() - plane.d();
  t /= plane.a() * v.x() + plane.b() * v.y() + plane.c() * v.z();
  return p1 + t * v;
}

} // namespace SFCGAL::algorithm
