// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/format.hpp>

#include "SFCGAL/Exception.h"
#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryVisitor.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/TriangulatedSurface.h"

#include <map>

namespace SFCGAL {

TriangulatedSurface::TriangulatedSurface() = default;

TriangulatedSurface::TriangulatedSurface(const std::vector<Triangle> &triangles)
{
  _triangles.reserve(triangles.size());
  for (const auto &triangle : triangles) {
    _triangles.emplace_back(triangle.clone());
  }
}

TriangulatedSurface::TriangulatedSurface(const TriangulatedSurface &other)
    : GeometryImpl(other)
{
  _triangles.reserve(other._triangles.size());
  for (const auto &triangle : other._triangles) {
    _triangles.emplace_back(triangle->clone());
  }
}

auto
TriangulatedSurface::operator=(TriangulatedSurface other)
    -> TriangulatedSurface &
{
  swap(other);
  return *this;
}

TriangulatedSurface::~TriangulatedSurface() = default;

auto
TriangulatedSurface::geometryType() const -> std::string
{
  return "TriangulatedSurface";
}

auto
TriangulatedSurface::geometryTypeId() const -> GeometryType
{
  return TYPE_TRIANGULATEDSURFACE;
}

auto
TriangulatedSurface::dimension() const -> int
{
  // surface
  return 2;
}

auto
TriangulatedSurface::coordinateDimension() const -> int
{
  if (_triangles.empty()) {
    return 0;
  }
  return _triangles[0]->coordinateDimension();
}

auto
TriangulatedSurface::isEmpty() const -> bool
{
  return _triangles.empty();
}

auto
TriangulatedSurface::is3D() const -> bool
{
  return !_triangles.empty() && _triangles.front()->is3D();
}

auto
TriangulatedSurface::isMeasured() const -> bool
{
  return !_triangles.empty() && _triangles.front()->isMeasured();
}

auto
TriangulatedSurface::dropZ() -> bool
{
  if (!is3D()) {
    return false;
  }

  for (auto &_triangle : _triangles) {
    _triangle->dropZ();
  }

  return true;
}

auto
TriangulatedSurface::dropM() -> bool
{
  if (!isMeasured()) {
    return false;
  }

  for (auto &_triangle : _triangles) {
    _triangle->dropM();
  }

  return true;
}

auto
TriangulatedSurface::swapXY() -> void
{
  for (auto &_triangle : _triangles) {
    _triangle->swapXY();
  }
}

void
TriangulatedSurface::addPatchs(const TriangulatedSurface &other)
{
  for (const auto &it : other) {
    addPatch(it);
  }
}

void
TriangulatedSurface::addTriangles(const TriangulatedSurface &other)
{
  addPatchs(other);
}

auto
TriangulatedSurface::patchN(size_t const &index) const -> const Triangle &
{
  if (index >= numPatches()) {
    BOOST_THROW_EXCEPTION(Exception(
        (boost::format("Cannot access geometry at position %s. "
                       "TriangulatedSurface has only %d geometries.") %
         index % numGeometries())
            .str()));
  }

  return *_triangles[index];
}

auto
TriangulatedSurface::patchN(size_t const &index) -> Triangle &
{
  if (index >= numPatches()) {
    BOOST_THROW_EXCEPTION(Exception(
        (boost::format("Cannot access geometry at position %s. "
                       "TriangulatedSurface has only %d geometries.") %
         index % numPatches())
            .str()));
  }

  return *_triangles[index];
}

void
TriangulatedSurface::setPatchN(std::unique_ptr<Triangle> triangle,
                               size_t const             &idx)
{
  BOOST_ASSERT(triangle != nullptr);

  if (idx >= numPatches()) {
    BOOST_THROW_EXCEPTION(Exception(
        (boost::format("Cannot set geometry at position %s. "
                       "TriangulatedSurface has only %d geometries.") %
         idx % numPatches())
            .str()));
  }

  _triangles[idx] = std::move(triangle);
}

void
TriangulatedSurface::setPatchN(Triangle *triangle, size_t const &idx)
{
  setPatchN(std::unique_ptr<Triangle>(triangle), idx);
}

void
TriangulatedSurface::setPatchN(const Triangle &triangle, size_t const &idx)
{
  setPatchN(triangle.clone(), idx);
}

void
TriangulatedSurface::setPatchN(std::unique_ptr<Geometry> geometry,
                               size_t const             &idx)
{
  if (geometry->geometryTypeId() != TYPE_TRIANGLE) {
    std::ostringstream oss;
    oss << "try to set a '" << geometry->geometryType()
        << "' in a TriangulatedSurface\n";
    BOOST_THROW_EXCEPTION(InappropriateGeometryException(oss.str()));
  }

  setPatchN(geom_unique_ptr_as<Triangle>(std::move(geometry)), idx);
}

void
TriangulatedSurface::setPatchN(Geometry *geometry, size_t const &idx)
{
  setPatchN(std::unique_ptr<Geometry>(geometry), idx);
}

void
TriangulatedSurface::setPatchN(const Geometry &geometry, size_t const &idx)
{
  setPatchN(geometry.clone(), idx);
}

void
TriangulatedSurface::reserve(const size_t &count)
{
  _triangles.reserve(count);
}

void
TriangulatedSurface::accept(GeometryVisitor &visitor)
{
  visitor.visit(*this);
}

void
TriangulatedSurface::accept(ConstGeometryVisitor &visitor) const
{
  visitor.visit(*this);
}

/// @{
/// @privatesection
#ifndef DOXYGEN_SHOULD_SKIP_THIS

// Private class
// A modifier creating triangles from a TriangulatedSurface with the incremental
// builder.
template <class HDS>
class Triangulated2Polyhedron : public CGAL::Modifier_base<HDS> {
public:
  Triangulated2Polyhedron(const TriangulatedSurface &surf) : surf(surf) {}

  using Vertex      = typename HDS::Vertex;
  using Point       = typename Vertex::Point;
  using K           = typename HDS::Traits;
  using PointMap    = std::map<Point, size_t>;
  using HalfedgeSet = std::set<std::pair<Point, Point>>;

  void
  operator()(HDS &hds) override
  {
    const size_t nrPatchs = surf.numPatches();
    // Postcondition: `hds' is a valid polyhedral surface.
    CGAL::Polyhedron_incremental_builder_3<HDS> builder(hds, true);
    builder.begin_surface(/* vertices */ nrPatchs * 3,
                          /* facets */ nrPatchs,
                          /* halfedges */ nrPatchs * 3);

    size_t vertex_idx = 0;

    // first pass: insert vertices, only if they are not shared between faces
    // thanks to a binary tree (PointMap)
    for (size_t i = 0; i < nrPatchs; i++) {
      for (size_t j = 0; j < 3; j++) {
        Point const point =
            surf.patchN(i).vertex(static_cast<int>(j)).toPoint_3();

        if (points.find(point) == points.end()) {
          builder.add_vertex(point);
          points[point] = vertex_idx++;
        }
      }
    }

    // second pass: adjacent triangles must be built with compliant orientations
    // the two halfedges of a shared edge must be of opposite orientation

    // Extract from CGAL's documentation
    // "The convention is that the halfedges are oriented counterclockwise
    // around facets as seen from the outside of the polyhedron"

    for (size_t i = 0; i < nrPatchs; i++) {
      builder.begin_facet();
      CGAL::Triangle_3<K> const tri(surf.patchN(i).toTriangle_3());
      CGAL::Point_3<K> const    pointA(tri[0]);
      CGAL::Point_3<K> const    pointB(tri[1]);
      CGAL::Point_3<K> const    pointC(tri[2]);

      if (edges.find(std::make_pair(pointA, pointB)) != edges.end() ||
          edges.find(std::make_pair(pointB, pointC)) != edges.end() ||
          edges.find(std::make_pair(pointC, pointA)) != edges.end()) {
        BOOST_THROW_EXCEPTION(
            Exception("When trying to build a CGAL::Polyhedron_3 from a "
                      "TriangulatedSurface: bad orientation for " +
                      surf.geometryN(i).asText() +
                      " consider using ConsistentOrientationBuilder first"));
      }

      builder.add_vertex_to_facet(points[pointA]);
      builder.add_vertex_to_facet(points[pointB]);
      builder.add_vertex_to_facet(points[pointC]);
      edges.insert(std::make_pair(pointA, pointB));
      edges.insert(std::make_pair(pointB, pointC));
      edges.insert(std::make_pair(pointC, pointA));
      builder.end_facet();
    }

    builder.end_surface();
  }

private:
  const TriangulatedSurface &surf;
  PointMap                   points;
  HalfedgeSet                edges;
};

template <typename Polyhedron>
struct Plane_from_facet {
  auto
  operator()(typename Polyhedron::Facet &facet) -> typename Polyhedron::Plane_3
  {
    typename Polyhedron::Halfedge_handle const halfedge = facet.halfedge();
    return typename Polyhedron::Plane_3(
        halfedge->vertex()->point(), halfedge->next()->vertex()->point(),
        halfedge->opposite()->vertex()->point());
  }
};

#endif // ifndef DOXYGEN_SHOULD_SKIP_THIS
/// @} end of private section

template <typename Polyhedron>
auto
TriangulatedSurface::toPolyhedron_3() const -> std::unique_ptr<Polyhedron>
{
  auto *poly = new Polyhedron();
  Triangulated2Polyhedron<typename Polyhedron::HalfedgeDS> converter(*this);
  poly->delegate(converter);

  // compute planes
  std::transform(poly->facets_begin(), poly->facets_end(), poly->planes_begin(),
                 Plane_from_facet<Polyhedron>());

  return std::unique_ptr<Polyhedron>(poly);
}

auto
TriangulatedSurface::toSurfaceMesh() const -> Surface_mesh_3
{
  Surface_mesh_3 mesh;

  if (isEmpty()) {
    return mesh;
  }

  // Map to track unique vertices
  std::map<Point, Surface_mesh_3::Vertex_index> vertexMap;

  // Add all triangles to the mesh
  for (size_t i = 0; i < numPatches(); ++i) {
    const Triangle                             &triangle = patchN(i);
    std::array<Surface_mesh_3::Vertex_index, 3> vertices;

    for (size_t j = 0; j < 3; ++j) {
      const Point &pt = triangle.vertex(static_cast<int>(j));

      // Check if vertex already exists
      auto it = vertexMap.find(pt);
      if (it != vertexMap.end()) {
        vertices[j] = it->second;
      } else {
        // Add new vertex
        Surface_mesh_3::Vertex_index vIdx = mesh.add_vertex(pt.toPoint_3());
        vertexMap[pt]                     = vIdx;
        vertices[j]                       = vIdx;
      }
    }

    // Add face to the mesh
    mesh.add_face(vertices[0], vertices[1], vertices[2]);
  }

  return mesh;
}

/**
 * @brief Explicit template instantiation of toPolyhedron_3 for
 * MarkedPolyhedron.
 *
 * Converts the triangulated surface to a MarkedPolyhedron.
 * Returns a std::unique_ptr owning the polyhedron.
 */
template SFCGAL_API std::unique_ptr<detail::MarkedPolyhedron>
TriangulatedSurface::toPolyhedron_3<detail::MarkedPolyhedron>() const;

/**
 * @brief Explicit template instantiation of toPolyhedron_3 for
 * CGAL::Polyhedron_3.
 *
 * Converts the triangulated surface to a CGAL Polyhedron_3.
 * Returns a std::unique_ptr owning the polyhedron.
 */
template SFCGAL_API std::unique_ptr<CGAL::Polyhedron_3<Kernel>>
TriangulatedSurface::toPolyhedron_3<CGAL::Polyhedron_3<Kernel>>() const;

} // namespace SFCGAL
