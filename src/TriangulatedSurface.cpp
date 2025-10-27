// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/format.hpp>

#include "SFCGAL/Exception.h"
#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryVisitor.h"
#include "SFCGAL/TriangulatedSurface.h"

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
TriangulatedSurface::patchN(size_t const &n) const -> const Triangle &
{
  if (n >= numPatches()) {
    BOOST_THROW_EXCEPTION(Exception(
        (boost::format("Cannot access geometry at position %s. "
                       "TriangulatedSurface has only %d geometries.") %
         n % numGeometries())
            .str()));
  }

  return *_triangles[n];
}

auto
TriangulatedSurface::patchN(size_t const &n) -> Triangle &
{
  if (n >= numPatches()) {
    BOOST_THROW_EXCEPTION(Exception(
        (boost::format("Cannot access geometry at position %s. "
                       "TriangulatedSurface has only %d geometries.") %
         n % numPatches())
            .str()));
  }

  return *_triangles[n];
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
TriangulatedSurface::reserve(const size_t &n)
{
  _triangles.reserve(n);
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
    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
    B.begin_surface(/* vertices */ nrPatchs * 3,
                    /* facets */ nrPatchs,
                    /* halfedges */ nrPatchs * 3);

    size_t vertex_idx = 0;

    // first pass: insert vertices, only if they are not shared between faces
    // thanks to a binary tree (PointMap)
    for (size_t i = 0; i < nrPatchs; i++) {
      for (size_t j = 0; j < 3; j++) {
        Point const p = surf.patchN(i).vertex(j).toPoint_3();

        if (points.find(p) == points.end()) {
          B.add_vertex(p);
          points[p] = vertex_idx++;
        }
      }
    }

    // second pass: adjacent triangles must be built with compliant orientations
    // the two halfedges of a shared edge must be of opposite orientation

    // Extract from CGAL's documentation
    // "The convention is that the halfedges are oriented counterclockwise
    // around facets as seen from the outside of the polyhedron"

    for (size_t i = 0; i < nrPatchs; i++) {
      B.begin_facet();
      CGAL::Triangle_3<K> const tri(surf.patchN(i).toTriangle_3());
      CGAL::Point_3<K> const    pa(tri[0]);
      CGAL::Point_3<K> const    pb(tri[1]);
      CGAL::Point_3<K> const    pc(tri[2]);

      if (edges.find(std::make_pair(pa, pb)) != edges.end() ||
          edges.find(std::make_pair(pb, pc)) != edges.end() ||
          edges.find(std::make_pair(pc, pa)) != edges.end()) {
        BOOST_THROW_EXCEPTION(
            Exception("When trying to build a CGAL::Polyhedron_3 from a "
                      "TriangulatedSurface: bad orientation for " +
                      surf.geometryN(i).asText() +
                      " consider using ConsistentOrientationBuilder first"));
      }

      B.add_vertex_to_facet(points[pa]);
      B.add_vertex_to_facet(points[pb]);
      B.add_vertex_to_facet(points[pc]);
      edges.insert(std::make_pair(pa, pb));
      edges.insert(std::make_pair(pb, pc));
      edges.insert(std::make_pair(pc, pa));
      B.end_facet();
    }

    B.end_surface();
  }

private:
  const TriangulatedSurface &surf;
  PointMap                   points;
  HalfedgeSet                edges;
};

template <typename Polyhedron>
struct Plane_from_facet {
  auto
  operator()(typename Polyhedron::Facet &f) -> typename Polyhedron::Plane_3
  {
    typename Polyhedron::Halfedge_handle const h = f.halfedge();
    return typename Polyhedron::Plane_3(h->vertex()->point(),
                                        h->next()->vertex()->point(),
                                        h->opposite()->vertex()->point());
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
