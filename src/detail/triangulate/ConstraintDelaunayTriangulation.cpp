// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/triangulate/ConstraintDelaunayTriangulation.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/TriangulatedSurface.h"

#include "SFCGAL/detail/triangulate/markDomains.h"

namespace SFCGAL::triangulate {

ConstraintDelaunayTriangulation::ConstraintDelaunayTriangulation() = default;

auto
ConstraintDelaunayTriangulation::addVertex(const Coordinate &position)
    -> ConstraintDelaunayTriangulation::Vertex_handle
{
  if (position.isEmpty()) {
    BOOST_THROW_EXCEPTION(Exception(
        "try to add empty position to ConstraintDelaunayTriangulation"));
  }

  Vertex_handle vertex =
      _projectionPlane
          ? _cdt.insert(_projectionPlane->to_2d(position.toPoint_3()))
          : _cdt.insert(position.toPoint_2());
  vertex->info().original = position;
  return vertex;
}

void
ConstraintDelaunayTriangulation::addConstraint(Vertex_handle source,
                                               Vertex_handle target)
{
  if (source == target) {
    return;
  }

  _cdt.insert_constraint(source, target);
}

void
ConstraintDelaunayTriangulation::clear()
{
  _cdt.clear();
}

auto
ConstraintDelaunayTriangulation::numVertices() const -> size_t
{
  return _cdt.number_of_vertices();
}

auto
ConstraintDelaunayTriangulation::numTriangles() const -> size_t
{
  return _cdt.number_of_faces();
}

void
ConstraintDelaunayTriangulation::setProjectionPlane(
    const Kernel::Plane_3 &projectionPlane)
{
  BOOST_ASSERT(!projectionPlane.is_degenerate());
  _projectionPlane = projectionPlane;
}

auto
ConstraintDelaunayTriangulation::projectionPlane() const -> Kernel::Plane_3
{
  if (_projectionPlane) {
    return *_projectionPlane;
  }
  return Kernel::Plane_3(Kernel::RT(0), Kernel::RT(0), Kernel::RT(1),
                         Kernel::RT(0));
}

/// adapted from CGAL example
void
ConstraintDelaunayTriangulation::markDomains()
{
  detail::markDomains(_cdt);
}

void
ConstraintDelaunayTriangulation::getTriangles(
    TriangulatedSurface &triangulatedSurface, bool filterExteriorParts) const
{
  triangulatedSurface.reserve(triangulatedSurface.numPatches() +
                              numTriangles());

  for (Finite_faces_iterator it = finite_faces_begin();
       it != finite_faces_end(); ++it) {
    if (filterExteriorParts && (it->info().nestingLevel % 2 == 0)) {
      continue;
    }

    const Coordinate &a = it->vertex(0)->info().original;

    const Coordinate &b = it->vertex(1)->info().original;

    const Coordinate &c = it->vertex(2)->info().original;

    if (!a.isEmpty() && !b.isEmpty() && !c.isEmpty()) {
      triangulatedSurface.addPatch(new Triangle(Point(a), Point(b), Point(c)));
    }
  }
}

auto
ConstraintDelaunayTriangulation::getTriangulatedSurface() const
    -> std::unique_ptr<TriangulatedSurface>
{
  std::unique_ptr<TriangulatedSurface> result(new TriangulatedSurface);
  getTriangles(*result, false);
  return result;
}

} // namespace SFCGAL::triangulate
