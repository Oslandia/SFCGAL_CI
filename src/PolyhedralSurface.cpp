// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryVisitor.h"

#include <CGAL/Cartesian_converter.h>

using namespace SFCGAL::detail;

namespace SFCGAL {

PolyhedralSurface::PolyhedralSurface() = default;

PolyhedralSurface::PolyhedralSurface(const std::vector<Polygon> &polygons)

{
  for (const auto &polygon : polygons) {
    _polygons.push_back(polygon.clone());
  }
}

PolyhedralSurface::PolyhedralSurface(const std::unique_ptr<Geometry> &geometry)
{
  if (geometry->is<PolyhedralSurface>()) {
    *this = static_cast<const PolyhedralSurface &>(*geometry);
  } else if (geometry->is<TriangulatedSurface>()) {
    const TriangulatedSurface &triangulatedSurface =
        geometry->as<TriangulatedSurface>();
    for (size_t i = 0; i < triangulatedSurface.numPatches(); ++i) {
      this->addPatch(triangulatedSurface.patchN(i));
    }
  } else if (geometry->is<Polygon>()) {
    this->addPatch(geometry->as<Polygon>());
  } else {
    throw std::invalid_argument("Cannot convert geometry to PolyhedralSurface");
  }
}

PolyhedralSurface::PolyhedralSurface(const PolyhedralSurface &other)

    = default;

PolyhedralSurface::PolyhedralSurface(const Mesh &sm)
{

  using vertex_descriptor = Mesh::Vertex_index;
  for (auto face : sm.faces()) {
    auto *new_face = new LineString();
    for (vertex_descriptor const vd :
         vertices_around_face(sm.halfedge(face), sm)) {
      new_face->addPoint(Point(sm.point(vd)));
    }

    new_face->addPoint(new_face->startPoint().clone());
    _polygons.push_back(new Polygon(new_face));
  }
}

PolyhedralSurface::PolyhedralSurface(const InexactMesh &inexactMesh)
{
  using inexact_to_exact  = CGAL::Cartesian_converter<InexactKernel, Kernel>;
  using vertex_descriptor = Mesh::Vertex_index;

  inexact_to_exact toExact;
  for (auto face : inexactMesh.faces()) {
    auto *new_face = new LineString();
    for (vertex_descriptor const vd :
         vertices_around_face(inexactMesh.halfedge(face), inexactMesh)) {
      new_face->addPoint(Point(toExact(inexactMesh.point(vd))));
    }

    new_face->addPoint(new_face->startPoint().clone());
    _polygons.push_back(new Polygon(new_face));
  }
}

auto
PolyhedralSurface::operator=(PolyhedralSurface other) -> PolyhedralSurface &
{
  swap(other);
  return *this;
}

PolyhedralSurface::~PolyhedralSurface() = default;

auto
PolyhedralSurface::clone() const -> PolyhedralSurface *
{
  return new PolyhedralSurface(*this);
}

auto
PolyhedralSurface::geometryType() const -> std::string
{
  return "PolyhedralSurface";
}

auto
PolyhedralSurface::geometryTypeId() const -> GeometryType
{
  return TYPE_POLYHEDRALSURFACE;
}

auto
PolyhedralSurface::dimension() const -> int
{
  return 2;
}

auto
PolyhedralSurface::coordinateDimension() const -> int
{
  if (isEmpty()) {
    return 0;
  }
  return _polygons.front().coordinateDimension();
}

auto
PolyhedralSurface::isEmpty() const -> bool
{
  return _polygons.empty();
}

auto
PolyhedralSurface::is3D() const -> bool
{
  if (isEmpty()) {
    return false;
  }
  return _polygons.front().is3D();
}

auto
PolyhedralSurface::isMeasured() const -> bool
{
  if (isEmpty()) {
    return false;
  }
  return _polygons.front().isMeasured();
}

auto
PolyhedralSurface::dropZ() -> bool
{
  if (!is3D()) {
    return false;
  }

  for (auto &_polygon : _polygons) {
    _polygon.dropZ();
  }

  return true;
}

auto
PolyhedralSurface::dropM() -> bool
{
  if (!isMeasured()) {
    return false;
  }

  for (auto &_polygon : _polygons) {
    _polygon.dropM();
  }

  return true;
}

auto
PolyhedralSurface::swapXY() -> void
{
  for (auto &_polygon : _polygons) {
    _polygon.swapXY();
  }
}

auto
PolyhedralSurface::toTriangulatedSurface() const -> TriangulatedSurface
{
  TriangulatedSurface result;
  triangulate::triangulatePolygon3D(*this, result);
  return result;
}

void
PolyhedralSurface::addPatch(const Polygon &patch)
{
  addPatch(patch.clone());
}

void
PolyhedralSurface::addPatch(Polygon *patch)
{
  BOOST_ASSERT(patch != NULL);
  _polygons.push_back(patch);
}

void
PolyhedralSurface::addPatchs(const PolyhedralSurface &polyhedralSurface)
{
  for (size_t i = 0; i < polyhedralSurface.numPatches(); i++) {
    addPatch(polyhedralSurface.patchN(i));
  }
}

void
PolyhedralSurface::addPolygon(const Polygon &polygon)
{
  return addPatch(polygon);
}

void
PolyhedralSurface::addPolygon(Polygon *polygon)
{
  return addPatch(polygon);
}

void
PolyhedralSurface::addPolygons(const PolyhedralSurface &polyhedralSurface)
{
  return addPatchs(polyhedralSurface);
}

void
PolyhedralSurface::setPatchN(Polygon *patch, size_t const &n)
{
  BOOST_ASSERT(patch != NULL);

  if (n >= numPatches()) {
    BOOST_THROW_EXCEPTION(
        Exception((boost::format("Cannot set geometry at position %s. "
                                 "PolyhedralSurface has only %d geometries.") %
                   n % numPatches())
                      .str()));
  }

  _polygons.replace(n, patch);
}

void
PolyhedralSurface::setPatchN(const Polygon &patch, size_t const &n)
{
  setPatchN(patch.clone(), n);
}

void
PolyhedralSurface::setPatchN(Geometry *geometry, size_t const &n)
{
  if (geometry->geometryTypeId() != TYPE_POLYGON) {
    std::ostringstream oss;
    oss << "try to set a '" << geometry->geometryType()
        << "' in a PolyhedralSurface\n";
    delete geometry; // we are responsible for the resource here
    BOOST_THROW_EXCEPTION(InappropriateGeometryException(oss.str()));
  }

  setPatchN(dynamic_cast<Polygon *>(geometry), n);
}

void
PolyhedralSurface::setPatchN(const Geometry &geometry, size_t const &n)
{
  setPatchN(geometry.clone(), n);
}

void
PolyhedralSurface::accept(GeometryVisitor &visitor)
{
  return visitor.visit(*this);
}

void
PolyhedralSurface::accept(ConstGeometryVisitor &visitor) const
{
  return visitor.visit(*this);
}

// Explicit instantiations
template PolyhedralSurface::PolyhedralSurface(const detail::MarkedPolyhedron &);
template PolyhedralSurface::PolyhedralSurface(
    const CGAL::Polyhedron_3<Kernel> &);

} // namespace SFCGAL
