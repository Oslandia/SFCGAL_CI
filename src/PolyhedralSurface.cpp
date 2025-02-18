// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
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
    for (size_t i = 0; i < triangulatedSurface.numTriangles(); ++i) {
      this->addPolygon(triangulatedSurface.triangleN(i));
    }
  } else if (geometry->is<Polygon>()) {
    this->addPolygon(geometry->as<Polygon>());
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
PolyhedralSurface::toTriangulatedSurface() const -> TriangulatedSurface
{
  TriangulatedSurface result;
  triangulate::triangulatePolygon3D(*this, result);
  return result;
}

void
PolyhedralSurface::addPolygon(const Polygon &polygon)
{
  addPolygon(polygon.clone());
}

void
PolyhedralSurface::addPolygon(Polygon *polygon)
{
  BOOST_ASSERT(polygon != NULL);
  _polygons.push_back(polygon);
}

void
PolyhedralSurface::addPolygons(const PolyhedralSurface &polyhedralSurface)
{
  for (size_t i = 0; i < polyhedralSurface.numPolygons(); i++) {
    addPolygon(polyhedralSurface.polygonN(i));
  }
}

auto
PolyhedralSurface::numGeometries() const -> size_t
{
  return _polygons.size();
}

auto
PolyhedralSurface::geometryN(size_t const &n) const -> const Polygon &
{
  if (n >= numGeometries()) {
    BOOST_THROW_EXCEPTION(
        Exception((boost::format("Cannot access geometry at position %s. "
                                 "PolyhedralSurface has only %d geometries.") %
                   n % numGeometries())
                      .str()));
  }

  return _polygons[n];
}

auto
PolyhedralSurface::geometryN(size_t const &n) -> Polygon &
{
  if (n >= numGeometries()) {
    BOOST_THROW_EXCEPTION(
        Exception((boost::format("Cannot access geometry at position %s. "
                                 "PolyhedralSurface has only %d geometries.") %
                   n % numGeometries())
                      .str()));
  }

  return _polygons[n];
}

void
PolyhedralSurface::setGeometryN(Polygon *polygon, size_t const &n)
{
  BOOST_ASSERT(polygon != NULL);

  if (n >= numGeometries()) {
    BOOST_THROW_EXCEPTION(
        Exception((boost::format("Cannot set geometry at position %s. "
                                 "PolyhedralSurface has only %d geometries.") %
                   n % numGeometries())
                      .str()));
  }

  _polygons.replace(n, polygon);
}

void
PolyhedralSurface::setGeometryN(const Polygon &polygon, size_t const &n)
{
  setGeometryN(polygon.clone(), n);
}

void
PolyhedralSurface::setGeometryN(Geometry *geometry, size_t const &n)
{
  if (geometry->geometryTypeId() != TYPE_POLYGON) {
    std::ostringstream oss;
    oss << "try to set a '" << geometry->geometryType()
        << "' in a PolyhedralSurface\n";
    delete geometry; // we are responsible for the resource here
    BOOST_THROW_EXCEPTION(InappropriateGeometryException(oss.str()));
  }

  setGeometryN(dynamic_cast<Polygon *>(geometry), n);
}

void
PolyhedralSurface::setGeometryN(const Geometry &geometry, size_t const &n)
{
  setGeometryN(geometry.clone(), n);
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
