// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/GeometryVisitor.h"

using namespace SFCGAL::detail;

namespace SFCGAL {

///
///
///
PolyhedralSurface::PolyhedralSurface() = default;

///
///
///
PolyhedralSurface::PolyhedralSurface(const std::vector<Polygon> &polygons)

{
  for (const auto &polygon : polygons) {
    _polygons.push_back(polygon.clone());
  }
}

///
///
///
PolyhedralSurface::PolyhedralSurface(const PolyhedralSurface &other)

    = default;

///
///
///
PolyhedralSurface::PolyhedralSurface(const MarkedPolyhedron &poly)
{
  for (MarkedPolyhedron::Facet_const_iterator fit = poly.facets_begin();
       fit != poly.facets_end(); ++fit) {
    auto *face = new LineString();
    MarkedPolyhedron::Halfedge_around_facet_const_circulator hit =
        fit->facet_begin();

    do {
      face->addPoint(hit->vertex()->point());
      ++hit;
    } while (hit != fit->facet_begin());

    // close the ring
    face->addPoint(hit->vertex()->point());
    _polygons.push_back(new Polygon(face));
  }
}

///
///
///
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

///
///
///
auto
PolyhedralSurface::operator=(PolyhedralSurface other) -> PolyhedralSurface &
{
  swap(other);
  return *this;
}

///
///
///
PolyhedralSurface::~PolyhedralSurface() = default;

///
///
///
auto
PolyhedralSurface::clone() const -> PolyhedralSurface *
{
  return new PolyhedralSurface(*this);
}

///
///
///
auto
PolyhedralSurface::geometryType() const -> std::string
{
  return "PolyhedralSurface";
}

///
///
///
auto
PolyhedralSurface::geometryTypeId() const -> GeometryType
{
  return TYPE_POLYHEDRALSURFACE;
}

///
///
///
auto
PolyhedralSurface::dimension() const -> int
{
  return 2;
}

///
///
///
auto
PolyhedralSurface::coordinateDimension() const -> int
{
  if (isEmpty()) {
    return 0;
  }
  return _polygons.front().coordinateDimension();
}

///
///
///
auto
PolyhedralSurface::isEmpty() const -> bool
{
  return _polygons.empty();
}

///
///
///
auto
PolyhedralSurface::is3D() const -> bool
{
  if (isEmpty()) {
    return false;
  }
  return _polygons.front().is3D();
}

///
///
///
auto
PolyhedralSurface::isMeasured() const -> bool
{
  if (isEmpty()) {
    return false;
  }
  return _polygons.front().isMeasured();
}

///
///
///
auto
PolyhedralSurface::toTriangulatedSurface() const -> TriangulatedSurface
{
  TriangulatedSurface result;
  triangulate::triangulatePolygon3D(*this, result);
  return result;
}

///
///
///
void
PolyhedralSurface::addPolygon(const Polygon &polygon)
{
  addPolygon(polygon.clone());
}

///
///
///
void
PolyhedralSurface::addPolygon(Polygon *polygon)
{
  BOOST_ASSERT(polygon != NULL);
  _polygons.push_back(polygon);
}

///
///
///
void
PolyhedralSurface::addPolygons(const PolyhedralSurface &polyhedralSurface)
{
  for (size_t i = 0; i < polyhedralSurface.numPolygons(); i++) {
    addPolygon(polyhedralSurface.polygonN(i));
  }
}

///
///
///
auto
PolyhedralSurface::numGeometries() const -> size_t
{
  return _polygons.size();
}

///
///
///
auto
PolyhedralSurface::geometryN(size_t const &n) const -> const Polygon &
{
  return _polygons[n];
}

///
///
///
auto
PolyhedralSurface::geometryN(size_t const &n) -> Polygon &
{
  return _polygons[n];
}

///
///
///
void
PolyhedralSurface::accept(GeometryVisitor &visitor)
{
  return visitor.visit(*this);
}

///
///
///
void
PolyhedralSurface::accept(ConstGeometryVisitor &visitor) const
{
  return visitor.visit(*this);
}
} // namespace SFCGAL
