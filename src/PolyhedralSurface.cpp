// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/GeometryVisitor.h>
#include <SFCGAL/PolyhedralSurface.h>

using namespace SFCGAL::detail;

namespace SFCGAL {

///
///
///
PolyhedralSurface::PolyhedralSurface() : Surface(), _polygons() {}

///
///
///
PolyhedralSurface::PolyhedralSurface(const std::vector<Polygon> &polygons)
    : Surface()
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
PolyhedralSurface::PolyhedralSurface(const MarkedPolyhedron &poly) : Surface()
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
  } else {
    return _polygons.front().coordinateDimension();
  }
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
  } else {
    return _polygons.front().is3D();
  }
}

///
///
///
auto
PolyhedralSurface::isMeasured() const -> bool
{
  if (isEmpty()) {
    return false;
  } else {
    return _polygons.front().isMeasured();
  }
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
