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
  for (size_t i = 0; i < polygons.size(); i++) {
    _polygons.push_back(polygons[i].clone());
  }
}

///
///
///
PolyhedralSurface::PolyhedralSurface(const PolyhedralSurface &other)
    : Surface(other), _polygons(other._polygons)
{
}

///
///
///
PolyhedralSurface::PolyhedralSurface(const MarkedPolyhedron &poly) : Surface()
{
  for (MarkedPolyhedron::Facet_const_iterator fit = poly.facets_begin();
       fit != poly.facets_end(); ++fit) {
    LineString *face = new LineString();
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
PolyhedralSurface &
PolyhedralSurface::operator=(PolyhedralSurface other)
{
  swap(other);
  return *this;
}

///
///
///
PolyhedralSurface::~PolyhedralSurface() {}

///
///
///
PolyhedralSurface *
PolyhedralSurface::clone() const
{
  return new PolyhedralSurface(*this);
}

///
///
///
std::string
PolyhedralSurface::geometryType() const
{
  return "PolyhedralSurface";
}

///
///
///
GeometryType
PolyhedralSurface::geometryTypeId() const
{
  return TYPE_POLYHEDRALSURFACE;
}

///
///
///
int
PolyhedralSurface::dimension() const
{
  return 2;
}

///
///
///
int
PolyhedralSurface::coordinateDimension() const
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
bool
PolyhedralSurface::isEmpty() const
{
  return _polygons.empty();
}

///
///
///
bool
PolyhedralSurface::is3D() const
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
bool
PolyhedralSurface::isMeasured() const
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
TriangulatedSurface
PolyhedralSurface::toTriangulatedSurface() const
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
size_t
PolyhedralSurface::numGeometries() const
{
  return _polygons.size();
}

///
///
///
const Polygon &
PolyhedralSurface::geometryN(size_t const &n) const
{
  return _polygons[n];
}

///
///
///
Polygon &
PolyhedralSurface::geometryN(size_t const &n)
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
