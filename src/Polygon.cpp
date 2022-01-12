// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/GeometryVisitor.h>
#include <SFCGAL/Polygon.h>

#include <SFCGAL/Triangle.h>
#include <SFCGAL/algorithm/orientation.h>

namespace SFCGAL {

///
///
///
Polygon::Polygon() : Surface() { _rings.push_back(new LineString()); }

///
///
///
Polygon::Polygon(const std::vector<LineString> &rings) : Surface()
{
  if (rings.empty()) {
    _rings.resize(1, new LineString());
  } else {
    for (size_t i = 0; i < rings.size(); i++) {
      _rings.push_back(rings[i].clone());
    }
  }
}

///
///
///
Polygon::Polygon(const LineString &exteriorRing) : Surface()
{
  _rings.push_back(exteriorRing.clone());
}

///
///
///
Polygon::Polygon(LineString *exteriorRing) : Surface()
{
  _rings.push_back(exteriorRing);
}

///
///
///
Polygon::Polygon(const Triangle &triangle) : Surface()
{
  _rings.push_back(new LineString());

  if (!triangle.isEmpty()) {
    for (size_t i = 0; i < 4; i++) {
      exteriorRing().addPoint(triangle.vertex(i));
    }
  }
}

///
///
///
Polygon::Polygon(const Polygon &other) : Surface(other)
{
  for (size_t i = 0; i < other.numRings(); i++) {
    _rings.push_back(other.ringN(i).clone());
  }
}

///
///
///
Polygon::Polygon(const CGAL::Polygon_2<Kernel> &other)
{
  _rings.push_back(new LineString());
  CGAL::Polygon_2<Kernel>::Edge_const_iterator ei;

  for (ei = other.edges_begin(); ei != other.edges_end(); ++ei) {
    _rings.back().addPoint(ei->source());
  }
}

///
///
///
Polygon::Polygon(const CGAL::Polygon_with_holes_2<Kernel> &poly)
{
  _rings.push_back(new LineString());
  CGAL::Polygon_2<Kernel>                      outer = poly.outer_boundary();
  CGAL::Polygon_2<Kernel>::Edge_const_iterator ei;

  for (ei = outer.edges_begin(); ei != outer.edges_end(); ++ei) {
    _rings.back().addPoint(ei->source());
  }

  _rings.back().addPoint(_rings.back().startPoint());

  for (CGAL::Polygon_with_holes_2<Kernel>::Hole_const_iterator hit =
           poly.holes_begin();
       hit != poly.holes_end(); ++hit) {
    _rings.push_back(new LineString());
    CGAL::Polygon_2<Kernel>::Edge_const_iterator ei;

    for (ei = hit->edges_begin(); ei != hit->edges_end(); ++ei) {
      _rings.back().addPoint(ei->source());
    }

    _rings.back().addPoint(_rings.back().startPoint());
  }
}

///
///
///
Polygon &
Polygon::operator=(Polygon other)
{
  swap(other);
  return *this;
}

///
///
///
Polygon::~Polygon() {}

///
///
///
int
Polygon::coordinateDimension() const
{
  return _rings[0].coordinateDimension();
}

///
///
///
std::string
Polygon::geometryType() const
{
  return "Polygon";
}

///
///
///
GeometryType
Polygon::geometryTypeId() const
{
  return TYPE_POLYGON;
}

///
///
///
Polygon *
Polygon::clone() const
{
  return new Polygon(*this);
}

///
///
///
bool
Polygon::isEmpty() const
{
  return exteriorRing().isEmpty();
}

///
///
///
bool
Polygon::is3D() const
{
  return exteriorRing().is3D();
}

///
///
///
bool
Polygon::isMeasured() const
{
  return exteriorRing().isMeasured();
}

///
///
///
void
Polygon::reverse()
{
  for (size_t i = 0; i < numRings(); i++) {
    ringN(i).reverse();
  }
}

///
///
///
void
Polygon::accept(GeometryVisitor &visitor)
{
  return visitor.visit(*this);
}

///
///
///
void
Polygon::accept(ConstGeometryVisitor &visitor) const
{
  return visitor.visit(*this);
}

///
///
///
bool
Polygon::isCounterClockWiseOriented() const
{
  return algorithm::isCounterClockWiseOriented(*this);
}

///
///
///
CGAL::Polygon_2<Kernel>
Polygon::toPolygon_2(bool fixOrientation) const
{
  return exteriorRing().toPolygon_2(fixOrientation);
}

///
///
///
CGAL::Polygon_with_holes_2<Kernel>
Polygon::toPolygon_with_holes_2(bool fixOrientation) const
{
  std::list<CGAL::Polygon_2<Kernel>> holes;

  for (size_t i = 0; i < numInteriorRings(); ++i) {
    // note that the orientation is fixed here to avoid double reverse for
    // interior rings
    CGAL::Polygon_2<Kernel> inner = interiorRingN(i).toPolygon_2(false);

    if (fixOrientation && inner.orientation() == CGAL::COUNTERCLOCKWISE) {
      inner.reverse_orientation();
    }

    holes.push_back(inner);
  }

  CGAL::Polygon_2<Kernel> outer = exteriorRing().toPolygon_2(fixOrientation);
  return CGAL::Polygon_with_holes_2<Kernel>(outer, holes.begin(), holes.end());
}

} // namespace SFCGAL
