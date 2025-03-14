// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/Polygon.h"
#include "SFCGAL/GeometryVisitor.h"

#include "SFCGAL/Triangle.h"
#include "SFCGAL/algorithm/orientation.h"

namespace SFCGAL {

Polygon::Polygon() { _rings.push_back(new LineString()); }

Polygon::Polygon(const std::vector<LineString> &rings)
{
  if (rings.empty()) {
    _rings.resize(1, new LineString());
  } else {
    for (const auto &ring : rings) {
      _rings.push_back(ring.clone());
    }
  }
}

Polygon::Polygon(const LineString &exteriorRing)
{
  _rings.push_back(exteriorRing.clone());
}

Polygon::Polygon(LineString *exteriorRing) { _rings.push_back(exteriorRing); }

Polygon::Polygon(const Triangle &triangle)
{
  _rings.push_back(new LineString());

  if (!triangle.isEmpty()) {
    for (size_t i = 0; i < 4; i++) {
      exteriorRing().addPoint(triangle.vertex(i));
    }
  }
}

Polygon::Polygon(const Polygon &other) : Surface(other)
{
  for (size_t i = 0; i < other.numRings(); i++) {
    _rings.push_back(other.ringN(i).clone());
  }
}

Polygon::Polygon(const CGAL::Polygon_2<Kernel> &other)
{
  _rings.push_back(new LineString());
  CGAL::Polygon_2<Kernel>::Edge_const_iterator ei;

  for (ei = other.edges_begin(); ei != other.edges_end(); ++ei) {
    _rings.back().addPoint(ei->source());
  }
}

Polygon::Polygon(const CGAL::Polygon_with_holes_2<Kernel> &poly)
{
  _rings.push_back(new LineString());
  const CGAL::Polygon_2<Kernel>               &outer = poly.outer_boundary();
  CGAL::Polygon_2<Kernel>::Edge_const_iterator ei;

  for (ei = outer.edges_begin(); ei != outer.edges_end(); ++ei) {
    _rings.back().addPoint(ei->source());
  }

  _rings.back().addPoint(_rings.back().startPoint());

  for (auto hit = poly.holes_begin(); hit != poly.holes_end(); ++hit) {
    _rings.push_back(new LineString());
    CGAL::Polygon_2<Kernel>::Edge_const_iterator ei;

    for (ei = hit->edges_begin(); ei != hit->edges_end(); ++ei) {
      _rings.back().addPoint(ei->source());
    }

    _rings.back().addPoint(_rings.back().startPoint());
  }
}

auto
Polygon::operator=(Polygon other) -> Polygon &
{
  swap(other);
  return *this;
}

Polygon::~Polygon() = default;

auto
Polygon::coordinateDimension() const -> int
{
  return _rings[0].coordinateDimension();
}

auto
Polygon::geometryType() const -> std::string
{
  return "Polygon";
}

auto
Polygon::geometryTypeId() const -> GeometryType
{
  return TYPE_POLYGON;
}

auto
Polygon::clone() const -> Polygon *
{
  return new Polygon(*this);
}

auto
Polygon::isEmpty() const -> bool
{
  return exteriorRing().isEmpty();
}

auto
Polygon::is3D() const -> bool
{
  return exteriorRing().is3D();
}

auto
Polygon::isMeasured() const -> bool
{
  return exteriorRing().isMeasured();
}

void
Polygon::reverse()
{
  for (size_t i = 0; i < numRings(); i++) {
    ringN(i).reverse();
  }
}

void
Polygon::accept(GeometryVisitor &visitor)
{
  return visitor.visit(*this);
}

void
Polygon::accept(ConstGeometryVisitor &visitor) const
{
  return visitor.visit(*this);
}

auto
Polygon::isCounterClockWiseOriented() const -> bool
{
  return algorithm::isCounterClockWiseOriented(*this);
}

auto
Polygon::toPolygon_2(bool fixOrientation) const -> CGAL::Polygon_2<Kernel>
{
  return exteriorRing().toPolygon_2(fixOrientation);
}

auto
Polygon::toPolygon_with_holes_2(bool fixOrientation) const
    -> CGAL::Polygon_with_holes_2<Kernel>
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

  CGAL::Polygon_2<Kernel> const outer =
      exteriorRing().toPolygon_2(fixOrientation);
  return CGAL::Polygon_with_holes_2<Kernel>(outer, holes.begin(), holes.end());
}

} // namespace SFCGAL
