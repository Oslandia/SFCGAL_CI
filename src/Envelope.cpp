// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/Envelope.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"

#include <memory>

namespace SFCGAL {

Envelope::Envelope()
{
  for (auto &_bound : _bounds) {
    _bound = detail::Interval();
  }
}

Envelope::Envelope(const double &xmin, const double &xmax, const double &ymin,
                   const double &ymax)
{
  _bounds[0] = detail::Interval(xmin, xmax);
  _bounds[1] = detail::Interval(ymin, ymax);
  _bounds[2] = detail::Interval();
}

Envelope::Envelope(const double &xmin, const double &xmax, const double &ymin,
                   const double &ymax, const double &zmin, const double &zmax)
{
  _bounds[0] = detail::Interval(xmin, xmax);
  _bounds[1] = detail::Interval(ymin, ymax);
  _bounds[2] = detail::Interval(zmin, zmax);
}

Envelope::Envelope(const Coordinate &coordinate)
{
  expandToInclude(coordinate);
}

Envelope::Envelope(const Coordinate &coordinate1, const Coordinate &coordinate2)
{
  expandToInclude(coordinate1);
  expandToInclude(coordinate2);
}

Envelope::Envelope(const Envelope &other)
{
  for (size_t i = 0; i < 3; i++) {
    _bounds[i] = other._bounds[i];
  }
}

auto
Envelope::operator=(const Envelope &other) -> Envelope &
{
  for (size_t i = 0; i < 3; i++) {
    _bounds[i] = other._bounds[i];
  }

  return *this;
}

auto
Envelope::isEmpty() const -> bool
{
  return _bounds[0].isEmpty() || _bounds[1].isEmpty();
}

auto
Envelope::is3D() const -> bool
{
  return !isEmpty() && !_bounds[2].isEmpty();
}

void
Envelope::expandToInclude(const Coordinate &coordinate)
{
  if (!coordinate.isEmpty()) {
    _bounds[0].expandToInclude(CGAL::to_double(coordinate.x()));
    _bounds[1].expandToInclude(CGAL::to_double(coordinate.y()));
  }

  if (coordinate.is3D()) {
    _bounds[2].expandToInclude(CGAL::to_double(coordinate.z()));
  }
}

auto
Envelope::contains(const Envelope &envelopeA, const Envelope &envelopeB) -> bool
{
  if (envelopeA.is3D()) {
    return envelopeB.xMin() >= envelopeA.xMin() &&
           envelopeB.xMax() <= envelopeA.xMax() &&
           envelopeB.yMin() >= envelopeA.yMin() &&
           envelopeB.yMax() <= envelopeA.yMax() &&
           envelopeB.zMin() >= envelopeA.zMin() &&
           envelopeB.zMax() <= envelopeA.zMax();
  }

  return envelopeB.xMin() >= envelopeA.xMin() &&
         envelopeB.xMax() <= envelopeA.xMax() &&
         envelopeB.yMin() >= envelopeA.yMin() &&
         envelopeB.yMax() <= envelopeA.yMax();
}

auto
Envelope::overlaps(const Envelope &envelopeA, const Envelope &envelopeB) -> bool
{
  if (envelopeA.is3D()) {
    CGAL::Bbox_3 const abox = envelopeA.toBbox_3();
    CGAL::Bbox_3 const bbox = envelopeB.toBbox_3();
    return CGAL::do_overlap(abox, bbox);
  }

  CGAL::Bbox_2 const abox = envelopeA.toBbox_2();
  CGAL::Bbox_2 const bbox = envelopeB.toBbox_2();
  return CGAL::do_overlap(abox, bbox);
}

auto
Envelope::toRing() const -> std::unique_ptr<LineString>
{
  auto ring = std::make_unique<LineString>();

  if (isEmpty()) {
    return ring;
  }

  ring->addPoint(std::make_unique<Point>(xMin(), yMin()));
  ring->addPoint(std::make_unique<Point>(xMax(), yMin()));
  ring->addPoint(std::make_unique<Point>(xMax(), yMax()));
  ring->addPoint(std::make_unique<Point>(xMin(), yMax()));
  ring->addPoint(ring->startPoint());

  return ring;
}

auto
Envelope::toPolygon() const -> std::unique_ptr<Polygon>
{
  return std::make_unique<Polygon>(toRing().release());
}

auto
Envelope::toShell() const -> std::unique_ptr<PolyhedralSurface>
{
  auto shell = std::make_unique<PolyhedralSurface>();

  if (!is3D()) {
    return shell;
  }

  Point const a(xMin(), yMin(), zMin());
  Point const b(xMax(), yMin(), zMin());
  Point const c(xMax(), yMax(), zMin());
  Point const d(xMin(), yMax(), zMin());

  Point const e(xMin(), yMin(), zMax());
  Point const f(xMax(), yMin(), zMax());
  Point const g(xMax(), yMax(), zMax());
  Point const h(xMin(), yMax(), zMax());

  // bottom : a,d,c,b
  {
    LineString ring;
    ring.addPoint(a);
    ring.addPoint(d);
    ring.addPoint(c);
    ring.addPoint(b);
    ring.addPoint(ring.startPoint());

    shell->addPatch(Polygon(ring));
  }

  // top : e,f,g,h
  {
    LineString ring;
    ring.addPoint(e);
    ring.addPoint(f);
    ring.addPoint(g);
    ring.addPoint(h);
    ring.addPoint(ring.startPoint());

    shell->addPatch(Polygon(ring));
  }

  // front : a,b,f,e
  {
    LineString ring;
    ring.addPoint(a);
    ring.addPoint(b);
    ring.addPoint(f);
    ring.addPoint(e);
    ring.addPoint(ring.startPoint());

    shell->addPatch(Polygon(ring));
  }

  // back : c,d,h,g
  {
    LineString ring;
    ring.addPoint(c);
    ring.addPoint(d);
    ring.addPoint(h);
    ring.addPoint(g);
    ring.addPoint(ring.startPoint());

    shell->addPatch(Polygon(ring));
  }

  // right : b,c,g,f
  {
    LineString ring;
    ring.addPoint(b);
    ring.addPoint(c);
    ring.addPoint(g);
    ring.addPoint(f);
    ring.addPoint(ring.startPoint());

    shell->addPatch(Polygon(ring));
  }

  // left : a,e,h,d
  {
    LineString ring;
    ring.addPoint(a);
    ring.addPoint(e);
    ring.addPoint(h);
    ring.addPoint(d);
    ring.addPoint(ring.startPoint());

    shell->addPatch(Polygon(ring));
  }

  return shell;
}

auto
Envelope::toSolid() const -> std::unique_ptr<Solid>
{
  return std::make_unique<Solid>(toShell().release());
}

auto
Envelope::print(std::ostream &ostr) const -> std::ostream &
{
  ostr << "[ " << xMin();
  ostr << ", " << xMax();
  ostr << ", " << yMin();
  ostr << ", " << yMax();

  if (is3D()) {
    ostr << ", " << zMin() << ", " << zMax();
  }

  ostr << " ]";
  return ostr;
}

/// @private
auto
operator==(const Envelope &lhs, const Envelope &rhs) -> bool
{
  if (lhs.is3D()) {
    return lhs.xMin() == rhs.xMin() && lhs.yMin() == rhs.yMin() &&
           lhs.zMin() == rhs.zMin() && lhs.xMax() == rhs.xMax() &&
           lhs.yMax() == rhs.yMax() && lhs.zMax() == rhs.zMax();
  }

  return lhs.xMin() == rhs.xMin() && lhs.yMin() == rhs.yMin() &&
         lhs.xMax() == rhs.xMax() && lhs.yMax() == rhs.yMax();
}
} // namespace SFCGAL
