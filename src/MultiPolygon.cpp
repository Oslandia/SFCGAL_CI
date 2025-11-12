// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/GeometryVisitor.h"
#include "SFCGAL/version.h"

namespace SFCGAL {

MultiPolygon::MultiPolygon() = default;

MultiPolygon::MultiPolygon(MultiPolygon const &other)

    = default;

#if SFCGAL_CGAL_VERSION_MAJOR >= 6
MultiPolygon::MultiPolygon(const CGAL::Multipolygon_with_holes_2<Kernel> &other)
{
  for (const auto &pwh : other.polygons_with_holes()) {
    auto polygon = std::make_unique<Polygon>();

    // Convert exterior ring
    if (!pwh.outer_boundary().is_empty()) {
      LineString exterior;
      for (auto vertex_it = pwh.outer_boundary().vertices_begin();
           vertex_it != pwh.outer_boundary().vertices_end(); ++vertex_it) {
        exterior.addPoint(Point(CGAL::to_double(vertex_it->x()),
                                CGAL::to_double(vertex_it->y())));
      }
      // Close the ring
      if (!exterior.isEmpty()) {
        exterior.addPoint(exterior.startPoint());
        polygon->setExteriorRing(exterior);
      }
    }

    // Convert holes
    for (auto hole_it = pwh.holes_begin(); hole_it != pwh.holes_end();
         ++hole_it) {
      LineString hole_ring;
      for (auto vertex_it = hole_it->vertices_begin();
           vertex_it != hole_it->vertices_end(); ++vertex_it) {
        hole_ring.addPoint(Point(CGAL::to_double(vertex_it->x()),
                                 CGAL::to_double(vertex_it->y())));
      }
      // Close the ring
      if (!hole_ring.isEmpty()) {
        hole_ring.addPoint(hole_ring.startPoint());
        polygon->addInteriorRing(hole_ring);
      }
    }

    addGeometry(polygon.release());
  }
}
#endif

auto
MultiPolygon::operator=(MultiPolygon other) -> MultiPolygon &
{
  swap(other);
  return *this;
}

MultiPolygon::~MultiPolygon() = default;

auto
MultiPolygon::geometryType() const -> std::string
{
  return "MultiPolygon";
}

auto
MultiPolygon::geometryTypeId() const -> GeometryType
{
  return TYPE_MULTIPOLYGON;
}

auto
MultiPolygon::isAllowed(Geometry const &geometry) -> bool
{
  return geometry.geometryTypeId() == TYPE_POLYGON;
}

void
MultiPolygon::accept(GeometryVisitor &visitor)
{
  return visitor.visit(*this);
}

void
MultiPolygon::accept(ConstGeometryVisitor &visitor) const
{
  return visitor.visit(*this);
}

#if SFCGAL_CGAL_VERSION_MAJOR >= 6
auto
MultiPolygon::toMultipolygon_with_holes_2(bool fixOrientation) const
    -> CGAL::Multipolygon_with_holes_2<Kernel>
{
  CGAL::Multipolygon_with_holes_2<Kernel> mp;

  for (size_t i = 0; i < numGeometries(); ++i) {
    const Polygon &polygon = polygonN(i);
    if (!polygon.isEmpty()) {
      mp.add_polygon_with_holes(polygon.toPolygon_with_holes_2(fixOrientation));
    }
  }

  return mp;
}
#endif

} // namespace SFCGAL
