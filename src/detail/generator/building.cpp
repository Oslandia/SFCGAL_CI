// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/detail/generator/building.h>

#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/MultiSolid.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>

#include <SFCGAL/algorithm/force3D.h>
#include <SFCGAL/algorithm/orientation.h>

#include <CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>

#include <boost/format.hpp>

namespace SFCGAL::generator {

using Point_2              = Kernel::Point_2;
using Point_3              = Kernel::Point_3;
using Polygon_2            = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;
using Straight_skeleton_2  = CGAL::Straight_skeleton_2<Kernel>;

/**
 * @brief Basic building generator relying on a straight skeleton
 */
auto
building(const Polygon &g, const Kernel::FT &wallHeight,
         const Kernel::FT &roofSlope) -> std::unique_ptr<Geometry>;
/**
 * @brief Basic building generator relying on a straight skeleton
 */
auto
building(const MultiPolygon &g, const Kernel::FT &wallHeight,
         const Kernel::FT &roofSlope) -> std::unique_ptr<Geometry>;

void
_buildingWall(const Polygon_2 &ring, const Kernel::FT &wallHeight,
              PolyhedralSurface &shell)
{
  size_t const npt = ring.size();

  for (size_t i = 0; i < npt; i++) {
    const Point_2 &a = ring.vertex(i);
    const Point_2 &b = ring.vertex((i + 1) % npt);

    LineString wallRing;
    wallRing.addPoint(new Point(a.x(), a.y(), Kernel::FT(0)));
    wallRing.addPoint(new Point(b.x(), b.y(), Kernel::FT(0)));
    wallRing.addPoint(new Point(b.x(), b.y(), wallHeight));
    wallRing.addPoint(new Point(a.x(), a.y(), wallHeight));
    wallRing.addPoint(new Point(a.x(), a.y(), Kernel::FT(0)));
    shell.addPolygon(Polygon(wallRing));
  }
}

///
///
///
auto
building(const Polygon &g, const Kernel::FT &wallHeight,
         const Kernel::FT &roofSlope) -> std::unique_ptr<Geometry>
{
  // typedef Straight_skeleton_2::Vertex_const_handle     Vertex_const_handle ;
  using Halfedge_const_handle = Straight_skeleton_2::Halfedge_const_handle;
  // typedef Straight_skeleton_2::Halfedge_const_iterator
  // Halfedge_const_iterator ;
  using Face_const_iterator = Straight_skeleton_2::Face_const_iterator;

  // convert to CGAL polygon and generate straight skeleton
  Polygon_with_holes_2 polygon = g.toPolygon_with_holes_2();

  // fix orientation
  algorithm::makeValidOrientation(polygon);

  boost::shared_ptr<Straight_skeleton_2> const skeleton =
      create_interior_straight_skeleton_2(
          polygon.outer_boundary().vertices_begin(),
          polygon.outer_boundary().vertices_end(), polygon.holes_begin(),
          polygon.holes_end(), Kernel());

  std::unique_ptr<PolyhedralSurface> shell(new PolyhedralSurface);
  // bottom part
  {
    Polygon bottom(polygon);
    bottom.reverse();
    algorithm::force3D(bottom);
    shell->addPolygon(bottom);
  }

  // walls
  {
    // exterior rings
    _buildingWall(polygon.outer_boundary(), wallHeight, *shell);

    // interior rings
    for (auto it = polygon.holes_begin();
         it != polygon.holes_end(); ++it) {
      _buildingWall(*it, wallHeight, *shell);
    }
  }

  // roof
  {
    for (Face_const_iterator it = skeleton->faces_begin();
         it != skeleton->faces_end(); ++it) {

      LineString            roofFaceRing;
      Halfedge_const_handle h                 = it->halfedge();
      Halfedge_const_handle const done(h);
      bool                  infiniteTimeFound = false;

      do {
        infiniteTimeFound = infiniteTimeFound || h->has_infinite_time();

        Point_2    const point  = h->vertex()->point();
        Kernel::FT const zPoint = wallHeight + h->vertex()->time() * roofSlope;

        roofFaceRing.addPoint(Point(point.x(), point.y(), zPoint));

        h = h->next();
      } while (h != done && !infiniteTimeFound);

      if (!infiniteTimeFound) {
        roofFaceRing.addPoint(roofFaceRing.startPoint());
        shell->addPolygon(Polygon(roofFaceRing));
      }
    }
  }

  return std::unique_ptr<Geometry>(new Solid(shell.release()));
}

///
///
///
auto
building(const MultiPolygon &g, const Kernel::FT &wallHeight,
         const Kernel::FT &roofSlope) -> std::unique_ptr<Geometry>
{
  std::unique_ptr<MultiSolid> multiSolid(new MultiSolid);

  for (size_t i = 0; i < g.numGeometries(); i++) {
    multiSolid->addGeometry(
        building(g.polygonN(i), wallHeight, roofSlope).release());
  }

  return std::unique_ptr<Geometry>(multiSolid.release());
}

///
///
///
auto
building(const Geometry &g, const Kernel::FT &wallHeight,
         const Kernel::FT &roofSlope) -> std::unique_ptr<Geometry>
{
  switch (g.geometryTypeId()) {
  case TYPE_POLYGON:
    return building(g.as<Polygon>(), wallHeight, roofSlope);

  case TYPE_MULTIPOLYGON:
    return building(g.as<MultiPolygon>(), wallHeight, roofSlope);

  default:
    BOOST_THROW_EXCEPTION(Exception(
        (boost::format("bad geometry type (%s) in generator::building") %
         g.geometryType())
            .str()));
    return std::unique_ptr<Geometry>(new GeometryCollection());
  }
}

} // namespace SFCGAL
