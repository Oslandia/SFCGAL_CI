// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/GeometrySet.h"

#include "SFCGAL/Curve.h"
#include "SFCGAL/Envelope.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/detail/TypeForDimension.h"

#include "SFCGAL/algorithm/connection.h"
#include "SFCGAL/algorithm/covers.h"
#include "SFCGAL/algorithm/volume.h"

#include <CGAL/Bbox_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <boost/graph/adjacency_list.hpp>

#include <limits>
#include <map>

auto
operator<(const CGAL::Segment_2<SFCGAL::Kernel> &sega,
          const CGAL::Segment_2<SFCGAL::Kernel> &segb) -> bool
{
  if (sega.source() == segb.source()) {
    return sega.target() < segb.target();
  }

  return sega.source() < segb.source();
}

auto
operator<(const CGAL::Segment_3<SFCGAL::Kernel> &sega,
          const CGAL::Segment_3<SFCGAL::Kernel> &segb) -> bool
{
  if (sega.source() == segb.source()) {
    return sega.target() < segb.target();
  }

  return sega.source() < segb.source();
}

namespace SFCGAL::detail {

/// \cond IGNORE
void
_decompose_triangle(const Triangle                    &tri,
                    GeometrySet<2>::SurfaceCollection &surfaces,
                    dim_t<2> /*unused*/)
{
  CGAL::Polygon_2<Kernel> outer;
  outer.push_back(tri.vertex(0).toPoint_2());
  outer.push_back(tri.vertex(1).toPoint_2());
  outer.push_back(tri.vertex(2).toPoint_2());

  if (outer.orientation() == CGAL::CLOCKWISE) {
    outer.reverse_orientation();
  }

  surfaces.emplace_back(CGAL::Polygon_with_holes_2<Kernel>(outer));
}
void
_decompose_triangle(const Triangle                    &tri,
                    GeometrySet<3>::SurfaceCollection &surfaces,
                    dim_t<3> /*unused*/)
{
  CGAL::Triangle_3<Kernel> const outtri(tri.vertex(0).toPoint_3(),
                                        tri.vertex(1).toPoint_3(),
                                        tri.vertex(2).toPoint_3());
  surfaces.emplace_back(outtri);
}

void
_decompose_polygon(const Polygon                     &poly,
                   GeometrySet<2>::SurfaceCollection &surfaces,
                   dim_t<2> /*unused*/)
{
  BOOST_ASSERT(!poly.isEmpty());
  surfaces.emplace_back(poly.toPolygon_with_holes_2());
}
void
_decompose_polygon(const Polygon                     &poly,
                   GeometrySet<3>::SurfaceCollection &surfaces,
                   dim_t<3> /*unused*/)
{
  BOOST_ASSERT(!poly.isEmpty());
  TriangulatedSurface surf;
  triangulate::triangulatePolygon3D(poly, surf);

  for (size_t i = 0; i < surf.numPatches(); ++i) {
    const Triangle &tri = surf.patchN(i);
    surfaces.emplace_back(CGAL::Triangle_3<Kernel>(tri.vertex(0).toPoint_3(),
                                                   tri.vertex(1).toPoint_3(),
                                                   tri.vertex(2).toPoint_3()));
  }
}

void
_decompose_solid(const Solid & /*unused*/,
                 GeometrySet<2>::VolumeCollection & /*unused*/,
                 dim_t<2> /*unused*/)
{
}
void
_decompose_solid(const Solid &solid, GeometrySet<3>::VolumeCollection &volumes,
                 dim_t<3> /*unused*/)
{
  BOOST_ASSERT(!solid.isEmpty());
  // volume orientation test
  // TODO: simplfiy ?
  MarkedPolyhedron p =
      *solid.exteriorShell().toPolyhedron_3<Kernel, MarkedPolyhedron>();

  if (algorithm::volume(solid) < 0) {
    // if the volume is "inverted", we reverse it
    // TODO: Once every boolean operations work with complement geometries, we
    // may want to keep the solid inverted
    p.inside_out();
  }

  volumes.emplace_back(p);
}
/// \endcond

template <int Dim>
GeometrySet<Dim>::GeometrySet() = default;

template <int Dim>
GeometrySet<Dim>::GeometrySet(const Geometry &g)
{
  _decompose(g);
}

template <int Dim>
GeometrySet<Dim>::GeometrySet(const typename TypeForDimension<Dim>::Point &g,
                              int /*flags*/)
{
  addPrimitive(g);
}

template <int Dim>
GeometrySet<Dim>::GeometrySet(const typename TypeForDimension<Dim>::Segment &g,
                              int /*flags*/)
{
  addPrimitive(g);
}

template <int Dim>
GeometrySet<Dim>::GeometrySet(const typename TypeForDimension<Dim>::Surface &g,
                              int /*flags*/)
{
  addPrimitive(g);
}

template <int Dim>
GeometrySet<Dim>::GeometrySet(const typename TypeForDimension<Dim>::Volume &g,
                              int /*flags*/)
{
  addPrimitive(g);
}

template <int Dim>
void
GeometrySet<Dim>::merge(const GeometrySet<Dim> &g)
{
  std::copy(g.points().begin(), g.points().end(),
            std::inserter(points(), points().end()));
  std::copy(g.segments().begin(), g.segments().end(),
            std::inserter(segments(), segments().end()));
  std::copy(g.surfaces().begin(), g.surfaces().end(),
            std::back_inserter(surfaces()));
  std::copy(g.volumes().begin(), g.volumes().end(),
            std::back_inserter(volumes()));
}

template <int Dim>
void
GeometrySet<Dim>::addGeometry(const Geometry &g)
{
  _decompose(g);
}

/**
 * Add primitive to 2D geometry set from primitive handle
 * @param p The primitive handle to add
 */
template <>
void
GeometrySet<2>::addPrimitive(const PrimitiveHandle<2> &p)
{
  switch (p.handle.which()) {
  case PrimitivePoint:
    _points.insert(*boost::get<const TypeForDimension<2>::Point *>(p.handle));
    break;

  case PrimitiveSegment:
    _segments.insert(
        *boost::get<const TypeForDimension<2>::Segment *>(p.handle));
    break;

  case PrimitiveSurface:
    _surfaces.emplace_back(
        *boost::get<const TypeForDimension<2>::Surface *>(p.handle));
    break;

  default:
    break;
  }
}

/**
 * Add primitive to 3D geometry set from primitive handle
 * @param p The primitive handle to add
 */
template <>
void
GeometrySet<3>::addPrimitive(const PrimitiveHandle<3> &p)
{
  switch (p.handle.which()) {
  case PrimitivePoint:
    _points.insert(*boost::get<const TypeForDimension<3>::Point *>(p.handle));
    break;

  case PrimitiveSegment:
    _segments.insert(
        *boost::get<const TypeForDimension<3>::Segment *>(p.handle));
    break;

  case PrimitiveSurface:
    _surfaces.emplace_back(
        *boost::get<const TypeForDimension<3>::Surface *>(p.handle));
    break;

  case PrimitiveVolume: {
    const TypeForDimension<3>::Volume &vol =
        *boost::get<const TypeForDimension<3>::Volume *>(p.handle);
    BOOST_ASSERT(!vol.empty());
    _volumes.emplace_back(vol);
    break;
  }
  }
}

/**
 * Add primitive to 3D geometry set from CGAL object
 * @param o The CGAL object to add as primitive
 * @param pointsAsRing If true, build polygon from point vector
 */
template <>
void
GeometrySet<3>::addPrimitive(const CGAL::Object &o, bool pointsAsRing)
{
  using TPoint   = TypeForDimension<3>::Point;
  using TSegment = TypeForDimension<3>::Segment;
  using TSurface = TypeForDimension<3>::Surface;
  using TVolume  = TypeForDimension<3>::Volume;

  if (const auto *p = CGAL::object_cast<TPoint>(&o)) {
    _points.insert(TPoint(*p));
  } else if (const auto *pts = CGAL::object_cast<std::vector<TPoint>>(&o)) {
    if (pointsAsRing) {
      // if pointsAsRing is true, build a polygon out of points
      // FIXME : we use triangulation here, which is not needed
      // We should have created a (planar) Polyhedron directly out of the points
      LineString ls;

      for (const auto &pt : *pts) {
        ls.addPoint(pt);
      }

      // close the ring
      ls.addPoint((*pts)[0]);
      Polygon const poly(ls);
      _decompose_polygon(poly, _surfaces, dim_t<3>());
    } else {
      std::copy(pts->begin(), pts->end(),
                std::inserter(_points, _points.end()));
    }
  } else if (const auto *p = CGAL::object_cast<TSegment>(&o)) {
    _segments.insert(TSegment(*p));
  } else if (const auto *p = CGAL::object_cast<TSurface>(&o)) {
    _surfaces.emplace_back(TSurface(*p));
  } else if (const auto *p = CGAL::object_cast<TVolume>(&o)) {
    BOOST_ASSERT(!p->empty());
    _volumes.emplace_back(TVolume(*p));
  }
}

/**
 * Add primitive to 2D geometry set from CGAL object
 * @param o The CGAL object to add as primitive
 * @param pointsAsRing If true, build polygon from point vector
 */
template <>
void
GeometrySet<2>::addPrimitive(const CGAL::Object &o, bool pointsAsRing)
{
  using TPoint   = TypeForDimension<2>::Point;
  using TSegment = TypeForDimension<2>::Segment;
  using TSurface = TypeForDimension<2>::Surface;
  using TVolume  = TypeForDimension<2>::Volume;

  if (const auto *p = CGAL::object_cast<TPoint>(&o)) {
    _points.insert(TPoint(*p));
  } else if (const auto *pts = CGAL::object_cast<std::vector<TPoint>>(&o)) {
    if (pointsAsRing) {
      // if pointsAsRing is true, build a polygon out of points
      CGAL::Polygon_2<Kernel> poly;

      for (const auto &pt : *pts) {
        poly.push_back(pt);
      }

      CGAL::Polygon_with_holes_2<Kernel> const polyh(poly);
      _surfaces.emplace_back(polyh);
    } else {
      std::copy(pts->begin(), pts->end(),
                std::inserter(_points, _points.end()));
    }
  } else if (const auto *tri =
                 CGAL::object_cast<CGAL::Triangle_2<Kernel>>(&o)) {
    // convert to a polygon
    CGAL::Polygon_2<Kernel> poly;
    poly.push_back(tri->vertex(0));
    poly.push_back(tri->vertex(1));
    poly.push_back(tri->vertex(2));
    CGAL::Polygon_with_holes_2<Kernel> const polyh(poly);
    _surfaces.emplace_back(polyh);
  } else if (const auto *p = CGAL::object_cast<TSegment>(&o)) {
    _segments.insert(TSegment(*p));
  } else if (const auto *p = CGAL::object_cast<TSurface>(&o)) {
    BOOST_ASSERT(!p->is_unbounded());
    _surfaces.emplace_back(TSurface(*p));
  } else if (const auto *p = CGAL::object_cast<TVolume>(&o)) {
    _volumes.emplace_back(TVolume(*p));
  }
}

template <int Dim>
void
GeometrySet<Dim>::addPrimitive(const typename TypeForDimension<Dim>::Point &p,
                               int flags)
{
  _points.insert(CollectionElement<typename Point_d<Dim>::Type>(p, flags));
}

template <int Dim>
void
GeometrySet<Dim>::addPrimitive(const typename TypeForDimension<Dim>::Segment &p,
                               int flags)
{
  _segments.insert(CollectionElement<typename Segment_d<Dim>::Type>(p, flags));
}

/**
 * Add 2D surface primitive to geometry set
 * @param p The surface to add
 * @param flags Optional flags for the surface
 */
template <>
void
GeometrySet<2>::addPrimitive(const TypeForDimension<2>::Surface &p, int flags)
{
  BOOST_ASSERT(!p.is_unbounded());
  _surfaces.emplace_back(p);
  _surfaces.back().setFlags(flags);
}
/**
 * Add 3D surface primitive to geometry set
 * @param p The surface to add
 * @param flags Optional flags for the surface
 */
template <>
void
GeometrySet<3>::addPrimitive(const TypeForDimension<3>::Surface &p, int flags)
{
  _surfaces.emplace_back(p);
  _surfaces.back().setFlags(flags);
}

/**
 * Add 2D volume primitive to geometry set (no-op for 2D)
 * @param volume The volume (unused in 2D)
 * @param flags The flags (unused)
 */
template <>
void
GeometrySet<2>::addPrimitive(
    [[maybe_unused]] const TypeForDimension<2>::Volume &volume,
    [[maybe_unused]] int                                flags)
{
}

/**
 * Add 3D volume primitive to geometry set
 * @param p The volume to add
 * @param flags Optional flags for the volume
 */
template <>
void
GeometrySet<3>::addPrimitive(const TypeForDimension<3>::Volume &p, int flags)
{
  BOOST_ASSERT(!p.empty());

  if (p.is_closed()) {
    _volumes.emplace_back(p, flags);
  } else {
    // it is an unclosed volume, i.e. a surface
    BOOST_ASSERT(p.is_pure_triangle());
    CGAL::Point_3<Kernel> p1;
    CGAL::Point_3<Kernel> p2;
    CGAL::Point_3<Kernel> p3;

    for (MarkedPolyhedron::Facet_const_iterator fit = p.facets_begin();
         fit != p.facets_end(); ++fit) {
      MarkedPolyhedron::Halfedge_around_facet_const_circulator cit =
          fit->facet_begin();
      p1 = cit->vertex()->point();
      cit++;
      p2 = cit->vertex()->point();
      cit++;
      p3 = cit->vertex()->point();
      CGAL::Triangle_3<Kernel> const tri(p1, p2, p3);
      _surfaces.emplace_back(tri);
    }
  }
}

template <int Dim>
auto
GeometrySet<Dim>::hasPoints() const -> bool
{
  return !points().empty();
}

template <int Dim>
auto
GeometrySet<Dim>::hasSegments() const -> bool
{
  return !segments().empty();
}

/**
 * Check if 2D geometry set has surfaces
 * @return True if the set contains surfaces
 */
template <>
auto
GeometrySet<2>::hasSurfaces() const -> bool
{
  return !surfaces().empty();
}
/**
 * Check if 3D geometry set has surfaces
 * @return True if the set contains surfaces or volumes
 */
template <>
auto
GeometrySet<3>::hasSurfaces() const -> bool
{
  if (!surfaces().empty()) {
    return true;
  }

  if (!volumes().empty()) {
    for (const auto &_volume : _volumes) {
      if (!_volume.primitive().is_closed()) {
        return true;
      }
    }
  }

  return false;
}

/**
 * Check if 2D geometry set has volumes
 * @return Always false for 2D sets
 */
template <>
auto
GeometrySet<2>::hasVolumes() const -> bool
{
  return false;
}
/**
 * Check if 3D geometry set has volumes
 * @return True if the set contains volumes
 */
template <>
auto
GeometrySet<3>::hasVolumes() const -> bool
{
  if (!volumes().empty()) {
    return true;
  }

  if (!volumes().empty()) {
    for (const auto &_volume : _volumes) {
      if (_volume.primitive().is_closed()) {
        return true;
      }
    }
  }

  return false;
}
template <int Dim>
void
GeometrySet<Dim>::_decompose(const Geometry &g)
{
  if (g.isEmpty()) {
    return;
  }

  if (g.is<GeometryCollection>()) {
    const auto &collect = g.as<GeometryCollection>();

    for (size_t i = 0; i < g.numGeometries(); ++i) {
      _decompose(collect.geometryN(i));
    }

    return;
  }

  switch (g.geometryTypeId()) {
  case TYPE_POINT:
    _points.insert(g.as<Point>().toPoint_d<Dim>());
    break;

  case TYPE_LINESTRING: {
    const auto &ls = g.as<LineString>();

    for (size_t i = 0; i < ls.numPoints() - 1; ++i) {
      typename TypeForDimension<Dim>::Segment const seg(
          ls.pointN(i).toPoint_d<Dim>(), ls.pointN(i + 1).toPoint_d<Dim>());
      _segments.insert(seg);
    }

    break;
  }

  case TYPE_NURBSCURVE: {
    // Convert NURBS curve to LineString approximation
    auto lineString = g.as<NURBSCurve>().toLineString(); // default parameters
    if (!lineString || lineString->isEmpty()) {
      return; // Cannot decompose empty curve
    }
    for (size_t i = 0; i < lineString->numPoints() - 1; ++i) {
      typename TypeForDimension<Dim>::Segment const seg(
          lineString->pointN(i).toPoint_d<Dim>(),
          lineString->pointN(i + 1).toPoint_d<Dim>());
      _segments.insert(seg);
    }
    break;
  }

  case TYPE_TRIANGLE: {
    _decompose_triangle(g.as<Triangle>(), _surfaces, dim_t<Dim>());
    break;
  }

  case TYPE_POLYGON: {
    _decompose_polygon(g.as<Polygon>(), _surfaces, dim_t<Dim>());
    break;
  }

  case TYPE_TRIANGULATEDSURFACE: {
    const auto &tri = g.as<TriangulatedSurface>();

    for (size_t i = 0; i < tri.numPatches(); ++i) {
      _decompose(tri.patchN(i));
    }

    break;
  }

  case TYPE_POLYHEDRALSURFACE: {
    const auto &tri = g.as<PolyhedralSurface>();

    for (size_t i = 0; i < tri.numPatches(); ++i) {
      _decompose(tri.patchN(i));
    }

    break;
  }

  case TYPE_SOLID: {
    const auto &solid = g.as<Solid>();
    _decompose_solid(solid, _volumes, dim_t<Dim>());
    break;
  }

  default:
    break;
  }
}

template <int Dim>
void
GeometrySet<Dim>::computeBoundingBoxes(
    typename HandleCollection<Dim>::Type &handles,
    typename BoxCollection<Dim>::Type    &boxes) const
{
  boxes.clear();

  for (auto it = _points.begin(); it != _points.end(); ++it) {
    const typename TypeForDimension<Dim>::Point *pt = &(it->primitive());
    PrimitiveHandle<Dim> const                   h(pt);
    handles.push_back(h);
    boxes.push_back(typename PrimitiveBox<Dim>::Type(it->primitive().bbox(),
                                                     &handles.back()));
  }

  for (auto it = _segments.begin(); it != _segments.end(); ++it) {
    handles.push_back(PrimitiveHandle<Dim>(&(it->primitive())));
    boxes.push_back(typename PrimitiveBox<Dim>::Type(it->primitive().bbox(),
                                                     &handles.back()));
  }

  for (auto it = _surfaces.begin(); it != _surfaces.end(); ++it) {
    handles.push_back(PrimitiveHandle<Dim>(&(it->primitive())));
    boxes.push_back(typename PrimitiveBox<Dim>::Type(it->primitive().bbox(),
                                                     &handles.back()));
  }

  for (auto it = _volumes.begin(); it != _volumes.end(); ++it) {
    handles.push_back(PrimitiveHandle<Dim>(&(it->primitive())));
    boxes.push_back(typename PrimitiveBox<Dim>::Type(
        compute_solid_bbox(it->primitive(), dim_t<Dim>()), &handles.back()));
  }
}

/// \cond IGNORE
template <int Dim>
void
recompose_points(const typename GeometrySet<Dim>::PointCollection &points,
                 std::vector<Geometry *> &rpoints, dim_t<Dim> /*unused*/)
{
  if (points.empty()) {
    return;
    //			rpoints.push_back( new Point() );
  }
  for (auto it = points.begin(); it != points.end(); ++it) {
    rpoints.push_back(new Point(it->primitive()));
  }
}
/// \endcond

/**
 * @brief Comparator for sorting points lexicographically
 */
struct ComparePoints {
  /**
   * @brief Compare 2D points lexicographically
   * @param lhs Left-hand side point
   * @param rhs Right-hand side point
   * @return True if lhs is lexicographically smaller than rhs
   */
  auto
  operator()(const CGAL::Point_2<Kernel> &lhs,
             const CGAL::Point_2<Kernel> &rhs) const -> bool
  {
    return lhs.x() == rhs.x() ? lhs.y() < rhs.y() : lhs.x() < rhs.x();
  }
  /**
   * @brief Compare 3D points lexicographically
   * @param lhs Left-hand side point
   * @param rhs Right-hand side point
   * @return True if lhs is lexicographically smaller than rhs
   */
  auto
  operator()(const CGAL::Point_3<Kernel> &lhs,
             const CGAL::Point_3<Kernel> &rhs) const -> bool
  {
    return lhs.x() == rhs.x()
               ? (lhs.y() == rhs.y() ? lhs.z() < rhs.z() : lhs.y() < rhs.y())
               : lhs.x() < rhs.x();
  }
};

/// \cond IGNORE
template <int Dim>
void
recompose_segments(const typename GeometrySet<Dim>::SegmentCollection &segments,
                   std::vector<Geometry *> &lines, dim_t<Dim> /*unused*/)
{
  if (segments.empty()) {
    //			lines.push_back( new LineString );
    return;
  }

  // what we need is a graph, we do a depth first traversal and stop a
  // linestring when more than one segment is connected first we need to label
  // vertices then build the graph and traverse depth first
  std::vector<Point> points;
  using Edge = std::pair<int, int>;
  std::vector<Edge> edges;
  {
    using PointMap = typename std::map<typename TypeForDimension<Dim>::Point,
                                       int, ComparePoints>;
    PointMap pointMap;

    for (auto it = segments.begin(); it != segments.end(); ++it) {
      const auto foundSource = pointMap.find(it->primitive().source());
      const auto foundTarget = pointMap.find(it->primitive().target());
      const int  sourceId =
          foundSource != pointMap.end() ? foundSource->second : points.size();

      if (foundSource == pointMap.end()) {
        points.push_back(it->primitive().source());
        pointMap[it->primitive().source()] = sourceId;
      }

      const int targetId =
          foundTarget != pointMap.end() ? foundTarget->second : points.size();

      if (foundTarget == pointMap.end()) {
        points.push_back(it->primitive().target());
        pointMap[it->primitive().target()] = targetId;
      }

      edges.emplace_back(sourceId, targetId);
    }
  }

  using Graph = boost::adjacency_list<
      boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property,
      boost::property<boost::edge_color_t, boost::default_color_type>>;
  Graph g(edges.begin(), edges.end(), edges.size());

  // now we find all branches without bifurcations,

  boost::graph_traits<Graph>::edge_iterator ei;
  boost::graph_traits<Graph>::edge_iterator ei_end;

  for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
    if (boost::get(boost::edge_color, g)[*ei] == boost::white_color) {
      // not already marked, find the first ancestor with multiple connections,
      // or no connections or self (in case of a loop)
      boost::graph_traits<Graph>::edge_descriptor root = *ei;
      {
        boost::graph_traits<Graph>::in_edge_iterator ej;
        boost::graph_traits<Graph>::in_edge_iterator ek;

        for (boost::tie(ej, ek) = boost::in_edges(boost::source(root, g), g);
             ek - ej == 1 && *ej != *ei;
             boost::tie(ej, ek) = boost::in_edges(boost::source(root, g), g)) {
          root = *ej;
        }
      }

      // now we go down
      auto *line = new LineString;
      lines.push_back(line);
      line->addPoint(points[boost::source(root, g)]);
      line->addPoint(points[boost::target(root, g)]);
      boost::get(boost::edge_color, g)[root] = boost::black_color;

      boost::graph_traits<Graph>::out_edge_iterator ej;
      boost::graph_traits<Graph>::out_edge_iterator ek;

      for (boost::tie(ej, ek) = boost::out_edges(boost::target(root, g), g);
           ek - ej == 1 && *ej != root;
           boost::tie(ej, ek) = boost::out_edges(boost::target(*ej, g), g)) {
        line->addPoint(points[boost::target(*ej, g)]);
        boost::get(boost::edge_color, g)[*ej] = boost::black_color;
      }
    }
  }
}

void
recompose_surfaces(const GeometrySet<2>::SurfaceCollection &surfaces,
                   std::vector<Geometry *> &output, dim_t<2> /*unused*/)
{
  for (const auto &surface : surfaces) {
    if (surface.primitive().holes_begin() == surface.primitive().holes_end() &&
        surface.primitive().outer_boundary().size() == 3) {
      auto vit = surface.primitive().outer_boundary().vertices_begin();
      CGAL::Point_2<Kernel> const p1(*vit++);
      CGAL::Point_2<Kernel> const p2(*vit++);
      CGAL::Point_2<Kernel> const p3(*vit++);
      output.push_back(new Triangle(CGAL::Triangle_2<Kernel>(p1, p2, p3)));
    } else {
      output.push_back(new Polygon(surface.primitive()));
    }
  }
}

void
recompose_surfaces(const GeometrySet<3>::SurfaceCollection &surfaces,
                   std::vector<Geometry *> &output, dim_t<3> /*unused*/)
{
  if (surfaces.empty()) {
    return;
  }

  // TODO : regroup triangles of the same mesh
  if (surfaces.size() == 1) {
    output.push_back(new Triangle(surfaces.begin()->primitive()));
    return;
  }

  std::unique_ptr<TriangulatedSurface> tri(new TriangulatedSurface);

  for (const auto &surface : surfaces) {
    tri->addPatch(new Triangle(surface.primitive()));
  }

  algorithm::SurfaceGraph const graph(*tri);
  std::vector<size_t> component(boost::num_vertices(graph.faceGraph()));
  BOOST_ASSERT(tri->numPatches() == component.size());
  const size_t numComponents =
      boost::connected_components(graph.faceGraph(), component.data());

  if (1 == numComponents) {
    output.push_back(tri.release());
  } else {
    std::vector<TriangulatedSurface *> sout(numComponents);

    for (unsigned c = 0; c < numComponents; c++) {
      sout[c] = new TriangulatedSurface;
      output.push_back(sout[c]);
    }

    const size_t numPatches = tri->numPatches();

    for (size_t t = 0; t != numPatches; ++t) {
      sout[component[t]]->addPatch(tri->patchN(t));
    }
  }
}

void
recompose_volumes(const GeometrySet<2>::VolumeCollection & /*unused*/,
                  std::vector<Geometry *> & /*unused*/, dim_t<2> /*unused*/)
{
}

void
recompose_volumes(const GeometrySet<3>::VolumeCollection &volumes,
                  std::vector<Geometry *> &output, dim_t<3> /*unused*/)
{
  if (volumes.empty()) {
    return;
  }

  for (const auto &volume : volumes) {
    if ((volume.flags() & FLAG_IS_PLANAR) != 0) {
      // extract the boundary
      std::list<CGAL::Point_3<Kernel>> boundary;

      for (MarkedPolyhedron::Halfedge_const_iterator it =
               volume.primitive().halfedges_begin();
           it != volume.primitive().halfedges_end(); ++it) {
        if (!it->is_border()) {
          continue;
        }

        CGAL::Point_3<Kernel> const p1 = it->prev()->vertex()->point();
        CGAL::Point_3<Kernel> const p2 = it->vertex()->point();

        // TODO: test for colinearity
        // Temporary vertice may have been introduced during triangulations
        // and since we expect here a planar surface, it is safe to simplify the
        // boundary by eliminating collinear points.
        if (boundary.empty()) {
          boundary.push_back(p1);
          boundary.push_back(p2);
        } else if (boundary.back() == p1) {
          boundary.push_back(p2);
        } else if (boundary.front() == p2) {
          boundary.push_front(p1);
        }
      }

      if (boundary.size() == 3) {
        // It is a triangle

        Point p[3];
        auto  it = boundary.begin();

        for (size_t i = 0; i < 3; ++i, ++it) {
          p[i] = *it;
        }

        output.push_back(new Triangle(p[0], p[1], p[2]));
      } else {
        // Else it is a polygon
        auto *ls = new LineString;

        for (auto &it : boundary) {
          ls->addPoint(it);
        }

        output.push_back(new Polygon(ls));
      }
    } else {

      auto *shell = new PolyhedralSurface(volume.primitive());
      // TODO: test open / closed
      output.push_back(new Solid(shell));
    }
  }
}
/// \endcond

template <int Dim>
auto
GeometrySet<Dim>::recompose() const -> std::unique_ptr<Geometry>
{
  std::vector<Geometry *> geometries;

  recompose_points(_points, geometries, dim_t<Dim>());
  recompose_segments(_segments, geometries, dim_t<Dim>());
  recompose_surfaces(_surfaces, geometries, dim_t<Dim>());
  recompose_volumes(_volumes, geometries, dim_t<Dim>());

  if (geometries.empty()) {
    return std::unique_ptr<Geometry>(new GeometryCollection);
  }

  if (geometries.size() == 1) {
    return std::unique_ptr<Geometry>(geometries[0]);
  }

  // else we have a mix of different types
  bool      hasCommonType = true;
  int const commonType    = geometries[0]->geometryTypeId();

  for (auto &geometrie : geometries) {
    if (geometrie->geometryTypeId() != commonType) {
      hasCommonType = false;
      break;
    }
  }

  GeometryCollection *ret = nullptr;

  if (hasCommonType) {
    if (commonType == TYPE_POINT) {
      ret = new MultiPoint;
    } else if (commonType == TYPE_LINESTRING) {
      ret = new MultiLineString;
    } else if (commonType == TYPE_POLYGON) {
      ret = new MultiPolygon;
    } else if (commonType == TYPE_SOLID) {
      ret = new MultiSolid;
    } else {
      // one common type, but no MULTI equivalent
      ret = new GeometryCollection;
    }
  } else {
    ret = new GeometryCollection;
  }

  BOOST_ASSERT(ret != 0);

  for (auto &geometrie : geometries) {
    ret->addGeometry(geometrie);
  }

  return std::unique_ptr<Geometry>(ret);
}

/// \cond IGNORE
void
_collect_points(const CGAL::Polygon_with_holes_2<Kernel> &poly,
                GeometrySet<2>::PointCollection          &points)
{
  for (auto vit = poly.outer_boundary().vertices_begin();
       vit != poly.outer_boundary().vertices_end(); ++vit) {
    points.insert(*vit);
  }

  for (auto hit = poly.holes_begin(); hit != poly.holes_end(); ++hit) {
    for (auto vit = hit->vertices_begin(); vit != hit->vertices_end(); ++vit) {
      points.insert(*vit);
    }
  }
}

void
_collect_points(const CGAL::Triangle_3<Kernel>  &tri,
                GeometrySet<3>::PointCollection &points)
{
  points.insert(tri.vertex(0));
  points.insert(tri.vertex(1));
  points.insert(tri.vertex(2));
}

void
_collect_points(const NoVolume & /*unused*/,
                GeometrySet<2>::PointCollection & /*unused*/)
{
}

void
_collect_points(const MarkedPolyhedron          &poly,
                GeometrySet<3>::PointCollection &points)
{
  for (MarkedPolyhedron::Vertex_const_iterator vit = poly.vertices_begin();
       vit != poly.vertices_end(); ++vit) {
    points.insert(vit->point());
  }
}
/// \endcond

template <int Dim>
void
GeometrySet<Dim>::collectPoints(const PrimitiveHandle<Dim> &pa)
{
  using TPoint   = typename TypeForDimension<Dim>::Point;
  using TSegment = typename TypeForDimension<Dim>::Segment;
  using TSurface = typename TypeForDimension<Dim>::Surface;
  using TVolume  = typename TypeForDimension<Dim>::Volume;

  switch (pa.handle.which()) {
  case PrimitivePoint: {
    const TPoint *pt = boost::get<const TPoint *>(pa.handle);
    _points.insert(*pt);
    break;
  }

  case PrimitiveSegment: {
    const TSegment *seg = boost::get<const TSegment *>(pa.handle);
    _points.insert(seg->source());
    _points.insert(seg->target());
    break;
  }

  case PrimitiveSurface: {
    _collect_points(*boost::get<const TSurface *>(pa.handle), _points);
    break;
  }

  case PrimitiveVolume: {
    _collect_points(*boost::get<const TVolume *>(pa.handle), _points);
    break;
  }
  }
}

/// \cond IGNORE
template <int Dim, class IT>
void
_filter_covered(IT ibegin, IT iend, GeometrySet<Dim> &output)
{
  for (IT it = ibegin; it != iend; ++it) {
    GeometrySet<Dim> v1;
    v1.addPrimitive(it->primitive());
    bool v1_covered = false;

    for (IT it2 = it; it2 != iend; ++it2) {
      if (it == it2) {
        continue;
      }

      GeometrySet<Dim> v2;
      v2.addPrimitive(it2->primitive());

      if (algorithm::covers(v2, v1)) {
        v1_covered = true;
        break;
      }
    }

    // if its not covered by another primitive
    if (!v1_covered) {
      // and not covered by another already inserted primitive
      bool const b = algorithm::covers(output, v1);

      if (!b) {
        output.addPrimitive(it->primitive(), it->flags());
      }
    }
  }
}
/// \endcond

/**
 * Add boundary segments from 2D surface
 * @param surface The surface whose boundary to add
 */
template <>
void
GeometrySet<2>::addBoundary(const TypeForDimension<2>::Surface &surface)
{
  addSegments(surface.outer_boundary().edges_begin(),
              surface.outer_boundary().edges_end());

  for (auto hit = surface.holes_begin(); hit != surface.holes_end(); ++hit) {
    addSegments(hit->edges_begin(), hit->edges_end());
  }
}

/**
 * Add boundary from 3D surface (not implemented)
 * @param surface The surface (unused)
 */
template <>
void
GeometrySet<3>::addBoundary(
    [[maybe_unused]] const TypeForDimension<3>::Surface &surface)
{
  // TODO
}

/**
 * Get maximum geometry dimension for 2D set
 * @return Maximum dimension (2 for surfaces, 1 for segments, 0 for points, -1
 * if empty)
 */
template <>
auto
GeometrySet<2>::dimension() const -> int
{
  if (!surfaces().empty()) {
    return 2;
  }

  if (!segments().empty()) {
    return 1;
  }

  if (!points().empty()) {
    return 0;
  }

  return -1;
}

/**
 * Get maximum geometry dimension for 3D set
 * @return Maximum dimension (3 for closed volumes, 2 for surfaces, 1 for
 * segments, 0 for points, -1 if empty)
 */
template <>
auto
GeometrySet<3>::dimension() const -> int
{
  if (!volumes().empty()) {
    for (const auto &it : volumes()) {
      if (it.primitive().is_closed()) {
        return 3;
      }
    }

    return 2;
  }

  if (!surfaces().empty()) {
    return 2;
  }

  if (!segments().empty()) {
    return 1;
  }

  if (!points().empty()) {
    return 0;
  }

  return -1;
}

template <int Dim>
void
GeometrySet<Dim>::filterCovered(GeometrySet<Dim> &output) const
{
  _filter_covered(_volumes.begin(), _volumes.end(), output);
  _filter_covered(_surfaces.begin(), _surfaces.end(), output);
  _filter_covered(_segments.begin(), _segments.end(), output);
  _filter_covered(_points.begin(), _points.end(), output);
}

/**
 * @brief Stream output operator for 2D GeometrySet
 * @param ostr Output stream
 * @param geomSet The 2D geometry set to output
 * @return The output stream
 */
auto
operator<<(std::ostream &ostr, const GeometrySet<2> &geomSet) -> std::ostream &
{
  ostr << "points: ";
  std::ostream_iterator<CollectionElement<Point_d<2>::Type>> const out_pt(ostr,
                                                                          ", ");
  std::copy(geomSet.points().begin(), geomSet.points().end(), out_pt);
  ostr << '\n' << "segments: ";
  std::ostream_iterator<CollectionElement<Segment_d<2>::Type>> const out_seg(
      ostr, ", ");
  std::copy(geomSet.segments().begin(), geomSet.segments().end(), out_seg);
  ostr << '\n' << "surfaces: ";
  std::ostream_iterator<CollectionElement<Surface_d<2>::Type>> const out_surf(
      ostr, ", ");
  std::copy(geomSet.surfaces().begin(), geomSet.surfaces().end(), out_surf);
  ostr << '\n';
  return ostr;
}

/**
 * @brief Stream output operator for 3D GeometrySet
 * @param ostr Output stream
 * @param geomSet The 3D geometry set to output
 * @return The output stream
 */
auto
operator<<(std::ostream &ostr, const GeometrySet<3> &geomSet) -> std::ostream &
{
  ostr << "points: ";
  std::ostream_iterator<CollectionElement<Point_d<3>::Type>> const out_pt(ostr,
                                                                          ", ");
  std::copy(geomSet.points().begin(), geomSet.points().end(), out_pt);
  ostr << '\n' << "segments: ";
  std::ostream_iterator<CollectionElement<Segment_d<3>::Type>> const out_seg(
      ostr, ", ");
  std::copy(geomSet.segments().begin(), geomSet.segments().end(), out_seg);
  ostr << '\n' << "surfaces: ";
  std::ostream_iterator<CollectionElement<Surface_d<3>::Type>> const out_surf(
      ostr, ", ");
  std::copy(geomSet.surfaces().begin(), geomSet.surfaces().end(), out_surf);
  ostr << '\n' << "volumes: ";
  std::ostream_iterator<CollectionElement<Volume_d<3>::Type>> const out_vol(
      ostr, ", ");
  std::copy(geomSet.volumes().begin(), geomSet.volumes().end(), out_vol);
  ostr << '\n';
  return ostr;
}

template class GeometrySet<2>;
template class GeometrySet<3>;
} // namespace SFCGAL::detail
