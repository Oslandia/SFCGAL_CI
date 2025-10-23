// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/isValid.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/algorithm/connection.h"
#include "SFCGAL/algorithm/distance.h"
#include "SFCGAL/algorithm/distance3d.h"
#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/algorithm/intersects.h"
#include "SFCGAL/algorithm/length.h"
#include "SFCGAL/algorithm/normal.h"
#include "SFCGAL/algorithm/orientation.h"
#include "SFCGAL/algorithm/plane.h"
#include "SFCGAL/detail/ForceValidityVisitor.h"
#include "SFCGAL/detail/GetPointsVisitor.h"
#include "SFCGAL/detail/algorithm/coversPoints.h"
#include "SFCGAL/detail/tools/Log.h"

#include <boost/format.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/undirected_dfs.hpp>
#include <boost/graph/visitors.hpp>

using namespace SFCGAL::detail::algorithm;

namespace SFCGAL {

/// @{
/// @privatesection
#ifndef DOXYGEN_SHOULD_SKIP_THIS

void
SFCGAL_ASSERT_GEOMETRY_VALIDITY_(const Geometry &g, const std::string &ctxt)
{
  if (!(g).hasValidityFlag()) {
    const Validity sfcgalAssertGeometryValidity = algorithm::isValid(g);
    if (!sfcgalAssertGeometryValidity) {
      throw GeometryInvalidityException(
          (boost::format(ctxt + "%s is invalid : %s : %s") %
           (g).geometryType() % sfcgalAssertGeometryValidity.reason() %
           (g).asText())
              .str());
    }
  }
}

#endif // ifndef DOXYGEN_SHOULD_SKIP_THIS
/// @} end of private section

void
SFCGAL_ASSERT_GEOMETRY_VALIDITY(const Geometry &g)
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_(g, "");
}

void
SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(const Geometry &g)
{
  if (!(g).hasValidityFlag()) {
    using namespace SFCGAL;
    if ((g).is3D()) {
      std::unique_ptr<SFCGAL::Geometry> sfcgalAssertGeometryValidityClone(
          (g).clone());
      algorithm::force2D(*sfcgalAssertGeometryValidityClone);
      SFCGAL_ASSERT_GEOMETRY_VALIDITY_((*sfcgalAssertGeometryValidityClone),
                                       "When converting to 2D - ");
    } else {
      SFCGAL_ASSERT_GEOMETRY_VALIDITY(g);
    }
  }
}

void
SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(const Geometry &g)
{
  if (!(g).hasValidityFlag()) {
    using namespace SFCGAL;
    if (!(g).is3D()) {
      std::unique_ptr<Geometry> sfcgalAssertGeometryValidityClone = g.clone();
      algorithm::force3D(*sfcgalAssertGeometryValidityClone);
      SFCGAL_ASSERT_GEOMETRY_VALIDITY_((*sfcgalAssertGeometryValidityClone),
                                       "When converting to 3D - ");
    } else {
      SFCGAL_ASSERT_GEOMETRY_VALIDITY(g);
    }
  }
}

void
SFCGAL_ASSERT_GEOMETRY_VALIDITY_ON_PLANE([[maybe_unused]] const Geometry &geom)
{
  throw NotImplementedException(
      "validation on geometry projected on arbitrary plane is not implemented");
}

namespace algorithm {

/// @{
/// @privatesection
#ifndef DOXYGEN_SHOULD_SKIP_THIS

struct LoopDetector : public boost::dfs_visitor<> {
  /** Constructor
   * @param hasLoop Reference to boolean flag to set when loop detected */
  LoopDetector(bool &hasLoop) : _hasLoop(hasLoop) {}

  /** Callback for back edge detection indicating a loop */
  template <class Edge, class Graph>
  void
  back_edge(Edge /*edge*/, const Graph & /*graph*/)
  {
    _hasLoop = true;
  }

private:
  bool &_hasLoop;
};

// Forward declarations for polygon validation helper functions
namespace detail {
auto
validatePolygonRingsBasic(const Polygon &polygon) -> Validity;
auto
validatePolygonRingsOrientation(const Polygon &polygon,
                                const double  &toleranceAbs) -> Validity;
auto
validatePolygonRingsIntersections(const Polygon &polygon) -> Validity;
auto
validatePolygonInteriorRings(const Polygon &polygon) -> Validity;
} // namespace detail

auto
isValid(const Point &point) -> Validity
{
  if (point.isEmpty()) {
    return Validity::valid();
  }

  (void)point;
  // return( isValid(point.coordinate() ) );
  return Validity::valid();
}

auto
isValid(const LineString &linestring, const double &toleranceAbs) -> Validity
{
  if (linestring.isEmpty()) {
    return Validity::valid();
  }

  return length3D(linestring) > toleranceAbs ? Validity::valid()
                                             : Validity::invalid("no length");
}

auto
isValid(const Polygon &polygon, const double &toleranceAbs) -> Validity
{
  if (polygon.isEmpty()) {
    return Validity::valid();
  }

  // Validate basic ring properties
  auto basicValidation = detail::validatePolygonRingsBasic(polygon);
  if (!basicValidation) {
    return basicValidation;
  }

  // Validate ring orientations
  auto orientationValidation =
      detail::validatePolygonRingsOrientation(polygon, toleranceAbs);
  if (!orientationValidation) {
    return orientationValidation;
  }

  // Validate ring intersections
  auto intersectionValidation =
      detail::validatePolygonRingsIntersections(polygon);
  if (!intersectionValidation) {
    return intersectionValidation;
  }

  // Validate interior rings
  auto interiorValidation = detail::validatePolygonInteriorRings(polygon);
  if (!interiorValidation) {
    return interiorValidation;
  }

  return Validity::valid();
}

namespace detail {

auto
validatePolygonRingsBasic(const Polygon &polygon) -> Validity
{
  const size_t numRings = polygon.numRings();

  for (size_t ring = 0; ring != numRings; ++ring) {
    if (polygon.ringN(ring).numPoints() < 4) {
      return Validity::invalid(
          (boost::format("not enough points in ring %d") % ring).str());
    }

    const double distanceToClose =
        polygon.is3D() ? distancePointPoint3D(polygon.ringN(ring).startPoint(),
                                              polygon.ringN(ring).endPoint())
                       : distancePointPoint(polygon.ringN(ring).startPoint(),
                                            polygon.ringN(ring).endPoint());

    if (distanceToClose > 0) {
      return Validity::invalid(
          (boost::format("ring %d is not closed") % ring).str());
    }

    if (polygon.is3D() ? selfIntersects3D(polygon.ringN(ring))
                       : selfIntersects(polygon.ringN(ring))) {
      return Validity::invalid(
          (boost::format("ring %d self intersects") % ring).str());
    }
  }

  // Check for degenerate rings (all points identical)
  for (size_t ring = 0; ring != numRings; ++ring) {
    const LineString &currentRing = polygon.ringN(ring);
    const Point      &start       = currentRing.startPoint();
    size_t            idx         = 0;
    for (; idx < currentRing.numPoints() && start == currentRing.pointN(idx);
         idx++) {
      ; // noop
    }
    if (idx == currentRing.numPoints()) {
      return Validity::invalid(
          (boost::format("ring %d degenerated to a point") % ring).str());
    }
  }

  return Validity::valid();
}

auto
validatePolygonRingsOrientation(const Polygon &polygon,
                                const double  &toleranceAbs) -> Validity
{
  // Orientation in 2D
  if (!polygon.is3D()) {
    const bool extCCWO = isCounterClockWiseOriented(polygon.exteriorRing());

    for (std::size_t ring = 0; ring < polygon.numInteriorRings(); ++ring) {
      if (extCCWO == isCounterClockWiseOriented(polygon.interiorRingN(ring))) {
        return Validity::invalid(
            (boost::format("exterior ring and interior ring %d have the same "
                           "orientation") %
             ring)
                .str());
      }
    }
  }
  // Orientation in 3D
  else {
    // Polygon must be planar
    if (!isPlane3D<Kernel>(polygon, toleranceAbs)) {
      return Validity::invalid("points don't lie in the same plane");
    }

    // Interior rings must be oriented opposite to exterior
    if (polygon.hasInteriorRings()) {
      const CGAL::Vector_3<Kernel> nExt =
          normal3D<Kernel>(polygon.exteriorRing());

      for (std::size_t ring = 0; ring < polygon.numInteriorRings(); ++ring) {
        const CGAL::Vector_3<Kernel> nInt =
            normal3D<Kernel>(polygon.interiorRingN(ring));

        if (nExt * nInt > 0) {
          return Validity::invalid(
              (boost::format("interior ring %d is oriented in the same "
                             "direction as exterior ring") %
               ring)
                  .str());
        }
      }
    }
  }

  return Validity::valid();
}

auto
validatePolygonRingsIntersections(const Polygon &polygon) -> Validity
{
  const size_t numRings = polygon.numRings();

  using Edge = std::pair<int, int>;
  std::vector<Edge> touchingRings;

  for (size_t ringIdx = 0; ringIdx < numRings; ++ringIdx) {
    for (size_t otherRingIdx = ringIdx + 1; otherRingIdx < numRings;
         ++otherRingIdx) {
      std::unique_ptr<Geometry> inter =
          polygon.is3D() ? intersection3D(polygon.ringN(ringIdx),
                                          polygon.ringN(otherRingIdx))
                         : intersection(polygon.ringN(ringIdx),
                                        polygon.ringN(otherRingIdx));

      if (!inter->isEmpty() && !inter->is<Point>()) {
        return Validity::invalid(
            (boost::format("intersection between ring %d and %d") % ringIdx %
             otherRingIdx)
                .str());
      }
      if (!inter->isEmpty() && inter->is<Point>()) {
        touchingRings.emplace_back(ringIdx, otherRingIdx);
      }
    }
  }

  // Check for unconnected interior using graph theory
  using namespace boost;
  using Graph    = adjacency_list<vecS, vecS, undirectedS, no_property,
                                  property<edge_color_t, default_color_type>>;
  using vertex_t = graph_traits<Graph>::vertex_descriptor;

  Graph graph(touchingRings.begin(), touchingRings.end(), numRings);

  bool               hasLoop = false;
  LoopDetector const vis(hasLoop);
  undirected_dfs(graph, root_vertex(vertex_t(0))
                            .visitor(vis)
                            .edge_color_map(get(edge_color, graph)));

  if (hasLoop) {
    return Validity::invalid("interior is not connected");
  }

  return Validity::valid();
}

auto
validatePolygonInteriorRings(const Polygon &polygon) -> Validity
{
  if (!polygon.hasInteriorRings()) {
    return Validity::valid();
  }

  // Interior rings must be interior to exterior ring
  for (size_t ring = 0; ring < polygon.numInteriorRings(); ++ring) {
    if (polygon.is3D() ? !coversPoints3D(Polygon(polygon.exteriorRing()),
                                         Polygon(polygon.interiorRingN(ring)))
                       : !coversPoints(Polygon(polygon.exteriorRing()),
                                       Polygon(polygon.interiorRingN(ring)))) {
      return Validity::invalid(
          (boost::format("exterior ring doesn't cover interior ring %d") % ring)
              .str());
    }
  }

  // Interior rings must not cover one another
  for (size_t ri = 0; ri < polygon.numInteriorRings(); ++ri) {
    for (size_t rj = ri + 1; rj < polygon.numInteriorRings(); ++rj) {
      if (polygon.is3D() ? coversPoints3D(Polygon(polygon.interiorRingN(ri)),
                                          Polygon(polygon.interiorRingN(rj)))
                         : coversPoints(Polygon(polygon.interiorRingN(ri)),
                                        Polygon(polygon.interiorRingN(rj)))) {
        return Validity::invalid(
            (boost::format("interior ring %d covers interior ring %d") % ri %
             rj)
                .str());
      }
    }
  }

  return Validity::valid();
}

} // namespace detail

auto
isValid(const Triangle &triangle, const double &toleranceAbs) -> Validity
{
  return isValid(triangle.toPolygon(), toleranceAbs);
}

auto
isValid(const MultiLineString &multilinestring, const double &toleranceAbs)
    -> Validity
{
  if (multilinestring.isEmpty()) {
    return Validity::valid();
  }

  const size_t numLineString = multilinestring.numGeometries();

  for (size_t l = 0; l != numLineString; ++l) {
    Validity const validity =
        isValid(multilinestring.lineStringN(l), toleranceAbs);

    if (!validity) {
      return Validity::invalid((boost::format("LineString %d is invalid: %s") %
                                l % validity.reason())
                                   .str());
    }
  }

  return Validity::valid();
}

auto
isValid(const MultiPolygon &multipolygon, const double &toleranceAbs)
    -> Validity
{
  if (multipolygon.isEmpty()) {
    return Validity::valid();
  }

  const size_t numPolygons = multipolygon.numGeometries();

  for (size_t poly = 0; poly != numPolygons; ++poly) {
    Validity const validity =
        isValid(multipolygon.polygonN(poly), toleranceAbs);

    if (!validity) {
      return Validity::invalid((boost::format("Polygon %d is invalid: %s") %
                                poly % validity.reason())
                                   .str());
    }
  }

  for (size_t pi = 0; pi != numPolygons; ++pi) {
    for (size_t pj = pi + 1; pj < numPolygons; ++pj) {
      std::unique_ptr<Geometry> inter =
          multipolygon.is3D() ? intersection3D(multipolygon.polygonN(pi),
                                               multipolygon.polygonN(pj))
                              : intersection(multipolygon.polygonN(pi),
                                             multipolygon.polygonN(pj));

      // intersection can be empty, a point, or a set of points
      if (!inter->isEmpty() && inter->dimension() != 0) {
        return Validity::invalid(
            (boost::format("intersection between Polygon %d and %d") % pi % pj)
                .str());
      }
    }
  }

  return Validity::valid();
}

auto
isValid(const GeometryCollection &geometrycollection,
        const double             &toleranceAbs) -> Validity
{
  if (geometrycollection.isEmpty()) {
    return Validity::valid();
  }

  const size_t numGeom = geometrycollection.numGeometries();

  for (size_t geom = 0; geom != numGeom; ++geom) {
    Validity const validity =
        isValid(geometrycollection.geometryN(geom), toleranceAbs);

    if (!validity) {
      return Validity::invalid(
          (boost::format("%s %d is invalid: %s") %
           geometrycollection.geometryN(geom).geometryType() % geom %
           validity.reason())
              .str());
    }
  }

  return Validity::valid();
}

auto
isValid(const TriangulatedSurface &triangulatedsurface,
        const SurfaceGraph &graph, const double &toleranceAbs) -> Validity
{
  if (triangulatedsurface.isEmpty()) {
    return Validity::valid();
  }

  size_t const numPatches = triangulatedsurface.numPatches();

  for (size_t tri = 0; tri != numPatches; ++tri) {
    Validity const validity =
        isValid(triangulatedsurface.patchN(tri), toleranceAbs);

    if (!validity) {
      return Validity::invalid((boost::format("Triangle %d is invalid: %s") %
                                tri % validity.reason())
                                   .str());
    }
  }

  if (!isConnected(graph)) {
    return Validity::invalid("not connected");
  }

  if (triangulatedsurface.is3D() ? selfIntersects3D(triangulatedsurface, graph)
                                 : selfIntersects(triangulatedsurface, graph)) {
    return Validity::invalid("self intersects");
  }

  return Validity::valid();
}

auto
isValid(const TriangulatedSurface &triangulatedsurface,
        const double              &toleranceAbs) -> Validity
{
  if (triangulatedsurface.isEmpty()) {
    return Validity::valid();
  }

  const SurfaceGraph graph(triangulatedsurface);
  return graph.isValid() ? isValid(triangulatedsurface, graph, toleranceAbs)
                         : graph.isValid();
}

auto
isValid(const PolyhedralSurface &polyhedralsurface, const SurfaceGraph &graph,
        const double &toleranceAbs) -> Validity
{
  if (polyhedralsurface.isEmpty()) {
    return Validity::valid();
  }

  size_t const numPatches = polyhedralsurface.numPatches();

  for (size_t patch = 0; patch != numPatches; ++patch) {
    Validity const validity =
        isValid(polyhedralsurface.patchN(patch), toleranceAbs);

    if (!validity) {
      return Validity::invalid((boost::format("Polygon %d is invalid: %s") %
                                patch % validity.reason())
                                   .str());
    }
  }

  if (!isConnected(graph)) {
    return Validity::invalid("not connected");
  }

  if (polyhedralsurface.is3D() ? selfIntersects3D(polyhedralsurface, graph)
                               : selfIntersects(polyhedralsurface, graph)) {
    return Validity::invalid("self intersects");
  }

  return Validity::valid();
}

auto
isValid(const PolyhedralSurface &polyhedralsurface, const double &toleranceAbs)
    -> Validity
{
  if (polyhedralsurface.isEmpty()) {
    return Validity::valid();
  }

  const SurfaceGraph graph(polyhedralsurface);
  return graph.isValid() ? isValid(polyhedralsurface, graph, toleranceAbs)
                         : graph.isValid();
}

auto
isValid(const Solid &solid, const double &toleranceAbs) -> Validity
{
  if (solid.isEmpty()) {
    return Validity::valid();
  }

  const size_t numShells = solid.numShells();

  for (size_t shell = 0; shell != numShells; ++shell) {
    const SurfaceGraph graph(solid.shellN(shell));
    Validity const validity = isValid(solid.shellN(shell), graph, toleranceAbs);

    if (!validity) {
      return Validity::invalid(
          (boost::format("PolyhedralSurface (shell) %d is invalid: %s") %
           shell % validity.reason())
              .str());
    }

    if (!isClosed(graph)) {
      return Validity::invalid(
          (boost::format("PolyhedralSurface (shell) %d is not closed") % shell)
              .str());
    }
  }

  if (solid.numInteriorShells() != 0U) {
    BOOST_THROW_EXCEPTION(
        Exception("function is not fully implemented (orientation, covering "
                  "and intersections of interior shells missing"));
  }

  return Validity::valid();
}

auto
isValid(const MultiSolid &multisolid, const double &toleranceAbs) -> Validity
{
  if (multisolid.isEmpty()) {
    return Validity::valid();
  }

  const size_t numMultiSolid = multisolid.numGeometries();

  for (size_t solid = 0; solid != numMultiSolid; ++solid) {
    Validity const validity = isValid(multisolid.solidN(solid), toleranceAbs);

    if (!validity) {
      return Validity::invalid(
          (boost::format("Solid %d is invalid: %s") % solid % validity.reason())
              .str());
    }
  }

  return Validity::valid();
}

auto
isValid(const NURBSCurve              &nurbsCurve,
        [[maybe_unused]] const double &toleranceAbs) -> Validity
{
  auto [valid, reason] = nurbsCurve.validateData();
  return valid ? Validity::valid() : Validity::invalid(reason);
}

#endif // ifndef DOXYGEN_SHOULD_SKIP_THIS
/// @} end of private section

auto
isValid(const Geometry &geometry, const double &toleranceAbs) -> Validity
{
  switch (geometry.geometryTypeId()) {
  case TYPE_POINT:
    return isValid(geometry.as<Point>());

  case TYPE_LINESTRING:
    return isValid(geometry.as<LineString>(), toleranceAbs);

  case TYPE_POLYGON:
    return isValid(geometry.as<Polygon>(), toleranceAbs);

  case TYPE_TRIANGLE:
    return isValid(geometry.as<Triangle>(), toleranceAbs);

  case TYPE_SOLID:
    return isValid(geometry.as<Solid>(), toleranceAbs);

  case TYPE_MULTIPOINT:
    return Validity::valid();

  case TYPE_MULTILINESTRING:
    return isValid(geometry.as<MultiLineString>(), toleranceAbs);

  case TYPE_MULTIPOLYGON:
    return isValid(geometry.as<MultiPolygon>(), toleranceAbs);

  case TYPE_MULTISOLID:
    return isValid(geometry.as<MultiSolid>(), toleranceAbs);

  case TYPE_GEOMETRYCOLLECTION:
    return isValid(geometry.as<GeometryCollection>(), toleranceAbs);

  case TYPE_TRIANGULATEDSURFACE:
    return isValid(geometry.as<TriangulatedSurface>(), toleranceAbs);

  case TYPE_POLYHEDRALSURFACE:
    return isValid(geometry.as<PolyhedralSurface>(), toleranceAbs);

  case TYPE_NURBSCURVE:
    return isValid(geometry.as<NURBSCurve>(), toleranceAbs);
  }

  BOOST_THROW_EXCEPTION(Exception(
      (boost::format("isValid( %s ) is not defined") % geometry.geometryType())
          .str()));
  return Validity::invalid(
      (boost::format("isValid( %s ) is not defined") % geometry.geometryType())
          .str()); // to avoid warning
}

void
propagateValidityFlag(Geometry &geometry, bool valid)
{
  SFCGAL::detail::ForceValidityVisitor visitor(valid);
  geometry.accept(visitor);
}

} // namespace algorithm
} // namespace SFCGAL
