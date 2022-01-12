// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _SFCGAL_GRAPH_GEOMETRYGRAPHBUILDER_H_
#define _SFCGAL_GRAPH_GEOMETRYGRAPHBUILDER_H_

#include <SFCGAL/LineString.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>

#include <SFCGAL/detail/graph/GeometryGraph.h>

namespace SFCGAL {
namespace graph {

/**
 * @brief [private]Convert Geometries to a GeometryGraph. Identifier in the
 * Graph are returned in order to allow the user to keep identify the geometry.
 *
 * @todo wrap vertex_descriptor, std::vector< edge_descriptor >, etc. in
 * SharedPoint, SharedLineString, SharedPolygon, etc. and add utility method on
 * the Graph?
 */
template <typename Graph>
class GeometryGraphBuilderT {
public:
  typedef Graph graph_t;

  typedef typename graph_t::vertex_properties vertex_properties;
  typedef typename graph_t::edge_properties   edge_properties;
  typedef typename graph_t::vertex_descriptor vertex_descriptor;
  typedef typename graph_t::edge_descriptor   edge_descriptor;

  /**
   * allows to match duplicates
   */
  typedef std::map<Coordinate, vertex_descriptor> coordinate_list;

  /**
   * default constructor
   */
  GeometryGraphBuilderT(graph_t &graph) : _graph(graph) {}

  /**
   * destructor
   */
  ~GeometryGraphBuilderT() {}

  /**
   * add a Point to the Graph
   */
  vertex_descriptor
  addPoint(const Point &point)
  {
    BOOST_ASSERT(!point.isEmpty());

    typename coordinate_list::const_iterator it =
        _vertices.find(point.coordinate());

    if (it != _vertices.end()) {
      return it->second;
    } else {
      vertex_descriptor vertex =
          _graph.addVertex(vertex_properties(point.coordinate()));
      _vertices.insert(std::make_pair(point.coordinate(), vertex));
      return vertex;
    }
  }

  /**
   * add a Point to the Graph
   * @return the edge inserted into the graph
   */
  edge_descriptor
  addLineSegment(const Point &a, const Point &b,
                 const edge_properties &edgeProperties = edge_properties())
  {
    BOOST_ASSERT(!a.isEmpty());
    BOOST_ASSERT(!b.isEmpty());

    return _graph.addEdge(addPoint(a), addPoint(b), edgeProperties);
  }

  /**
   * add a LineString to the graph
   * @return the list of edges inserted into the graph
   */
  std::vector<edge_descriptor>
  addLineString(const LineString      &lineString,
                const edge_properties &edgeProperties = edge_properties())
  {
    BOOST_ASSERT(!lineString.isEmpty());

    std::vector<edge_descriptor> edges;

    for (size_t i = 0; i < lineString.numPoints() - 1; i++) {
      edges.push_back(addLineSegment(lineString.pointN(i),
                                     lineString.pointN(i + 1), edgeProperties));
    }

    return edges;
  }

  /**
   * add a Triangle to the graph
   * @return the list of edges inserted into the graph
   */
  std::vector<edge_descriptor>
  addTriangle(const Triangle        &triangle,
              const edge_properties &edgeProperties = edge_properties())
  {
    BOOST_ASSERT(!triangle.isEmpty());

    std::vector<edge_descriptor> edges;

    for (size_t i = 0; i < 3; i++) {
      edges.push_back(addLineSegment(triangle.vertex(i), triangle.vertex(i + 1),
                                     edgeProperties));
    }

    return edges;
  }

  /**
   * add a Polygon to the graph
   * @returns the list of rings inserted into the graph
   */
  std::vector<std::vector<edge_descriptor>>
  addPolygon(const Polygon         &polygon,
             const edge_properties &edgeProperties = edge_properties())
  {
    BOOST_ASSERT(!polygon.isEmpty());

    std::vector<std::vector<edge_descriptor>> rings;

    for (size_t i = 0; i < polygon.numRings(); i++) {
      rings.push_back(addLineString(polygon.ringN(i), edgeProperties));
    }

    return rings;
  }

  /**
   * add a TriangulatedSurface to the graph
   * @returns the list of rings inserted into the graph
   */
  std::vector<std::vector<edge_descriptor>>
  addTriangulatedSurface(
      const TriangulatedSurface &triangulatedSurface,
      const edge_properties     &edgeProperties = edge_properties())
  {
    BOOST_ASSERT(!triangulatedSurface.isEmpty());

    std::vector<std::vector<edge_descriptor>> triangles;

    for (size_t i = 0; i < triangulatedSurface.numGeometries(); i++) {
      triangles.push_back(
          addTriangle(triangulatedSurface.geometryN(i), edgeProperties));
    }

    return triangles;
  }

  /**
   * add a PolyhedralSurface to the graph
   * @returns the list of rings inserted into the graph
   */
  std::vector<std::vector<std::vector<edge_descriptor>>>
  addPolyhedralSurface(
      const PolyhedralSurface &polyhedralSurface,
      const edge_properties   &edgeProperties = edge_properties())
  {
    BOOST_ASSERT(!polyhedralSurface.isEmpty());

    std::vector<std::vector<std::vector<edge_descriptor>>> polygons;

    for (size_t i = 0; i < polyhedralSurface.numPolygons(); i++) {
      polygons.push_back(
          addPolygon(polyhedralSurface.polygonN(i), edgeProperties));
    }

    return polygons;
  }

private:
  graph_t        &_graph;
  coordinate_list _vertices;

  /**
   * no copy constructor
   */
  GeometryGraphBuilderT(const GeometryGraphBuilderT &other);
};

typedef GeometryGraphBuilderT<GeometryGraph> GeometryGraphBuilder;

} // namespace graph
} // namespace SFCGAL

#endif
