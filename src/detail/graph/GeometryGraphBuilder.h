// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_GRAPH_GEOMETRYGRAPHBUILDER_H_
#define SFCGAL_GRAPH_GEOMETRYGRAPHBUILDER_H_

#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"

#include "SFCGAL/detail/graph/GeometryGraph.h"

namespace SFCGAL::graph {

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
  using graph_t = Graph; ///< Graph type being used

  typedef typename graph_t::vertex_properties
      vertex_properties; ///< Vertex properties type
  typedef typename graph_t::edge_properties
      edge_properties; ///< Edge properties type
  typedef typename graph_t::vertex_descriptor
      vertex_descriptor; ///< Vertex descriptor type
  typedef typename graph_t::edge_descriptor
      edge_descriptor; ///< Edge descriptor type

  /**
   * allows to match duplicates
   */
  using coordinate_list = std::map<Coordinate, vertex_descriptor>;

  /**
   * @brief Default constructor
   * @param graph The graph to build into
   */
  GeometryGraphBuilderT(graph_t &graph) : _graph(graph) {}

  /**
   * destructor
   */
  ~GeometryGraphBuilderT() {}

  /**
   * @brief Add a Point to the Graph
   * @param point The point to add
   * @return The vertex descriptor for the added point
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
   * @brief Add a line segment to the Graph
   * @param a First point of the segment
   * @param b Second point of the segment
   * @param edgeProperties Properties for the edge
   * @return The edge inserted into the graph
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
   * @brief Add a LineString to the graph
   * @param lineString The linestring to add
   * @param edgeProperties Properties for the edges
   * @return The list of edges inserted into the graph
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
   * @brief Add a Triangle to the graph
   * @param triangle The triangle to add
   * @param edgeProperties Properties for the edges
   * @return The list of edges inserted into the graph
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
   * @brief Add a Polygon to the graph
   * @param polygon The polygon to add
   * @param edgeProperties Properties for the edges
   * @return The list of rings inserted into the graph
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
   * @brief Add a TriangulatedSurface to the graph
   * @param triangulatedSurface The triangulated surface to add
   * @param edgeProperties Properties for the edges
   * @return The list of triangles inserted into the graph
   */
  std::vector<std::vector<edge_descriptor>>
  addTriangulatedSurface(
      const TriangulatedSurface &triangulatedSurface,
      const edge_properties     &edgeProperties = edge_properties())
  {
    BOOST_ASSERT(!triangulatedSurface.isEmpty());

    std::vector<std::vector<edge_descriptor>> triangles;

    for (size_t i = 0; i < triangulatedSurface.numPatches(); i++) {
      triangles.push_back(
          addTriangle(triangulatedSurface.patchN(i), edgeProperties));
    }

    return triangles;
  }

  /**
   * @brief Add a PolyhedralSurface to the graph
   * @param polyhedralSurface The polyhedral surface to add
   * @param edgeProperties Properties for the edges
   * @return The list of polygons inserted into the graph
   */
  std::vector<std::vector<std::vector<edge_descriptor>>>
  addPolyhedralSurface(
      const PolyhedralSurface &polyhedralSurface,
      const edge_properties   &edgeProperties = edge_properties())
  {
    BOOST_ASSERT(!polyhedralSurface.isEmpty());

    std::vector<std::vector<std::vector<edge_descriptor>>> polygons;

    for (size_t i = 0; i < polyhedralSurface.numPatches(); i++) {
      polygons.push_back(
          addPolygon(polyhedralSurface.patchN(i), edgeProperties));
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

typedef GeometryGraphBuilderT<GeometryGraph>
    GeometryGraphBuilder; ///< Geometry graph builder type

} // namespace SFCGAL::graph

#endif
