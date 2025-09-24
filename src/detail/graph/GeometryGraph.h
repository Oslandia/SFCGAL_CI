// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_GRAPH_GEOMETRYGRAPH_H_
#define SFCGAL_GRAPH_GEOMETRYGRAPH_H_

#include <boost/graph/adjacency_list.hpp>

#include "SFCGAL/detail/graph/Edge.h"
#include "SFCGAL/detail/graph/Vertex.h"

namespace SFCGAL {
namespace graph {

/**
 * @enum EdgeDirection
 * @brief Specifies the direction of an edge in a graph.
 */
typedef enum { DIRECT = 0, REVERSE = 1 } EdgeDirection;

/**
 * @brief Returns the reverse of the given EdgeDirection.
 *
 * Converts DIRECT to REVERSE and REVERSE to DIRECT.
 *
 * @param direction The current edge direction.
 * @return The opposite edge direction.
 */
inline EdgeDirection
reverse(const EdgeDirection &direction)
{
  return (EdgeDirection)(1 - direction);
}

/**
 *
 * @brief [private]Represents the vertices and edges for a list of geometries.
 *
 * A boost::adjancency_list is wrapped in order to be able to annex some
 * information and to provide basic functionalities.
 *
 * @warning duplicate matching is performed in GeometryGraphBuilder (allows to
 * modify position once it's done)
 *
 */
template <typename VertexProperties, typename EdgeProperties>
class GeometryGraphT {
public:
  typedef VertexProperties vertex_properties; ///< Vertex properties type
  typedef EdgeProperties   edge_properties;   ///< Edge properties type

  /**
   * the wrapped graphEdgeProperties
   */
  typedef boost::adjacency_list<
      boost::listS, /* stable identifiers */
      boost::listS, /* parallel edges allowed + stable identifiers */
      boost::bidirectionalS, vertex_properties, edge_properties>
      graph_t;

  typedef typename boost::graph_traits<graph_t>::vertex_descriptor
      vertex_descriptor; ///< Vertex descriptor type from boost graph
  typedef typename boost::graph_traits<graph_t>::edge_descriptor
      edge_descriptor; ///< Edge descriptor type from boost graph

  typedef typename boost::graph_traits<graph_t>::vertex_iterator
      vertex_iterator; ///< Vertex iterator type
  typedef typename boost::graph_traits<graph_t>::edge_iterator
      edge_iterator; ///< Edge iterator type
  /**
   * An edge descriptor, with a direction.
   *
   * From the vertex point of view, out edges are DIRECT, in edges are REVERSE.
   */
  typedef std::pair<edge_descriptor, EdgeDirection> directed_edge_descriptor;

  typedef typename boost::graph_traits<graph_t>::in_edge_iterator
      in_edge_iterator; ///< Incoming edge iterator type
  typedef typename boost::graph_traits<graph_t>::out_edge_iterator
      out_edge_iterator; ///< Outgoing edge iterator type

  /**
   * @brief Returns the number of vertices
   * @return The number of vertices in the graph
   */
  inline size_t
  numVertices() const
  {
    return boost::num_vertices(_graph);
  }

  /**
   * @brief Return vertex iterator pair
   * @return Pair of vertex iterators (begin, end)
   *
   * @code
   * typename GeometryGraph<V,E>::vertex_iterator it, end ;
   * for ( boost::tie( it, end ) = g.vertices(); it != end; ++it ){
   * 		typename GeometryGraph<V,E>::vertex_descriptor vertex = *it ;
   * 		std::cout << g[ vertex ].coordinate << std::endl;
   * }
   * @endcode
   */
  inline std::pair<vertex_iterator, vertex_iterator>
  vertices() const
  {
    return boost::vertices(_graph);
  }

  /**
   * @brief Add a vertex to the graph
   * @param properties The properties for the new vertex
   * @return The identifier of the vertex
   */
  vertex_descriptor
  addVertex(const vertex_properties &properties = vertex_properties())
  {
    return boost::add_vertex(properties, _graph);
  }
  /**
   * @brief Remove a vertex (and all its adjacent edges)
   * @param vertex The vertex descriptor to remove
   */
  void
  removeVertex(const vertex_descriptor &vertex)
  {
    boost::clear_vertex(vertex);
    boost::remove_vertex(vertex);
  }

  /**
   * @brief Returns the number of edges
   * @return The number of edges in the graph
   */
  inline size_t
  numEdges() const
  {
    return boost::num_edges(_graph);
  }

  /**
   * @brief Return edge iterator pair
   * @return Pair of edge iterators (begin, end)
   *
   * @code
   * typename GeometryGraph<V,E>::edge_iterator it, end ;
   * for ( boost::tie( it, end ) = g.edges(); it != end; ++it ){
   * 		typename GeometryGraph<V,E>::edge_descriptor edge = *it ;
   * 		std::cout << g.source(edge) << "," << g.target(edge) << std::endl;
   * }
   * @endcode
   */
  inline std::pair<edge_iterator, edge_iterator>
  edges() const
  {
    return boost::edges(_graph);
  }

  /**
   * @brief Add an Edge to the Graph
   * @param source The source vertex descriptor
   * @param target The target vertex descriptor
   * @param properties The properties for the new edge
   * @return The identifier of the edge
   */
  edge_descriptor
  addEdge(const vertex_descriptor &source, const vertex_descriptor &target,
          const EdgeProperties &properties = EdgeProperties())
  {
    return boost::add_edge(source, target, properties, _graph).first;
  }

  /**
   * @brief Get the source vertex for an edge
   * @param edge The edge descriptor
   * @return The source vertex descriptor
   */
  vertex_descriptor
  source(const edge_descriptor &edge) const
  {
    return boost::source(edge, _graph);
  }
  /**
   * @brief Get the source vertex for an edge with direction
   * @param edge The edge descriptor
   * @param direction The edge direction
   * @return The source vertex descriptor
   */
  vertex_descriptor
  source(const edge_descriptor &edge, const EdgeDirection &direction) const
  {
    return direction == DIRECT ? boost::source(edge, _graph)
                               : boost::target(edge, _graph);
  }
  /**
   * @brief Get the source vertex for a directed edge
   * @param edge The directed edge descriptor
   * @return The source vertex descriptor
   */
  vertex_descriptor
  source(const directed_edge_descriptor &edge) const
  {
    return source(edge.first, edge.second);
  }

  /**
   * @brief Get the target vertex for an edge
   * @param edge The edge descriptor
   * @return The target vertex descriptor
   */
  vertex_descriptor
  target(const edge_descriptor &edge) const
  {
    return boost::target(edge, _graph);
  }
  /**
   * @brief Get the target vertex for an edge with direction
   * @param edge The edge descriptor
   * @param direction The edge direction
   * @return The target vertex descriptor
   */
  vertex_descriptor
  target(const edge_descriptor &edge, const EdgeDirection &direction) const
  {
    return direction == DIRECT ? boost::target(edge, _graph)
                               : boost::source(edge, _graph);
  }
  /**
   * @brief Get the target vertex for a directed edge
   * @param edge The directed edge descriptor
   * @return The target vertex descriptor
   */
  vertex_descriptor
  target(const directed_edge_descriptor &edge) const
  {
    return target(edge.first, edge.second);
  }
  /**
   * @brief Remove an edge
   * @param edge The edge descriptor to remove
   */
  void
  removeEdge(const edge_descriptor &edge)
  {
    boost::remove_edge(edge, _graph);
  }

  /**
   * @brief Get edges from a to b and from b to a
   * @param a First vertex descriptor
   * @param b Second vertex descriptor
   * @return Vector of directed edge descriptors connecting the vertices
   */
  std::vector<directed_edge_descriptor>
  edges(const vertex_descriptor &a, const vertex_descriptor &b) const
  {
    std::vector<directed_edge_descriptor> result;

    // out_edges from a targeting b
    {
      out_edge_iterator it, end;

      for (boost::tie(it, end) = boost::out_edges(a, _graph); it != end; ++it) {
        if (target(*it) != b) {
          continue;
        }

        result.push_back(std::make_pair(*it, DIRECT));
      }
    }
    // out_edges from b targeting a
    {
      out_edge_iterator it, end;

      for (boost::tie(it, end) = boost::out_edges(b, _graph); it != end; ++it) {
        if (target(*it) != a) {
          continue;
        }

        result.push_back(std::make_pair(*it, REVERSE));
      }
    }
    return result;
  }

  /**
   * @brief Returns the degree of a vertex
   * @param vertex The vertex descriptor
   * @return The degree of the vertex
   */
  inline size_t
  degree(const vertex_descriptor &vertex) const
  {
    return boost::degree(vertex, _graph);
  }

  /**
   * @brief Get in edges for a vertex
   * @param vertex The vertex descriptor
   * @return Vector of incoming edge descriptors
   */
  std::vector<edge_descriptor>
  inEdges(const vertex_descriptor &vertex)
  {
    std::vector<edge_descriptor> edges;

    in_edge_iterator it, end;

    for (boost::tie(it, end) = boost::in_edges(vertex, _graph); it != end;
         ++it) {
      edges.push_back(*it);
    }

    return edges;
  }
  /**
   * @brief Get out edges for a vertex
   * @param vertex The vertex descriptor
   * @return Vector of outgoing edge descriptors
   */
  std::vector<edge_descriptor>
  outEdges(const vertex_descriptor &vertex) const
  {
    std::vector<edge_descriptor> edges;

    out_edge_iterator it, end;

    for (boost::tie(it, end) = boost::out_edges(vertex, _graph); it != end;
         ++it) {
      edges.push_back(*it);
    }

    return edges;
  }
  /**
   * @brief Get in/out edges for a vertex
   * @param vertex The vertex descriptor
   * @return Vector of directed edge descriptors (in and out edges)
   */
  std::vector<directed_edge_descriptor>
  inOutEdges(const vertex_descriptor &vertex) const
  {
    std::vector<directed_edge_descriptor> edges;
    {
      in_edge_iterator it, end;

      for (boost::tie(it, end) = boost::in_edges(vertex, _graph); it != end;
           ++it) {
        edges.push_back(std::make_pair(*it, REVERSE));
      }
    }

    // out edges
    {
      out_edge_iterator it, end;

      for (boost::tie(it, end) = boost::out_edges(vertex, _graph); it != end;
           ++it) {
        edges.push_back(std::make_pair(*it, DIRECT));
      }
    }
    return edges;
  }

  /**
   * @brief Returns the list of the adjacent vertices using both DIRECT and
   * REVERSE direction
   * @param vertex Input vertex
   * @param withReverseDirection Indicates if in_edges are used
   * @return Set of adjacent vertex descriptors
   */
  std::set<vertex_descriptor>
  adjacentVertices(const vertex_descriptor &vertex,
                   bool                     withReverseDirection = true)
  {
    std::set<vertex_descriptor> vertices;
    // out edges
    {
      out_edge_iterator it, end;

      for (boost::tie(it, end) = boost::out_edges(vertex, _graph); it != end;
           ++it) {
        vertex_descriptor reached = target(*it);

        if (reached != vertex) {
          vertices.insert(reached);
        }
      }
    }

    // in edges
    if (withReverseDirection) {
      in_edge_iterator it, end;

      for (boost::tie(it, end) = boost::in_edges(vertex, _graph); it != end;
           ++it) {
        vertex_descriptor reached = source(*it);

        if (reached != vertex) {
          vertices.insert(reached);
        }
      }
    }

    return vertices;
  }

  /**
   * @brief Indicates if edges are opposite
   * @param a First edge descriptor
   * @param b Second edge descriptor
   * @return True if edges are opposite
   */
  inline bool
  areOpposite(const edge_descriptor &a, const edge_descriptor &b) const
  {
    return source(a) == target(b) && target(a) == source(b);
  }
  /**
   * @brief Indicates if edges are parallel
   * @param a First edge descriptor
   * @param b Second edge descriptor
   * @return True if edges are parallel
   */
  inline bool
  areParallel(const edge_descriptor &a, const edge_descriptor &b) const
  {
    return source(a) == source(b) && target(a) == target(b);
  }

  /**
   * @brief Revert the order of a list of edges. Old edges are removed from
   * the graph, new ones are created.
   * @param edges Vector of edge descriptors to reverse
   *
   * @warning properties are kept but oriented one (left face, right face, etc.)
   * are lost.
   */
  void
  reverse(std::vector<edge_descriptor> &edges)
  {
    std::vector<edge_descriptor> result;

    for (typename std::vector<edge_descriptor>::reverse_iterator it =
             edges.rbegin();
         it != edges.rend(); ++it) {
      edge_descriptor newEdge = addEdge(target(*it), source(*it), (*this)[*it]);
      result.push_back(newEdge);
      removeEdge(*it);
    }

    edges = result;
  }

  /**
   * @brief Returns the VertexProperties attached to a Vertex
   * @param vertex The vertex descriptor
   * @return Const reference to vertex properties
   */
  inline const vertex_properties &
  operator[](const vertex_descriptor &vertex) const
  {
    return _graph[vertex];
  }
  /**
   * @brief Returns the VertexProperties attached to a Vertex
   * @param vertex The vertex descriptor
   * @return Reference to vertex properties
   */
  inline vertex_properties &
  operator[](const vertex_descriptor &vertex)
  {
    return _graph[vertex];
  }

  /**
   * @brief Returns the EdgeProperties attached to an Edge
   * @param edge The edge descriptor
   * @return Const reference to edge properties
   */
  inline const edge_properties &
  operator[](const edge_descriptor &edge) const
  {
    return _graph[edge];
  }
  /**
   * @brief Returns the EdgeProperties attached to an Edge
   * @param edge The edge descriptor
   * @return Reference to edge properties
   */
  inline edge_properties &
  operator[](const edge_descriptor &edge)
  {
    return _graph[edge];
  }

  /**
   * @brief Returns the wrapped boost::graph
   * @return Reference to the underlying graph
   */
  inline graph_t &
  graph()
  {
    return _graph;
  }
  /**
   * @brief Returns the wrapped boost::graph
   * @return Const reference to the underlying graph
   */
  inline const graph_t &
  graph() const
  {
    return _graph;
  }

  /**
   * implicit cast to the wrapped boost graph in order to keep boost graph
   * interface
   */
  operator graph_t &(void) { return _graph; }
  /**
   * implicit cast to the wrapped boost graph in order to keep boost graph
   * interface
   */
  operator const graph_t &(void) const { return _graph; }

private:
  /**
   * a wrapped boost::graph
   */
  graph_t _graph;
};

/**
 * Default GeometryGraph with predefined Vertex and Edge properties for general
 * usage
 */
typedef GeometryGraphT<Vertex, Edge> GeometryGraph;

} // namespace graph
} // namespace SFCGAL

#endif
