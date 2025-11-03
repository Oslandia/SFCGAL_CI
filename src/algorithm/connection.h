// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_GRAPH_ISCONNECTED_H_
#define SFCGAL_GRAPH_ISCONNECTED_H_

#include "SFCGAL/Coordinate.h"
#include "SFCGAL/Geometry.h"
#include "SFCGAL/Validity.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/noncopyable.hpp>
#include <map>

namespace SFCGAL::algorithm {

/**
 * Represents a polyhedral surface as a graph where faces are nodes and egde are
 * graph edges
 * @pre the polygons are valid
 * @todo unittest
 */

class SFCGAL_API SurfaceGraph : boost::noncopyable {
public:
  /// Vertex index type for identifying vertices in the surface graph
  using VertexIndex = size_t;
  /// Face index type for identifying faces in the surface graph
  using FaceIndex = size_t;
  /// Map type for associating coordinates with vertex indices
  using CoordinateMap = std::map<Coordinate, VertexIndex>;
  /// Constant representing an invalid index value
  static const size_t INVALID_INDEX;
  // an edge is inserted with vtx ordered by the first polygon we treat,
  // we search the edge with reverse ordered vtx indexes.
  // as a result, an inconsistent orientation between polygons can be spotted by
  // finding the edge in the same order
  // note that this situation may be caused if a face is duplicated
  /// Map type for storing edges and their adjacent faces
  using EdgeMap = std::map<std::pair<VertexIndex, VertexIndex>,
                           std::pair<FaceIndex, FaceIndex>>;
  /// Graph type representing face adjacency relationships
  using FaceGraph =
      boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;
  /**
   * @brief Construct from PolyhedralSurface
   * @param surface The polyhedral surface to construct the graph from
   * @throw Exception if surface is not connected
   */
  SurfaceGraph(const PolyhedralSurface &surface);

  /**
   * @brief Construct from TriangulatedSurface
   * @param triangulatedSurface The triangulated surface to construct the graph
   * from
   * @throw Exception if surface is not connected
   */
  SurfaceGraph(const TriangulatedSurface &triangulatedSurface);

  /**
   * @brief Get the edge map
   * @return Const reference to the edge map
   */
  const EdgeMap &
  edgeMap() const
  {
    return _edgeMap;
  }
  /**
   * @brief Get the face adjacency graph
   * @return Const reference to the face graph
   */
  const FaceGraph &
  faceGraph() const
  {
    return _graph;
  }
  // const CoordinateMap & coordMap() const { return _coordinateMap ; }
  /**
   * @brief Check if the surface graph is valid
   * @return Validity status of the surface graph
   */
  const Validity
  isValid() const
  {
    return _isValid;
  }

private:
  CoordinateMap _coordinateMap;
  EdgeMap       _edgeMap;
  FaceGraph     _graph;
  VertexIndex   _numVertices;

  Validity _isValid;

  void
  addRing(const LineString &ring, FaceIndex faceIndex); // helper for ctor
};

/**
 * test if a surface is connected, the graph should be build beforehand
 * @param graph the surface graph to test
 * @return true if the surface is connected, false otherwise
 */
SFCGAL_API bool
isConnected(const SurfaceGraph &graph);

/**
 * test if a surface is closed, the graph should be build beforehand
 * @param graph the surface graph to test
 * @return true if the surface is closed, false otherwise
 * @note the surface may not be connected, eg. two spheres will yield a true
 * result
 */
SFCGAL_API bool
isClosed(const SurfaceGraph &graph);

} // namespace SFCGAL::algorithm
#endif
