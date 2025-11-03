// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_CONSISTENTORIENTATIONBUILDER_H_
#define SFCGAL_ALGORITHM_CONSISTENTORIENTATIONBUILDER_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Geometry.h"

#include "SFCGAL/detail/graph/GeometryGraph.h"
#include "SFCGAL/detail/graph/GeometryGraphBuilder.h"

namespace SFCGAL::algorithm {

/**
 * Make orientation consistent in a triangle set
 */
class SFCGAL_API ConsistentOrientationBuilder {
public:
  /// Vertex descriptor type for geometry graph
  using vertex_descriptor = graph::GeometryGraph::vertex_descriptor;
  /// Edge descriptor type for geometry graph
  using edge_descriptor = graph::GeometryGraph::edge_descriptor;
  /// Directed edge descriptor type for geometry graph
  using directed_edge_descriptor =
      graph::GeometryGraph::directed_edge_descriptor;

  /**
   * default constructor
   */
  ConsistentOrientationBuilder();

  /**
   * @brief Add a Triangle to the builder
   * @param triangle The triangle to add
   */
  void
  addTriangle(const Triangle &triangle);
  /**
   * @brief Add a TriangulatedSurface to the builder
   * @param triangulatedSurface The triangulated surface to add
   */
  void
  addTriangulatedSurface(const TriangulatedSurface &triangulatedSurface);

  /**
   * @brief Get the resulting TriangulatedSurface where each connected part
   * has consistent orientation.
   *
   * @return TriangulatedSurface with consistent orientation
   * @throw SFCGAL::Exception if such a TriangulatedSurface can't be built
   */
  TriangulatedSurface
  buildTriangulatedSurface();

  /**
   * @brief Returns the number of triangles
   * @return The number of triangles in the builder
   */
  inline size_t
  numTriangles() const
  {
    return _triangles.size();
  }
  /**
   * @brief Returns the n-th triangle
   * @param n Index of the triangle to retrieve
   * @return The triangle at index n
   */
  Triangle
  triangleN(const size_t &n) const;

  /**
   * @brief Get the neighbors of the n-th triangle (advanced usage)
   * @param n Index of the triangle
   * @return Set of indices of neighboring triangles
   * @note Use after buildTriangulatedSurface
   */
  const std::set<size_t> &
  neighbors(const size_t &n) const;

private:
  graph::GeometryGraph                      _graph;
  graph::GeometryGraphBuilder               _graphBuilder;
  std::vector<std::vector<edge_descriptor>> _triangles;

  std::vector<bool>             _visited;
  std::vector<bool>             _oriented;
  std::vector<std::set<size_t>> _neighbors;

  /**
   * make triangle orientation consistent
   */
  void
  _makeOrientationConsistent();

  /**
   * compute neighbors for each triangles
   */
  void
  _computeNeighbors();

  /**
   * find the next triangle to visit (may select a new reference triangle)
   */
  int
  _findNextTriangle();
};

} // namespace SFCGAL::algorithm

#endif
