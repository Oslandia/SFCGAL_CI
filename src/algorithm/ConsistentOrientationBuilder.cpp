// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/ConsistentOrientationBuilder.h"
#include "SFCGAL/detail/graph/algorithm/orientation.h"

namespace SFCGAL::algorithm {

///
///
///
ConsistentOrientationBuilder::ConsistentOrientationBuilder()
    : _graph(), _graphBuilder(_graph)
{
}

///
///
///
void
ConsistentOrientationBuilder::addTriangle(const Triangle &triangle)
{
  _triangles.push_back(
      _graphBuilder.addTriangle(triangle, graph::Edge(_triangles.size())));
}

///
///
///
void
ConsistentOrientationBuilder::addTriangulatedSurface(
    const TriangulatedSurface &triangulatedSurface)
{
  for (size_t i = 0; i < triangulatedSurface.numGeometries(); i++) {
    addTriangle(triangulatedSurface.geometryN(i));
  }
}

///
///
///
auto
ConsistentOrientationBuilder::buildTriangulatedSurface() -> TriangulatedSurface
{
  _makeOrientationConsistent();
  TriangulatedSurface triangulatedSurface;

  for (size_t i = 0; i < numTriangles(); i++) {
    triangulatedSurface.addTriangle(triangleN(i));
  }

  return triangulatedSurface;
}

///
///
///
auto
ConsistentOrientationBuilder::triangleN(const size_t &n) const -> Triangle
{
  const edge_descriptor &ab = _triangles[n][0];
  const edge_descriptor &bc = _triangles[n][1];
  const edge_descriptor &ca = _triangles[n][2];

  return Triangle(Point(_graph[_graph.source(ab)].coordinate),
                  Point(_graph[_graph.source(bc)].coordinate),
                  Point(_graph[_graph.source(ca)].coordinate));
}

///
///
///
void
ConsistentOrientationBuilder::_makeOrientationConsistent()
{
  if (_triangles.empty()) {
    return;
  }

  /*
   * mark all triangles as not oriented and not visited
   */
  _visited.resize(numTriangles());
  _oriented.resize(numTriangles());

  for (size_t i = 0; i < numTriangles(); i++) {
    _visited[i]  = false;
    _oriented[i] = false;
  }

  _computeNeighbors();

  // mark first one as oriented (reference)
  int currentTriangle = -1;

  while ((currentTriangle = _findNextTriangle()) != -1) {
    // mark triangle as visited
    _visited[currentTriangle] = true;

    // orient neighbors
    const std::set<size_t> &neighbors = _neighbors[currentTriangle];

    for (unsigned long const neighbor : neighbors) {
      bool hasOppositeEdge = false;
      bool hasParallelEdge = false;
      graph::algorithm::studyOrientation(_graph, _triangles[currentTriangle],
                                         _triangles[neighbor], hasOppositeEdge,
                                         hasParallelEdge);

      // orientation is consistent
      if (!hasParallelEdge) {
        _oriented[neighbor] = true;
        continue;
      }

      // orientation can't be consistent
      if (hasOppositeEdge && hasParallelEdge) {
        BOOST_THROW_EXCEPTION(
            Exception("can't build consistent orientation from triangle set"));
      }

      // orientation has already been fixed (moebius)
      if (hasParallelEdge && _oriented[neighbor]) {
        BOOST_THROW_EXCEPTION(
            Exception("can't build consistent orientation from triangle set, "
                      "inconsistent orientation for triangle"));
      }

      // here, neighbor triangle should be reversed
      _graph.reverse(_triangles[neighbor]);
      _oriented[neighbor] = true;
    }
  }
}

///
///
///
void
ConsistentOrientationBuilder::_computeNeighbors()
{
  _neighbors.clear();
  _neighbors.resize(numTriangles());

  for (size_t i = 0; i < _triangles.size(); i++) {
    const std::vector<edge_descriptor> &triangle = _triangles[i];

    for (const auto &j : triangle) {
      vertex_descriptor source = _graph.source(j);
      vertex_descriptor target = _graph.target(j);

      // get neighbor edges
      std::vector<directed_edge_descriptor> const neighborEdges =
          _graph.edges(source, target);

      // use marker to fill neighborGraph
      for (const auto &neighborEdge : neighborEdges) {
        auto idOtherTriangle = (size_t)_graph[neighborEdge.first].face;

        if (idOtherTriangle == i) {
          continue;
        }

        _neighbors[i].insert(idOtherTriangle);
      }
    }
  }
}

///
///
///
auto
ConsistentOrientationBuilder::_findNextTriangle() -> int
{
  int result = -1;

  /*
   * find an oriented triangle (reached) and not visited
   */
  for (size_t i = 0; i < numTriangles(); i++) {
    if (!_oriented[i] || _visited[i]) {
      continue;
    }

    result = i;
    break;
  }

  // triangle found
  if (result != -1) {
    return result;
  }

  /*
   * here, a new connected part begins
   */
  for (size_t i = 0; i < numTriangles(); i++) {
    if (!_oriented[i]) {
      _oriented[i] = true;
      return i;
    }
  }

  BOOST_ASSERT(result == -1);
  return result;
}

} // namespace SFCGAL::algorithm
