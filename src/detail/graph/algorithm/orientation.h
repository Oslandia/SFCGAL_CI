// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_GRAPH_ALGORITHM_MAKECONSISTENTORIENTATION_H_
#define SFCGAL_GRAPH_ALGORITHM_MAKECONSISTENTORIENTATION_H_

#include "SFCGAL/detail/graph/GeometryGraph.h"
#include <map>

namespace SFCGAL {
namespace graph {
namespace algorithm {

/**
 * @brief [private]Study orientation between two EdgeStrings
 * @tparam Graph The graph type
 * @param graph The geometry graph
 * @param reference Reference edge vector
 * @param target Target edge vector
 * @param hasOppositeEdge Output flag for opposite edge detection
 * @param hasParallelEdge Output flag for parallel edge detection
 */
template <typename Graph>
void
studyOrientation(Graph                                        &graph,
                 std::vector<typename Graph::edge_descriptor> &reference,
                 std::vector<typename Graph::edge_descriptor> &target,
                 bool &hasOppositeEdge, bool &hasParallelEdge)
{
  /*
   * look for opposite or parallel edges in reference and target
   */
  hasOppositeEdge = false;
  hasParallelEdge = false;

  for (size_t i = 0; i < reference.size(); i++) {
    for (size_t j = 0; j < target.size(); j++) {
      if (graph.areOpposite(reference[i], target[j])) {
        hasOppositeEdge = true;
      }

      if (graph.areParallel(reference[i], target[j])) {
        hasParallelEdge = true;
      }
    }
  }
}

/**
 * @brief Try to build consistent orientation between two edge string
 * @tparam Graph The graph type
 * @param graph The geometry graph
 * @param reference Reference edge vector
 * @param target Target edge vector
 * @return true on success
 */
template <typename Graph>
auto
makeConsistentOrientation(
    Graph &graph, std::vector<typename Graph::edge_descriptor> &reference,
    std::vector<typename Graph::edge_descriptor> &target) -> bool
{
  /*
   * look for opposite or parallel edges in "reference" and "target" edge sets
   */
  bool hasOppositeEdge, hasParallelEdge;
  studyOrientation(graph, reference, target, hasOppositeEdge, hasParallelEdge);

  /*
   * if both opposite and parallel edge are found, there is no possible
   * consistent orientation
   */
  if (hasOppositeEdge && hasParallelEdge) {
    BOOST_THROW_EXCEPTION(
        Exception("can't make consistent orientation between EdgeStrings (both "
                  "opposite and parallel edge found)"));
  }

  /*
   * Only parallel edge found, lets revert the orientation of the "target" edge
   * set
   */
  if (hasParallelEdge) {
    graph.reverse(target);
    return true;
  } else {
    return false;
  }
}

} // namespace algorithm
} // namespace graph
} // namespace SFCGAL

#endif
