// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_GRAPH_EDGE_H_
#define SFCGAL_GRAPH_EDGE_H_

#include "SFCGAL/config.h"

namespace SFCGAL {
namespace graph {

/**
 * @brief [private]An edge in a GeometryGraph with minimal requirements (some
 * algorithms could need more information)
 */
struct SFCGAL_API Edge {
  /**
   * @brief Constructor for Edge
   * @param face_ The face identifier (default -1)
   */
  Edge(const int &face_ = -1);

  int face; ///< Face identifier for this edge
};

} // namespace graph
} // namespace SFCGAL

#endif
