// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_GRAPH_VERTEX_H_
#define SFCGAL_GRAPH_VERTEX_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Coordinate.h"

namespace SFCGAL::graph {

/**
 * @brief [private]A vertex in a GeometryGraph with minimal requirements (some
 * algorithms could need a richer class)
 */
struct SFCGAL_API Vertex {
  /**
   * @brief Constructor with coordinate
   * @param coordinate_ The coordinate for this vertex (default empty
   * coordinate)
   */
  Vertex(const Coordinate &coordinate_ = Coordinate());

  Coordinate coordinate; ///< The coordinate of this vertex
};

} // namespace SFCGAL::graph

#endif
