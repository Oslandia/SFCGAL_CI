// Copyright (c) 2025-2025, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_SURFACE_SIMPLIFICATION_H
#define SFCGAL_ALGORITHM_SURFACE_SIMPLIFICATION_H

#include "SFCGAL/Geometry.h"
#include "SFCGAL/config.h"

#include <memory>

namespace SFCGAL::algorithm {
struct NoValidityCheck;

/**
 * @brief Strategy for surface mesh simplification cost and placement
 *
 * Defines the method used to calculate edge collapse cost and vertex placement
 * during surface mesh simplification. Only strategies compatible with exact
 * kernel are supported.
 */
enum class SimplificationStrategy {
  /**
   * @brief Edge Length strategy
   *
   * Uses edge length as cost function and midpoint placement for vertex
   * positioning. This strategy is compatible with exact kernels and provides
   * good simplification results while maintaining geometric accuracy.
   *
   * @see https://doc.cgal.org/latest/Surface_mesh_simplification/index.html
   */
  EDGE_LENGTH
};

/**
 * @brief Stop predicate for surface mesh simplification
 *
 * Defines when the simplification process should stop.
 */
struct SimplificationStopPredicate {
  /**
   * @brief Type of stop predicate
   */
  enum class Type {
    /**
     * @brief Stop when edge count reaches the specified value
     */
    EDGE_COUNT,

    /**
     * @brief Stop when edge count ratio reaches the specified value
     *
     * The ratio is calculated as: remaining_edges / original_edges
     */
    EDGE_COUNT_RATIO
  };

  Type   type;  ///< Type of stop predicate
  double value; ///< Target value (count or ratio)

  /**
   * @brief Create an edge count stop predicate
   * @param count Target number of edges to keep
   * @return Stop predicate
   */
  static auto
  edgeCount(size_t count) -> SimplificationStopPredicate
  {
    return {Type::EDGE_COUNT, static_cast<double>(count)};
  }

  /**
   * @brief Create an edge count ratio stop predicate
   * @param ratio Target ratio of edges to keep (0.0 to 1.0)
   * @return Stop predicate
   */
  static auto
  edgeCountRatio(double ratio) -> SimplificationStopPredicate
  {
    return {Type::EDGE_COUNT_RATIO, ratio};
  }
};

/**
 * @brief Simplify a surface mesh using CGAL edge collapse algorithm
 *
 * Simplifies a 3D surface mesh by iteratively collapsing edges while
 * minimizing geometric error. Supports TriangulatedSurface, PolyhedralSurface,
 * Solid, and MultiSolid geometries.
 *
 * The algorithm uses CGAL's Surface_mesh_simplification package to perform
 * edge collapse operations. The cost and placement of edge collapses are
 * determined by the specified strategy.
 *
 * @param geometry The input geometry to simplify (must be a surface or solid)
 * @param stopPredicate When to stop the simplification process
 * @param strategy The cost and placement strategy to use (default: EDGE_LENGTH)
 *
 * @return A simplified copy of the input geometry
 *
 * @pre The input geometry must be valid and non-empty
 * @pre For EDGE_COUNT_RATIO, the ratio must be in the range (0.0, 1.0)
 * @pre The geometry must be 3-dimensional
 *
 * @warning Surface geometries (polygon soups) will be converted to
 * CGAL::Surface_mesh, which may modify the topology if the input is not
 * manifold
 *
 * @see
 * https://doc.cgal.org/latest/Surface_mesh_simplification/index.html#Chapter_Triangulated_Surface_Mesh_Simplification
 *
 * @example
 * @code
 * auto simplified = surfaceSimplification(
 *     polyhedralSurface,
 *     SimplificationStopPredicate::edgeCountRatio(0.5),
 *     SimplificationStrategy::EDGE_LENGTH
 * );
 * @endcode
 */
SFCGAL_API auto
surfaceSimplification(const Geometry                 &geometry,
                      const SimplificationStopPredicate &stopPredicate,
                      SimplificationStrategy strategy =
                          SimplificationStrategy::EDGE_LENGTH)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Simplify a surface mesh using CGAL edge collapse algorithm (no
 * validity check)
 *
 * @param geometry The input geometry to simplify
 * @param stopPredicate When to stop the simplification process
 * @param strategy The cost and placement strategy to use
 * @param noCheck Validity check bypass tag
 *
 * @return A simplified copy of the input geometry
 *
 * @warning No validity check is performed on the input geometry
 *
 * @see surfaceSimplification(const Geometry&, const
 * SimplificationStopPredicate&, SimplificationStrategy)
 */
SFCGAL_API auto
surfaceSimplification(const Geometry                 &geometry,
                      const SimplificationStopPredicate &stopPredicate,
                      SimplificationStrategy             strategy,
                      NoValidityCheck                    noCheck)
    -> std::unique_ptr<Geometry>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_SURFACE_SIMPLIFICATION_H
