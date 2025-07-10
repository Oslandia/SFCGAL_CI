// Copyright (c) 2012-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_PARTITION_2_H_
#define SFCGAL_ALGORITHM_PARTITION_2_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Geometry.h"

namespace SFCGAL {
namespace algorithm {
struct NoValidityCheck;

/**
 * Partition algorithm available
 *
 * @since 1.4.2
 */
enum PartitionAlgorithm {
  y_monotone,    /*!< Y Monotone Partition:
                    https://doc.cgal.org/latest/Partition_2/index.html#secpartition_2_monotone
                  */
  approx_convex, /*!< Simple approximation algorithm of Hertel and Mehlhorn
                    https://doc.cgal.org/latest/Partition_2/index.html#secpartition_2_convex
                  */
  greene_approx_convex, /*!< Sweep-line approximation algorithm of Greene
                           https://doc.cgal.org/latest/Partition_2/index.html#secpartition_2_convex
                         */
  optimal_convex        /*!< Optimal convex partition
                           https://doc.cgal.org/latest/Partition_2/index.html#secpartition_2_convex
                         */
};

/**
 * Compute the partition of a 2D polygon
 * https://doc.cgal.org/latest/Partition_2/index.html#Chapter_2D_Polygon_Partitioning
 * @since 1.4.2
 */
SFCGAL_API auto
partition_2(const Geometry &g, PartitionAlgorithm alg = y_monotone)
    -> std::unique_ptr<Geometry>;

/**
 * Compute the partition of a 2D polygon
 * https://doc.cgal.org/latest/Partition_2/index.html#Chapter_2D_Polygon_Partitioning
 * @pre g is a valid geometry
 * @warning No actual validity check is done
 */
SFCGAL_API auto
partition_2(const Geometry &g, PartitionAlgorithm alg, NoValidityCheck)
    -> std::unique_ptr<Geometry>;

} // namespace algorithm
} // namespace SFCGAL

#endif
