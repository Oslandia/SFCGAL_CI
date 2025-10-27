// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_PARTITION_2_H_
#define SFCGAL_ALGORITHM_PARTITION_2_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Geometry.h"

namespace SFCGAL::algorithm {
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
 * @param geometry input geometry
 * @param algorithm partition algorithm to use
 * @return partitioned geometry
 * @since 1.4.2
 */
SFCGAL_API auto
partition_2(const Geometry &geometry, PartitionAlgorithm algorithm = y_monotone)
    -> std::unique_ptr<Geometry>;

/**
 * Compute the partition of a 2D polygon
 * https://doc.cgal.org/latest/Partition_2/index.html#Chapter_2D_Polygon_Partitioning
 * @param geometry input geometry
 * @param algorithm partition algorithm to use
 * @return partitioned geometry
 * @pre geometry is a valid geometry
 * @warning No actual validity check is done
 */
SFCGAL_API auto
partition_2(const Geometry &geometry, PartitionAlgorithm algorithm,
            NoValidityCheck) -> std::unique_ptr<Geometry>;

} // namespace SFCGAL::algorithm

#endif
