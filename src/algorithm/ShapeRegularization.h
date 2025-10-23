// Copyright (c) 2025-2025, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_SHAPEREGULARIZATION_H
#define SFCGAL_ALGORITHM_SHAPEREGULARIZATION_H

#include "SFCGAL/config.h"
#include "SFCGAL/Geometry.h"
#include <vector>

namespace SFCGAL {
namespace algorithm {

struct NoValidityCheck;

/**
 * @brief Shape regularization API based on CGAL Shape Regularization package
 *
 * Provides algorithms to regularize geometric shapes by aligning segments and contours
 * to principal directions, improving orthogonality and removing small deviations.
 *
 * Supported operations:
 * - **Segment regularization**: align segments to be parallel, orthogonal, or collinear
 * - **Open contour regularization**: regularize polylines
 * - **Closed contour regularization**: regularize polygons while preserving topology
 *
 * The implementation preserves Z/M coordinates through interpolation when possible.
 * All regularization operations work in 2D (XY plane).
 *
 * @note Closed contour regularization uses segment-based regularization to avoid
 *       known stability issues in CGAL's contour regularization algorithms.
 */
class SFCGAL_API ShapeRegularization {
public:
  /**
   * @brief Neighbor query strategies for segment regularization
   */
  enum class NeighborQuery {
    DELAUNAY,       ///< Delaunay triangulation-based (default, most robust)
    NEAREST,        ///< K-nearest neighbors (not yet implemented)
    HIERARCHICAL    ///< Hierarchical clustering (not yet implemented)
  };

  /**
   * @brief Direction estimators for contour regularization
   */
  enum class DirectionEstimator {
    LONGEST,        ///< Use longest edge direction as principal direction
    MULTIPLE        ///< Detect multiple principal directions
  };

  // ============================================================================
  // SEGMENT REGULARIZATION
  // ============================================================================

  /**
   * @brief Regularize 2D segments extracted from geometry
   *
   * Extracts all segments from the input geometry and applies regularization to:
   * - Align segments to be parallel or orthogonal (angle regularization)
   * - Align parallel segments to be collinear (offset regularization)
   *
   * @param geometry Source geometry (LineString, MultiLineString, Polygon, 
   *                 MultiPolygon, Triangle, or GeometryCollection)
   * @param applyAngles Enable angle regularization (parallelism/orthogonality)
   * @param applyOffsets Enable offset regularization (collinearity)
   * @param angleBoundDeg Maximum angle deviation in degrees (0-90)
   * @param offsetBound Maximum offset distance for collinearity
   * @param angleTolerance Tolerance for considering angles equal (degrees)
   * @param query Neighbor query strategy (only DELAUNAY is currently supported)
   * @return MultiLineString containing regularized 2-point segments
   * @throws std::invalid_argument if geometry is not 2D
   */
  static auto
  regularizeSegments(const Geometry &geometry, 
                     bool applyAngles = true,
                     bool applyOffsets = false, 
                     double angleBoundDeg = 25.0,
                     double offsetBound = 0.5,
                     double angleTolerance = 5.0,
                     NeighborQuery query = NeighborQuery::DELAUNAY) 
      -> std::unique_ptr<Geometry>;

  /**
   * @brief Regularize segments without validity check
   * @copydoc regularizeSegments
   * @warning No validity check performed
   */
  static auto
  regularizeSegments(const Geometry &geometry, bool applyAngles,
                     bool applyOffsets, double angleBoundDeg, double offsetBound,
                     double angleTolerance, NeighborQuery query,
                     NoValidityCheck) -> std::unique_ptr<Geometry>;


  // ============================================================================
  // CONTOUR REGULARIZATION
  // ============================================================================

  /**
   * @brief Regularize a closed contour (polygon or closed linestring)
   *
   * Regularizes closed contours by aligning segments to principal directions.
   * Uses segment-based regularization to avoid stability issues in CGAL.
   *
   * @param geometry Source geometry (Polygon with optional holes, closed LineString,
   *                 or MultiPolygon)
   * @param estimator Direction estimation strategy (currently not used,
   *                  segment-based approach is applied instead)
   * @param angleToleranceDeg Angle tolerance for direction detection (degrees)
   * @param minSegmentLength Minimum segment length to consider
   * @return Regularized geometry of the same type as input
   * @throws std::invalid_argument if geometry is not 2D or not closed
   * @note The estimator parameter is provided for API compatibility but
   *       the implementation uses segment-based regularization.
   */
  static auto
  regularizeClosedContour(const Geometry &geometry,
                          DirectionEstimator estimator = DirectionEstimator::MULTIPLE,
                          double angleToleranceDeg = 10.0,
                          double minSegmentLength = 0.01)
      -> std::unique_ptr<Geometry>;

  /**
   * @brief Regularize closed contour without validity check
   * @copydoc regularizeClosedContour
   * @warning No validity check performed
   */
  static auto
  regularizeClosedContour(const Geometry &geometry,
                          DirectionEstimator estimator,
                          double angleToleranceDeg,
                          double minSegmentLength,
                          NoValidityCheck)
      -> std::unique_ptr<Geometry>;

  /**
   * @brief Regularize an open contour (polyline)
   *
   * Regularizes open contours by aligning them to principal directions.
   * Detects dominant directions and adjusts vertices accordingly.
   *
   * @param geometry Source geometry (open LineString or MultiLineString)
   * @param estimator Direction estimation strategy (LONGEST or MULTIPLE)
   * @param angleToleranceDeg Angle tolerance for direction detection (degrees)
   * @param minSegmentLength Minimum segment length to consider
   * @return Regularized geometry of the same type as input
   * @throws std::invalid_argument if geometry is not 2D
   */
  static auto
  regularizeOpenContour(const Geometry &geometry,
                        DirectionEstimator estimator = DirectionEstimator::LONGEST,
                        double angleToleranceDeg = 10.0,
                        double minSegmentLength = 0.01)
      -> std::unique_ptr<Geometry>;

  /**
   * @brief Regularize open contour without validity check
   * @copydoc regularizeOpenContour
   * @warning No validity check performed
   */
  static auto
  regularizeOpenContour(const Geometry &geometry, 
                        DirectionEstimator estimator,
                        double angleToleranceDeg,
                        double minSegmentLength,
                        NoValidityCheck) -> std::unique_ptr<Geometry>;

};

} // namespace algorithm
} // namespace SFCGAL

#endif // SFCGAL_ALGORITHM_SHAPEREGULARIZATION_H
