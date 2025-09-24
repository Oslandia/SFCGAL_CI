// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_CENTROID_H_
#define SFCGAL_ALGORITHM_CENTROID_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"

namespace SFCGAL {
namespace algorithm {

/**
 * Holds weighted data to compute the centroid for a Geometry
 * @ingroup detail
 */
class SFCGAL_API WeightedCentroid {
public:
  /// total area or length
  SFCGAL::Kernel::FT area;
  /// 3D centroid
  CGAL::Vector_3<SFCGAL::Kernel> centroid;
  /// weighted average of m
  SFCGAL::Kernel::FT m;

  /**
   * @brief Constructs a WeightedCentroid object.
   *
   * Initializes a weighted centroid with the given area, centroid vector, and
   * mass.
   *
   * @param a Initial area (default = 0.0).
   * @param c Initial centroid vector (default = zero vector).
   * @param mTmp Initial mass or weight (default = 0.0).
   */
  WeightedCentroid(
      SFCGAL::Kernel::FT             a    = 0.0,
      CGAL::Vector_3<SFCGAL::Kernel> c    = CGAL::Vector_3<SFCGAL::Kernel>(),
      SFCGAL::Kernel::FT             mTmp = 0.0)
      : area(a), centroid(c), m(mTmp)
  {
  }
};

/**
 * @brief Returns the 2D centroid for a Geometry
 *
 * The result is the weighted centroid of a geometry. The implementation follows
 * PostGIS one (https://postgis.net/docs/ST_Centroid.html). The weigth is
 * computed in the XY space.
 *
 * @warning Z component is ignored, geometries must be valid when projected in
 * the XY plane. Vertical geometries will generate an error.
 * @param geom the input geometry
 * @pre geom is a valid geometry in 2D
 * @return A Point representing the 2D geometry centroid
 */
SFCGAL_API std::unique_ptr<Point>
           centroid(const Geometry &geom);

/**
 * @brief Returns the 3D centroid for a Geometry
 *
 * The result is the weighted centroid of a geometry. The implementation follows
 * PostGIS one (https://postgis.net/docs/ST_Centroid.html). The weigth is
 * computed in the 3D space.
 *
 * @param geom the input geometry
 * @pre geom is a valid geometry
 * @return A Point representing the 3D geometry centroid
 */
SFCGAL_API std::unique_ptr<Point>
           centroid3D(const Geometry &geom);

/**
 * @brief Computes the weighted centroid of a geometry.
 *
 * @param geom The input geometry for which the weighted centroid will be
 * computed.
 * @param enable3DComputation Optional flag (default: false). If set to true,
 *        the function will take the Z-coordinate into account when computing
 *        the centroid for 3D geometries.
 *
 * @return WeightedCentroid The weighted centroid of the input geometry.
 *
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const Geometry &geom, bool enable3DComputation = false);

/**
 * @brief Computes the weighted centroid of a Triangle.
 *
 * @param triangle The input Triangle for which the weighted centroid will be
 * computed.
 * @param enable3DComputation Optional flag (default: false). If set to true,
 *        the function will take the Z-coordinate into account when computing
 *        the centroid for 3D geometries.
 *
 * @return WeightedCentroid The weighted centroid of the input geometry.
 *
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const Triangle &triangle, bool enable3DComputation = false);

/**
 * @brief Computes the weighted centroid of a Triangle from 3 Point.
 *
 * @param pta The first Point of the Triangle
 * @param ptb The second Point of the Triangle
 * @param ptc The third Point of the Triangle
 * @param enable3DComputation Optional flag (default: false). If set to true,
 *        the function will take the Z-coordinate into account when computing
 *        the centroid for 3D geometries.
 *
 * @return WeightedCentroid The weighted centroid of the Triangle.
 *
 * @ingroup detail
 */
SFCGAL_API auto
weightedCentroid(const Point &pta, const Point &ptb, const Point &ptc,
                 bool enable3DComputation = false) -> WeightedCentroid;

/**
 * @brief Computes the weighted centroid of a LineString.
 *
 * @param lineString The input LineString for which the weighted centroid will
 * be computed.
 * @param enable3DComputation Optional flag (default: false). If set to true,
 *        the function will take the Z-coordinate into account when computing
 *        the centroid for 3D geometries.
 *
 * @return WeightedCentroid The weighted centroid of the input geometry.
 *
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const LineString &lineString,
                 bool              enable3DComputation = false);

/**
 * @brief Computes the weighted centroid of a Polygon.
 *
 * @param polygon The input polygon for which the weighted centroid will be
 * computed.
 * @param enable3DComputation Optional flag (default: false). If set to true,
 *        the function will take the Z-coordinate into account when computing
 *        the centroid for 3D geometries.
 *
 * @return WeightedCentroid The weighted centroid of the input geometry.
 *
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const Polygon &polygon, bool enable3DComputation = false);

/**
 * @brief Computes the weighted centroid of a GeometryCollection.
 *
 *
 * @param collection The input collection for which the weighted centroid will
 * be computed.
 * @param enable3DComputation Optional flag (default: false). If set to true,
 *        the function will take the Z-coordinate into account when computing
 *        the centroid for 3D geometries.
 *
 * @return WeightedCentroid The weighted centroid of the input geometry.
 *
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const GeometryCollection &collection,
                 bool                      enable3DComputation = false);

/**
 * @brief Computes the weighted centroid of a TriangulatedSurface.
 *
 * @param tin The input TriangulatedSurface for which the weighted centroid will
 * be computed.
 * @param enable3DComputation Optional flag (default: false). If set to true,
 *        the function will take the Z-coordinate into account when computing
 *        the centroid for 3D geometries.
 *
 * @return WeightedCentroid The weighted centroid of the input geometry.
 *
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const TriangulatedSurface &tin,
                 bool                       enable3DComputation = false);

/**
 * @brief Computes the weighted centroid of a PolyhedralSurface.
 *
 * @param surface The input PolyhedralSurface for which the weighted centroid
 * will be computed.
 * @param enable3DComputation Optional flag (default: false). If set to true,
 *        the function will take the Z-coordinate into account when computing
 *        the centroid for 3D geometries.
 *
 * @return WeightedCentroid The weighted centroid of the input geometry.
 *
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const PolyhedralSurface &surface,
                 bool                     enable3DComputation = false);

/**
 * @brief Computes the weighted centroid of a Solid.
 *
 * @param solid The input Solid for which the weighted centroid will be
 * computed.
 * @param enable3DComputation Optional flag (default: false). If set to true,
 *        the function will take the Z-coordinate into account when computing
 *        the centroid for 3D geometries.
 *
 * @return WeightedCentroid The weighted centroid of the input geometry.
 *
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const Solid &solid, bool enable3DComputation = false);

} // namespace algorithm
} // namespace SFCGAL

#endif
