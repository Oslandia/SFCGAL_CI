// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_CENTROID_H_
#define SFCGAL_ALGORITHM_CENTROID_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"

namespace SFCGAL {
class Curve;
}

namespace SFCGAL {
namespace algorithm {

/**
 * Holds weighted data used to compute a centroid for a Geometry.
 *
 * Contains accumulated weight (area or length), centroid vector, and
 * weighted average of measured values for centroid computation.
 *
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
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const Geometry &geom, bool enable3DComputation = false);

/**
 * @brief Computes the weighted centroid of a Triangle.
 *
 * @param triangle Input triangle
 * @param enable3DComputation Whether to use 3D computation
 * @return WeightedCentroid result
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const Triangle &triangle, bool enable3DComputation = false);

/**
 * @brief Computes the weighted centroid of a Triangle from 3 Point.
 *
 * @ingroup detail
 */
SFCGAL_API auto
weightedCentroid(const Point &pta, const Point &ptb, const Point &ptc,
                 bool enable3DComputation = false) -> WeightedCentroid;

/**
 * @brief Computes the weighted centroid of a LineString.
 *
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const LineString &lineString,
                 bool              enable3DComputation = false);

/**
 * @brief Computes the weighted centroid for a Curve
 *
 * @ingroup detail
 */
SFCGAL_API auto
weightedCentroid(const Curve &g, bool enable3DComputation = false)
    -> WeightedCentroid;

/**
 * @brief Computes the weighted centroid of a Polygon.
 *
 * @param polygon Input polygon
 * @param enable3DComputation Whether to use 3D computation
 * @return WeightedCentroid result
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const Polygon &polygon, bool enable3DComputation = false);

/**
 * @brief Computes the weighted centroid of a GeometryCollection.
 *
 * @param collection Input collection
 * @param enable3DComputation Whether to use 3D computation
 * @return WeightedCentroid result
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const GeometryCollection &collection,
                 bool                      enable3DComputation = false);

/**
 * @brief Computes the weighted centroid of a TriangulatedSurface.
 *
 * @param tin Input triangulated surface
 * @param enable3DComputation Whether to use 3D computation
 * @return WeightedCentroid result
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const TriangulatedSurface &tin,
                 bool                       enable3DComputation = false);

/**
 * @brief Computes the weighted centroid of a PolyhedralSurface.
 *
 * @param surface Input polyhedral surface
 * @param enable3DComputation Whether to use 3D computation
 * @return WeightedCentroid result
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const PolyhedralSurface &surface,
                 bool                     enable3DComputation = false);

/**
 * @brief Computes the weighted centroid of a Solid.
 *
 * @param solid Input solid
 * @param enable3DComputation Whether to use 3D computation
 * @return WeightedCentroid result
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const Solid &solid, bool enable3DComputation = false);

} // namespace algorithm
} // namespace SFCGAL

#endif
