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
   * Default constructor
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
 * @pre g is a valid geometry in 2D
 */
SFCGAL_API std::unique_ptr<Point>
           centroid(const Geometry &g);

/**
 * @brief Returns the 3D centroid for a Geometry
 *
 * The result is the weighted centroid of a geometry. The implementation follows
 * PostGIS one (https://postgis.net/docs/ST_Centroid.html). The weigth is
 * computed in the 3D space.
 *
 * @pre g is a valid geometry
 */
SFCGAL_API std::unique_ptr<Point>
           centroid3D(const Geometry &g);

/**
 * Returns the weighted centroid for a Geometry
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const Geometry &g, bool enable3DComputation = false);

/**
 * Returns the weighted centroid for a Triangle
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const Triangle &g, bool enable3DComputation = false);

/**
 * Returns the weighted centroid for a Triangle
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const Point &a, const Point &b, const Point &c,
                 bool enable3DComputation = false);

/**
 * Returns the weighted centroid for a LineString
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const LineString &g, bool enable3DComputation = false);

/**
 * Returns the weighted centroid for a Polygon
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const Polygon &g, bool enable3DComputation = false);

/**
 * Returns the weighted centroid for a GeometryCollection
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const GeometryCollection &g, bool enable3DComputation = false);

/**
 * Returns the weighted centroid for a TriangulatedSurface
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const TriangulatedSurface &g,
                 bool                       enable3DComputation = false);

/**
 * Returns the weighted centroid for a PolyhedralSurface
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const PolyhedralSurface &g, bool enable3DComputation = false);

/**
 * Returns the weighted centroid for a Solid
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const Solid &g, bool enable3DComputation = false);

} // namespace algorithm
} // namespace SFCGAL

#endif
