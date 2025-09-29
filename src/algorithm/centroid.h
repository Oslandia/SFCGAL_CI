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
 * @brief Compute a weighted centroid for a Geometry by dispatching to the
 * appropriate overload for the geometry's concrete type.
 *
 * Dispatches on g.geometryTypeId() and delegates computation to the
 * corresponding weightedCentroid overload (Point, LineString, Curve,
 * Polygon, Triangle, GeometryCollection, TriangulatedSurface,
 * PolyhedralSurface, Solid). For a Point returns a WeightedCentroid with
 * zero area, the point's vector as centroid, and the point's M value if
 * measured.
 *
 * @param geom Geometry to compute the weighted centroid for.
 * @param enable3DComputation When true, area/length computations prefer 3D
 *        formulas (when available); otherwise 2D computations are used.
 * @return WeightedCentroid Aggregated area, centroid vector, and M value
 *         for the input geometry.
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
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const Triangle &triangle, bool enable3DComputation = false);

/**
 * @brief Compute the weighted centroid contribution of a triangle defined by
 * three points.
 *
 * Computes the triangle's area, its centroid vector, and an M-value (measure)
 * used for weighted aggregation across geometries. When enable3DComputation is
 * true the area is computed in 3D (positive); otherwise a 2D signed area is
 * computed from the XY coordinates.
 *
 * The returned WeightedCentroid contains:
 * - area: triangle area (3D: non-negative; 2D: signed area),
 * - centroid: the arithmetic mean of the three point position vectors,
 * - m: the average of the three points' M values if all three are measured,
 * otherwise 0.
 *
 * @param pta First triangle vertex.
 * @param ptb Second triangle vertex.
 * @param ptc Third triangle vertex.
 * @param enable3DComputation If true, compute area in 3D; if false, compute 2D
 * signed area.
 * @return WeightedCentroid Area, centroid vector, and M value for the triangle.
 *
 * @ingroup detail
 */
SFCGAL_API auto
weightedCentroid(const Point &pta, const Point &ptb, const Point &ptc,
                 bool enable3DComputation = false) -> WeightedCentroid;

/**
 * @brief Computes the weighted centroid for a LineString (open polyline or
 * closed ring).
 *
 * Computes an area/length-weighted centroid and the associated weighted measure
 * (M).
 * - If lineString.isClosed() is true, treats the LineString as a polygonal ring
 * and decomposes it into triangles fan-based from the first point; each
 * triangle's area/centroid is accumulated (triangle area computation honors
 * enable3DComputation).
 * - If lineString is open, treats it as a polyline and accumulates segment
 * contributions using segment midpoints weighted by segment length.
 *
 * @param lineString The LineString to process; when closed this represents a
 * polygonal ring.
 * @param enable3DComputation When true, triangle area computations use 3D
 * geometry; otherwise 2D signed area is used for polygonal contributions.
 * @return WeightedCentroid A struct containing the total area/length, the
 * computed centroid (as a Vector_3), and the averaged M value (weighted by
 * area/length).
 * @throws InappropriateGeometryException Thrown when the aggregated total
 * area/length is zero (invalid LineString for centroid computation).
 *
 * @ingroup detail
 */
SFCGAL_API WeightedCentroid
weightedCentroid(const LineString &lineString,
                 bool              enable3DComputation = false);

/**
 * @brief Computes the area/length-weighted centroid of a parametric Curve by
 * approximating it as a LineString.
 *
 * Attempts an adaptive sampling of the curve to produce a LineString
 * approximation; if that fails or yields an empty result it falls back to a
 * denser uniform sampling (256 points). If no valid approximation can be
 * produced, returns an empty WeightedCentroid. Otherwise delegates to the
 * LineString overload to compute the weighted centroid.
 *
 * @param g Curve to approximate and evaluate.
 * @param enable3DComputation If true, use 3D measures (length/area and 3D
 * centroids); otherwise use 2D computations.
 *
 * @note When applied to NURBSCurve geometries, the centroid
 *   is internally computed on a LineString obtained via
 *   toLineString() with its default parameters.
 *
 * @return WeightedCentroid Empty if the curve cannot be approximated as a
 * non-empty LineString, otherwise the computed weighted centroid.
 *
 * @ingroup detail
 */
SFCGAL_API auto
weightedCentroid(const Curve &curve, bool enable3DComputation = false)
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
