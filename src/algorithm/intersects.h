// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_INTERSECTS_ALGORITHM
#define SFCGAL_INTERSECTS_ALGORITHM

#include "SFCGAL/config.h"

namespace SFCGAL {
class Geometry;
class LineString;
class PolyhedralSurface;
class TriangulatedSurface;
namespace detail {
template <int Dim>
class GeometrySet;
template <int Dim>
struct PrimitiveHandle;
} // namespace detail

namespace algorithm {
class SurfaceGraph;
// defined in isValid.h
struct NoValidityCheck;

/**
 * Robust intersection test on 2D geometries. Force projection to z=0 if needed
 * @pre geometry1 and geometry2 are valid geometries
 * @param geometry1 first geometry to operate
 * @param geometry2 second geometry to operate
 * @return true when @p geometry1 intersects @p geometry2
 */
SFCGAL_API auto
intersects(const Geometry &geometry1, const Geometry &geometry2) -> bool;

/**
 * Robust intersection test on 3D geometries. Assume z = 0 if needed
 * @pre geometry1 and geometry2 are valid geometries
 * @param geometry1 first geometry to operate
 * @param geometry2 second geometry to operate
 * @return true when @p geometry1 intersects @p geometry2
 */
SFCGAL_API auto
intersects3D(const Geometry &geometry1, const Geometry &geometry2) -> bool;

/**
 * Intersection test on 2D geometries. Force projection to z=0 if needed
 * @pre geometry1 and geometry2 are valid geometries
 * @param geometry1 first geometry to operate
 * @param geometry2 second geometry to operate
 * @return true when @p geometry1 intersects @p geometry2
 *
 * @warning the validity is assumed, no actual check is done
 */
SFCGAL_API auto
intersects(const Geometry &geometry1, const Geometry &geometry2,
           NoValidityCheck) -> bool;

/**
 * Intersection test on 3D geometries. Assume z = 0 if needed
 * @pre geometry1 and geometry2 are valid geometries
 * @param geometry1 first geometry to operate
 * @param geometry2 second geometry to operate
 * @return true when @p geometry1 intersects @p geometry2
 *
 * @warning the validity is assumed, no actual check is done
 */
SFCGAL_API auto
intersects3D(const Geometry &geometry1, const Geometry &geometry2,
             NoValidityCheck) -> bool;

/**
 * Intersection test on GeometrySet
 * @param geometrySet1 first geometry set to operate
 * @param geometrySet2 second geometry set to operate
 * @return true when @p geometrySet1 intersects @p geometrySet2
 *
 */
template <int Dim>
auto
intersects(const detail::GeometrySet<Dim> &geometrySet1,
           const detail::GeometrySet<Dim> &geometrySet2) -> bool;

/**
 * Intersection test on a PrimitiveHandle
 * @param handle1 first primitive handle to operate
 * @param handle2 second primitive handle to operate
 * @return true when @p handle1 intersects @p handle2
 *
 */
template <int Dim>
auto
intersects(const detail::PrimitiveHandle<Dim> &handle1,
           const detail::PrimitiveHandle<Dim> &handle2) -> bool;

/**
 * Self intersection test for 2D LineString (false if only endpoint touch)
 * @param lineString geometry to check self intersection with
 * @return true when @p lineString self intersects
 *
 */
auto
selfIntersects(const LineString &lineString) -> bool;

/**
 * Self intersection test for 3D LineString (false if only endpoints touch)
 * @param lineString geometry to check self intersection with
 * @return true when @p lineString self intersects
 *
 */
auto
selfIntersects3D(const LineString &lineString) -> bool;

/**
 * Self intersection test for 2D PolyhedralSurface (false if only point touch)
 * @param surface geometry to check self intersection with
 * @param graph surface graph matching @p surface
 * @return true when @p surface self intersects
 *
 */
auto
selfIntersects(const PolyhedralSurface &surface, const SurfaceGraph &graph)
    -> bool;

/**
 * Self intersection test for 3D PolyhedralSurface (false if only point touch)
 * @param surface geometry to check self intersection with
 * @param graph surface graph matching @p surface
 * @return true when @p surface self intersects
 *
 */

auto
selfIntersects3D(const PolyhedralSurface &surface, const SurfaceGraph &graph)
    -> bool;

/**
 * Self intersection test for 2D TriangulatedSurface (false if only point touch)
 * @param triangulatedSurface geometry to check self intersection with
 * @param graph surface graph matching @p triangulatedSurface
 * @return true when @p triangulatedSurface self intersects
 *
 */
auto
selfIntersects(const TriangulatedSurface &triangulatedSurface,
               const SurfaceGraph        &graph) -> bool;

/**
 * Self intersection test for 3D TriangulatedSurface (false if only point touch)
 * @param triangulatedSurface geometry to check self intersection with
 * @param graph surface graph matching @p triangulatedSurface
 * @return true when @p triangulatedSurface self intersects
 *
 */
auto
selfIntersects3D(const TriangulatedSurface &triangulatedSurface,
                 const SurfaceGraph        &graph) -> bool;
} // namespace algorithm
} // namespace SFCGAL
#endif
