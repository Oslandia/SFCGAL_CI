// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_INTERSECTION_ALGORITHM
#define SFCGAL_INTERSECTION_ALGORITHM

#include "SFCGAL/config.h"

#include <memory>

namespace SFCGAL {
class Geometry;
namespace detail {
template <int Dim>
class GeometrySet;
template <int Dim>
struct PrimitiveHandle;
} // namespace detail

namespace algorithm {
struct NoValidityCheck;

/**
 * Intersection on 2D geometries.
 * @param geometry1 the first geometry
 * @param geometry2 the second geometry
 * @return the intersection of the two geometries as a unique_ptr<Geometry>
 * @pre geometry1 and geometry2 are valid geometries
 */
SFCGAL_API auto
intersection(const Geometry &geometry1, const Geometry &geometry2)
    -> std::unique_ptr<Geometry>;

/**
 * Intersection on 2D geometries. No validity check variant
 * @param geometry1 the first geometry
 * @param geometry2 the second geometry
 * @return the intersection of the two geometries as a unique_ptr<Geometry>
 * @pre geometry1 and geometry2 are valid geometries
 * @warning No actual validity check is done.
 */
SFCGAL_API auto
intersection(const Geometry &geometry1, const Geometry &geometry2,
             NoValidityCheck) -> std::unique_ptr<Geometry>;

/**
 * Intersection on 3D geometries. Assume z = 0 if needed
 * @param geometry1 the first geometry
 * @param geometry2 the second geometry
 * @return the 3D intersection of the two geometries as a unique_ptr<Geometry>
 * @pre geometry1 and geometry2 are valid geometries
 */
SFCGAL_API auto
intersection3D(const Geometry &geometry1, const Geometry &geometry2)
    -> std::unique_ptr<Geometry>;

/**
 * Intersection on 3D geometries. Assume z = 0 if needed
 * @param geometry1 the first geometry
 * @param geometry2 the second geometry
 * @return the 3D intersection of the two geometries as a unique_ptr<Geometry>
 * @pre geometry1 and geometry2 are valid geometries
 * @warning No actual validity check is done
 */
SFCGAL_API auto
intersection3D(const Geometry &geometry1, const Geometry &geometry2,
               NoValidityCheck) -> std::unique_ptr<Geometry>;

/**
 * Intersection between two geometry sets
 * @param geometrySet1 the first geometry set
 * @param geometrySet2 the second geometry set
 * @param result the result geometry set (output parameter)
 */
template <int Dim>
auto
intersection(const ::SFCGAL::detail::GeometrySet<Dim> &geometrySet1,
             const ::SFCGAL::detail::GeometrySet<Dim> &geometrySet2,
             ::SFCGAL::detail::GeometrySet<Dim>       &result) -> void;

/**
 * Intersection between two primitive handles
 * @param primitiveHandle1 the first primitive handle
 * @param primitiveHandle2 the second primitive handle
 * @param result the result geometry set (output parameter)
 */
template <int Dim>
auto
intersection(const ::SFCGAL::detail::PrimitiveHandle<Dim> &primitiveHandle1,
             const ::SFCGAL::detail::PrimitiveHandle<Dim> &primitiveHandle2,
             ::SFCGAL::detail::GeometrySet<Dim>           &result) -> void;

} // namespace algorithm
} // namespace SFCGAL

#endif
