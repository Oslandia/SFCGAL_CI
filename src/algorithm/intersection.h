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
 * @param ga First geometry
 * @param gb Second geometry
 * @pre ga and gb are valid geometries
 * @return The intersection of ga and gb
 */
SFCGAL_API std::unique_ptr<Geometry>
           intersection(const Geometry &ga, const Geometry &gb);

/**
 * Intersection on 2D geometries. No validity check variant
 * @param ga First geometry
 * @param gb Second geometry
 * @param NoValidityCheck Tag type to skip validity checks
 * @pre ga and gb are valid geometries
 * @warning No actual validity check is done.
 * @return The intersection of ga and gb
 */
SFCGAL_API std::unique_ptr<Geometry>
intersection(const Geometry &ga, const Geometry &gb, NoValidityCheck);

/**
 * Intersection on 3D geometries. Assume z = 0 if needed
 * @param ga First geometry
 * @param gb Second geometry
 * @pre ga and gb are valid geometries
 * @return The 3D intersection of ga and gb
 */
SFCGAL_API std::unique_ptr<Geometry>
           intersection3D(const Geometry &ga, const Geometry &gb);

/**
 * Intersection on 3D geometries. Assume z = 0 if needed
 * @param ga First geometry
 * @param gb Second geometry
 * @param NoValidityCheck Tag type to skip validity checks
 * @pre ga and gb are valid geometries
 * @warning No actual validity check is done
 * @return The 3D intersection of ga and gb
 */
SFCGAL_API std::unique_ptr<Geometry>
intersection3D(const Geometry &ga, const Geometry &gb, NoValidityCheck);

/**
 * @ingroup detail
 * @brief Compute intersection between two GeometrySets
 * @tparam Dim Dimension (2 or 3)
 * @param a First geometry set
 * @param b Second geometry set
 * @param output Output geometry set to store results
 */
template <int Dim>
void
intersection(const detail::GeometrySet<Dim> &a,
             const detail::GeometrySet<Dim> &b,
             detail::GeometrySet<Dim>       &output);

/**
 * @ingroup detail
 * @brief Compute intersection between two primitive handles
 * @tparam Dim Dimension (2 or 3)
 * @param a First primitive handle
 * @param b Second primitive handle
 * @param output Output geometry set to store results
 */
template <int Dim>
void
intersection(const detail::PrimitiveHandle<Dim> &a,
             const detail::PrimitiveHandle<Dim> &b,
             detail::GeometrySet<Dim>           &output);

/**
 * @brief Clear cached data used by 3D intersection operations.
 *
 * Should be called periodically to free memory used by internal caches
 * for triangulations and spatial query structures.
 */
SFCGAL_API void
clearIntersectionCaches();

} // namespace algorithm
} // namespace SFCGAL

#endif
