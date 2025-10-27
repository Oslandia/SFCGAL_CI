// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_DIFFERENCE_ALGORITHM
#define SFCGAL_DIFFERENCE_ALGORITHM

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
 * Difference on 2D geometries.
 * @param geometry1 the first geometry
 * @param geometry2 the second geometry
 * @return the difference of the two geometries as a unique_ptr<Geometry>
 * @pre geometry1 and geometry2 are valid geometries
 */
SFCGAL_API std::unique_ptr<Geometry>
           difference(const Geometry &geometry1, const Geometry &geometry2);

/**
 * Difference on 2D geometries. No validity check variant
 * @param geometry1 the first geometry
 * @param geometry2 the second geometry
 * @return the difference of the two geometries as a unique_ptr<Geometry>
 * @pre geometry1 and geometry2 are valid geometries
 * @warning No actual validity check is done.
 */
SFCGAL_API std::unique_ptr<Geometry>
           difference(const Geometry &geometry1, const Geometry &geometry2,
                      NoValidityCheck);

/**
 * Difference on 3D geometries. Assume z = 0 if needed
 * @param geometry1 the first geometry
 * @param geometry2 the second geometry
 * @return the 3D difference of the two geometries as a unique_ptr<Geometry>
 * @pre geometry1 and geometry2 are valid geometries
 */
SFCGAL_API std::unique_ptr<Geometry>
           difference3D(const Geometry &geometry1, const Geometry &geometry2);

/**
 * Difference on 3D geometries. Assume z = 0 if needed
 * @param geometry1 the first geometry
 * @param geometry2 the second geometry
 * @return the 3D difference of the two geometries as a unique_ptr<Geometry>
 * @pre geometry1 and geometry2 are valid geometries
 * @warning No actual validity check is done
 */
SFCGAL_API std::unique_ptr<Geometry>
           difference3D(const Geometry &geometry1, const Geometry &geometry2,
                        NoValidityCheck);

/**
 * Difference between two geometry sets
 * @param geometrySet1 the first geometry set
 * @param geometrySet2 the second geometry set
 * @param result the result geometry set (output parameter)
 */
template <int Dim>
void
difference(const detail::GeometrySet<Dim> &geometrySet1,
           const detail::GeometrySet<Dim> &geometrySet2,
           detail::GeometrySet<Dim>       &result);

} // namespace algorithm
} // namespace SFCGAL

#endif
