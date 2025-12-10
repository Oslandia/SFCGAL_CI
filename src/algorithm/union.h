// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_UNION_ALGORITHM
#define SFCGAL_UNION_ALGORITHM

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
 * Union on 2D geometries.
 * @param geometry1 first geometry
 * @param geometry2 second geometry
 * @return union of the two geometries
 * @pre geometry1 and geometry2 are valid geometries
 */
SFCGAL_API auto
union_(const Geometry &geometry1, const Geometry &geometry2)
    -> std::unique_ptr<Geometry>;

/**
 * Union on 2D geometries. No validity check variant
 * @param geometry1 first geometry
 * @param geometry2 second geometry
 * @return union of the two geometries
 * @pre geometry1 and geometry2 are valid geometries
 * @warning No actual validity check is done.
 */
SFCGAL_API auto
union_(const Geometry &geometry1, const Geometry &geometry2, NoValidityCheck)
    -> std::unique_ptr<Geometry>;

/**
 * Union on 3D geometries. Assume z = 0 if needed
 * @param geometry1 first geometry
 * @param geometry2 second geometry
 * @return union of the two geometries
 * @pre geometry1 and geometry2 are valid geometries
 */
SFCGAL_API auto
union3D(const Geometry &geometry1, const Geometry &geometry2)
    -> std::unique_ptr<Geometry>;

/**
 * Union on 3D geometries. Assume z = 0 if needed
 * @param geometry1 first geometry
 * @param geometry2 second geometry
 * @return union of the two geometries
 * @pre geometry1 and geometry2 are valid geometries
 * @warning No actual validity check is done
 */
SFCGAL_API auto
union3D(const Geometry &geometry1, const Geometry &geometry2, NoValidityCheck)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Compute union of two GeometrySet objects
 * @param geometrySet1 First geometry set
 * @param geometrySet2 Second geometry set
 * @param result Output geometry set containing the union
 */
template <int Dim>
auto
union_(const ::SFCGAL::detail::GeometrySet<Dim> &geometrySet1,
       const ::SFCGAL::detail::GeometrySet<Dim> &geometrySet2,
       ::SFCGAL::detail::GeometrySet<Dim>       &result) -> void;

/**
 * @brief Compute union of two PrimitiveHandle objects
 * @param primitiveHandle1 First primitive handle
 * @param primitiveHandle2 Second primitive handle
 * @param result Output geometry set containing the union
 */
template <int Dim>
auto
union_(const ::SFCGAL::detail::PrimitiveHandle<Dim> &primitiveHandle1,
       const ::SFCGAL::detail::PrimitiveHandle<Dim> &primitiveHandle2,
       ::SFCGAL::detail::GeometrySet<Dim>           &result) -> void;

} // namespace algorithm
} // namespace SFCGAL

#endif
