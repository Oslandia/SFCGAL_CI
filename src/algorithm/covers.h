// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_COVERS_ALGORITHM
#define SFCGAL_COVERS_ALGORITHM

#include "SFCGAL/config.h"

#include <vector>

namespace SFCGAL {
class Geometry;
class Solid;
class Point;
namespace detail {
template <int Dim>
class GeometrySet;
template <int Dim>
struct PrimitiveHandle;
} // namespace detail

namespace algorithm {
/**
 * @brief Cover test on 2D geometries. Checks if gA covers gB. Force projection
 * to z=0 if needed
 * @param geometry1 The first geometry
 * @param geometry2 The second geometry
 * @return true if geometry1 covers geometry2, false otherwise
 */
SFCGAL_API auto
covers(const Geometry &geometry1, const Geometry &geometry2) -> bool;

/**
 * @brief Cover test on 3D geometries. Checks if gA covers gB. Assume z = 0 if
 * needed
 * @param geometry1 The first geometry
 * @param geometry2 The second geometry
 * @return true if geometry1 covers geometry2, false otherwise
 */
SFCGAL_API auto
covers3D(const Geometry &geometry1, const Geometry &geometry2) -> bool;

/**
 * @brief Cover test on GeometrySet objects
 * @param geometrySet1 The first geometry set
 * @param geometrySet2 The second geometry set
 * @return true if geometrySet1 covers geometrySet2, false otherwise
 */
template <int Dim>
auto
covers(const detail::GeometrySet<Dim> &geometrySet1,
       const detail::GeometrySet<Dim> &geometrySet2) -> bool;

/**
 * @brief Cover test on PrimitiveHandle objects
 * @param handle1 The first primitive handle
 * @param handle2 The second primitive handle
 * @return true if handle1 covers handle2, false otherwise
 */
template <int Dim>
auto
covers(const detail::PrimitiveHandle<Dim> &handle1,
       const detail::PrimitiveHandle<Dim> &handle2) -> bool;
} // namespace algorithm
} // namespace SFCGAL

#endif
