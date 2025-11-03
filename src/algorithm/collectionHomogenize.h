// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_COLLECTION_HOMOGENIZE_ALGORITHM
#define SFCGAL_COLLECTION_HOMOGENIZE_ALGORITHM

#include "SFCGAL/config.h"

#include <memory>

#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryCollection.h"

namespace SFCGAL::algorithm {
/**
 * Given a geometry collection, returns the "simplest" representation of the
 * contents. Singletons will be returned as singletons. Collections that are
 * homogeneous will be returned as the appropriate multi-type.
 * @param g the geometry collection to homogenize
 * @return the homogenized geometry as a unique_ptr<Geometry>
 * @warning Ownership is taken from the parameter
 */
SFCGAL_API std::unique_ptr<Geometry>
           collectionHomogenize(std::unique_ptr<Geometry> g);
} // namespace SFCGAL::algorithm

#endif
