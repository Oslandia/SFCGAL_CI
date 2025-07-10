// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_COLLECTION_TO_MULTI_ALGORITHM
#define SFCGAL_COLLECTION_TO_MULTI_ALGORITHM

#include "SFCGAL/config.h"

#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryCollection.h"

namespace SFCGAL {
namespace algorithm {
/**
 * Given a geometry collection of triangles, TINs and polygons
 * returns a MultiPolygon
 *
 * @warning Ownership is taken from the parameter
 */
SFCGAL_API std::unique_ptr<Geometry>
           collectionToMulti(std::unique_ptr<Geometry> coll);
} // namespace algorithm
} // namespace SFCGAL

#endif
