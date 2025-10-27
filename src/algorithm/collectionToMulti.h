// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_COLLECTION_TO_MULTI_ALGORITHM
#define SFCGAL_COLLECTION_TO_MULTI_ALGORITHM

#include "SFCGAL/config.h"

#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryCollection.h"

namespace SFCGAL::algorithm {
/**
 * Given a geometry collection of triangles, TINs and polygons
 * returns a MultiPolygon
 * @param g the geometry collection to convert to multi-geometry
 * @return a MultiPolygon as a unique_ptr<Geometry>
 * @warning Ownership is taken from the parameter
 */
SFCGAL_API std::unique_ptr<Geometry>
           collectionToMulti(std::unique_ptr<Geometry> g);
} // namespace SFCGAL::algorithm

#endif
