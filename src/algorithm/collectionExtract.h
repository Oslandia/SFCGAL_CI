// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_COLLECTION_EXTRACT_ALGORITHM
#define SFCGAL_COLLECTION_EXTRACT_ALGORITHM

#include "SFCGAL/config.h"

#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryCollection.h"

namespace SFCGAL::algorithm {
/**
 * Given a geometry collection
 * returns a MultiPolygon from triangles, polygons, polyhedral and polygons
 * @param g the geometry collection to extract polygons from
 * @return a MultiPolygon containing extracted polygons as a
 * unique_ptr<Geometry>
 * @warning Ownership is taken from the parameter
 */
SFCGAL_API std::unique_ptr<Geometry>
           collectionExtractPolygons(std::unique_ptr<Geometry> g);
} // namespace SFCGAL::algorithm

#endif
