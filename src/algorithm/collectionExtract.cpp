// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/collectionExtract.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"

#include "SFCGAL/detail/transform/ForceZ.h"

namespace SFCGAL::algorithm {

// If nothing has to be built, g will be moved to the result without
// copying and a new allocation. Otherwise, a new geometry is built and
// the old one is deleted.
/// @private
auto
collectionExtractPolygons(std::unique_ptr<Geometry> g)
    -> std::unique_ptr<Geometry>
{
  if (!g->is<GeometryCollection>()) {
    // not a collection, nothing to do
    return g;
  }

  const GeometryCollection &coll = g->as<GeometryCollection>();

  // if it is empty, do not do anything
  if (coll.isEmpty()) {
    return g;
  }

  auto ret_geo = std::make_unique<MultiPolygon>();

  // copy each geometry
  for (size_t i = 0; i < coll.numGeometries(); ++i) {

    std::unique_ptr<Geometry> gi = coll.geometryN(i).clone();

    switch (gi->geometryTypeId()) {
    case TYPE_POLYGON:
    case TYPE_TRIANGLE:
      ret_geo->addGeometry(std::move(gi));
      break;

    case TYPE_POLYHEDRALSURFACE:
    case TYPE_TRIANGULATEDSURFACE: {
      for (size_t j = 0; j < gi->numGeometries(); ++j) {
        ret_geo->addGeometry(gi->geometryN(j));
      }
    } break;

    default:
      // nothing
      break;
    }
  }

  return ret_geo;
}

} // namespace SFCGAL::algorithm
