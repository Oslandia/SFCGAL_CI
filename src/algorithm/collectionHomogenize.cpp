// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/collectionHomogenize.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"

namespace SFCGAL::algorithm {

// If nothing has to be built, g will be moved to the result without
// copying and a new allocation. Otherwise, a new geometry is built and
// the old one is deleted.
auto
collectionHomogenize(std::unique_ptr<Geometry> g) -> std::unique_ptr<Geometry>
{
  // unknown type
  int common_type = 0;

  if (!g->is<GeometryCollection>()) {
    // not a collection, nothing to do
    return g;
  }

  const GeometryCollection &coll = g->as<GeometryCollection>();

  // test if it is a singleton
  if (coll.numGeometries() == 1) {
    return std::unique_ptr<Geometry>(coll.geometryN(0).clone());
  }

  for (size_t i = 0; i < coll.numGeometries(); ++i) {
    const Geometry &gi = coll.geometryN(i);

    if (common_type == 0) {
      common_type = gi.geometryTypeId();
      continue;
    }

    if (gi.geometryTypeId() != common_type) {
      common_type = 0;
      break;
    }
  }

  if (common_type == 0) {
    // not an homogeneous collection, give back
    return g;
  }

  GeometryCollection *ret_geo = nullptr;

  if (common_type == TYPE_POINT) {
    ret_geo = new MultiPoint;
  } else if (common_type == TYPE_LINESTRING) {
    ret_geo = new MultiLineString;
  } else if (common_type == TYPE_POLYGON) {
    ret_geo = new MultiPolygon;
  } else if (common_type == TYPE_SOLID) {
    ret_geo = new MultiSolid;
  }

  // copy each geometry
  for (size_t i = 0; i < coll.numGeometries(); ++i) {
    ret_geo->addGeometry(coll.geometryN(i));
  }

  return std::unique_ptr<Geometry>(ret_geo);
}

} // namespace SFCGAL::algorithm
