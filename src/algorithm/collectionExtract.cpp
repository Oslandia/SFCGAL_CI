// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/algorithm/collectionExtract.h>

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/MultiSolid.h>

#include <SFCGAL/detail/transform/ForceZ.h>

namespace SFCGAL {
namespace algorithm {

// If nothing has to be built, g will be moved to the result without
// copying and a new allocation. Otherwise, a new geometry is built and
// the old one is deleted.
std::unique_ptr<Geometry>
collectionExtractPolygons(std::unique_ptr<Geometry> g)
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

  MultiPolygon *ret_geo = new MultiPolygon;

  // copy each geometry
  for (size_t i = 0; i < coll.numGeometries(); ++i) {

    Geometry *gi = coll.geometryN(i).clone();

    switch (gi->geometryTypeId()) {
    case TYPE_TRIANGLE:
      ret_geo->addGeometry(Polygon(gi->as<Triangle>()));
      break;

    case TYPE_TRIANGULATEDSURFACE: {
      for (size_t j = 0; j < gi->numGeometries(); ++j) {
        ret_geo->addGeometry(Polygon(gi->geometryN(j).as<Triangle>()));
      }
    } break;

    case TYPE_POLYHEDRALSURFACE: {
      for (size_t j = 0; j < gi->numGeometries(); ++j) {
        ret_geo->addGeometry(gi->geometryN(j));
      }
    } break;

    case TYPE_POLYGON:
      ret_geo->addGeometry(*gi);
      break;

    default:
      // nothing
      break;
    }
  }

  return std::unique_ptr<Geometry>(ret_geo);
}

} // namespace algorithm
} // namespace SFCGAL
