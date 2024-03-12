// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/collectionToMulti.h"

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
auto
collectionToMulti(std::unique_ptr<Geometry> g) -> std::unique_ptr<Geometry>
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

  bool has2d = false;
  bool has3d = false;

  for (size_t i = 0; i < coll.numGeometries(); ++i) {
    const Geometry &gi = coll.geometryN(i);

    if (!has3d && gi.is3D()) {
      has3d = true;
    }

    if (!has2d && !gi.is3D()) {
      has2d = true;
    }

    if (!gi.isEmpty() && (gi.geometryTypeId() != TYPE_POLYGON) &&
        (gi.geometryTypeId() != TYPE_TRIANGLE) &&
        (gi.geometryTypeId() != TYPE_POLYHEDRALSURFACE) &&
        (gi.geometryTypeId() != TYPE_TRIANGULATEDSURFACE)) {
      // it contains a bad type, abort
      return g;
    }
  }

  bool const force3d = has2d && has3d;

  auto *ret_geo = new MultiPolygon;

  // copy each geometry
  for (size_t i = 0; i < coll.numGeometries(); ++i) {

    Geometry *gi = coll.geometryN(i).clone();

    if (force3d && !gi->is3D()) {
      transform::ForceZ forceZ;
      gi->accept(forceZ);
    }

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

    case TYPE_GEOMETRYCOLLECTION:

      // do not include empty geometrycollection
      if (gi->isEmpty()) {
        continue;
      }
      ret_geo->addGeometry(*gi);
      break;

    default:
      ret_geo->addGeometry(*gi);
    }
  }

  return std::unique_ptr<Geometry>(ret_geo);
}

} // namespace SFCGAL::algorithm
