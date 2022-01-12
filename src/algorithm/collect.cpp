// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/algorithm/collect.h>

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/MultiSolid.h>

namespace SFCGAL {
namespace algorithm {
std::unique_ptr<Geometry>
collect(const Geometry &ga, const Geometry &gb)
{
  if (ga.geometryTypeId() == gb.geometryTypeId()) {
    if (ga.geometryTypeId() == TYPE_POINT) {
      MultiPoint *mp = new MultiPoint;
      mp->addGeometry(ga);
      mp->addGeometry(gb);
      return std::unique_ptr<Geometry>(mp);
    } else if (ga.geometryTypeId() == TYPE_LINESTRING) {
      MultiLineString *mls = new MultiLineString();
      mls->addGeometry(ga);
      mls->addGeometry(gb);
      return std::unique_ptr<Geometry>(mls);
    } else if (ga.geometryTypeId() == TYPE_POLYGON) {
      MultiPolygon *mp = new MultiPolygon();
      mp->addGeometry(ga);
      mp->addGeometry(gb);
      return std::unique_ptr<Geometry>(mp);
    } else if (ga.geometryTypeId() == TYPE_SOLID) {
      MultiSolid *mp = new MultiSolid();
      mp->addGeometry(ga);
      mp->addGeometry(gb);
      return std::unique_ptr<Geometry>(mp);
    }
  }

  // else
  GeometryCollection *coll = new GeometryCollection();
  coll->addGeometry(ga);
  coll->addGeometry(gb);
  return std::unique_ptr<Geometry>(coll);
}
} // namespace algorithm
} // namespace SFCGAL
