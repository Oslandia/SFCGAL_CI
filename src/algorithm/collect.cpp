// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/algorithm/collect.h>

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/MultiSolid.h>

namespace SFCGAL::algorithm {
auto
collect(const Geometry &ga, const Geometry &gb) -> std::unique_ptr<Geometry>
{
  if (ga.geometryTypeId() == gb.geometryTypeId()) {
    if (ga.geometryTypeId() == TYPE_POINT) {
      auto *mp = new MultiPoint;
      mp->addGeometry(ga);
      mp->addGeometry(gb);
      return std::unique_ptr<Geometry>(mp);
    } if (ga.geometryTypeId() == TYPE_LINESTRING) {
      auto *mls = new MultiLineString();
      mls->addGeometry(ga);
      mls->addGeometry(gb);
      return std::unique_ptr<Geometry>(mls);
    } if (ga.geometryTypeId() == TYPE_POLYGON) {
      auto *mp = new MultiPolygon();
      mp->addGeometry(ga);
      mp->addGeometry(gb);
      return std::unique_ptr<Geometry>(mp);
    } else if (ga.geometryTypeId() == TYPE_SOLID) {
      auto *mp = new MultiSolid();
      mp->addGeometry(ga);
      mp->addGeometry(gb);
      return std::unique_ptr<Geometry>(mp);
    }
  }

  // else
  auto *coll = new GeometryCollection();
  coll->addGeometry(ga);
  coll->addGeometry(gb);
  return std::unique_ptr<Geometry>(coll);
}
} // namespace SFCGAL
