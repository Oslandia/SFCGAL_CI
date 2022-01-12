/**
 *   SFCGAL
 *
 *   Copyright (C) 2012-2013 Oslandia <infos@oslandia.com>
 *   Copyright (C) 2012-2013 IGN (http://www.ign.fr)
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Library General Public
 *   License as published by the Free Software Foundation; either
 *   version 2 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Library General Public License for more details.

 *   You should have received a copy of the GNU Library General Public
 *   License along with this library; if not, see
 <http://www.gnu.org/licenses/>.
 */

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
