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

#include <SFCGAL/detail/polygonSetToMultiPolygon.h>

#include <CGAL/Polygon_with_holes_2.h>

#include <list>

namespace SFCGAL {
namespace detail {

///
///
///
std::unique_ptr<MultiPolygon>
polygonSetToMultiPolygon(const CGAL::Polygon_set_2<Kernel> &polygonSet)
{
  typedef CGAL::Polygon_with_holes_2<Kernel> Polygon_with_holes_2;

  std::list<Polygon_with_holes_2> res;
  polygonSet.polygons_with_holes(std::back_inserter(res));

  std::unique_ptr<MultiPolygon> result(new MultiPolygon);

  for (std::list<Polygon_with_holes_2>::const_iterator it = res.begin();
       it != res.end(); ++it) {
    result->addGeometry(new Polygon(*it));
  }

  return result;
}

} // namespace detail
} // namespace SFCGAL
