// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/polygonSetToMultiPolygon.h"

#include <CGAL/Polygon_with_holes_2.h>

#include <list>

namespace SFCGAL::detail {

auto
polygonSetToMultiPolygon(const CGAL::Polygon_set_2<Kernel> &polygonSet)
    -> std::unique_ptr<MultiPolygon>
{
  std::list<Polygon_with_holes_2> res;
  polygonSet.polygons_with_holes(std::back_inserter(res));

  std::unique_ptr<MultiPolygon> result(new MultiPolygon);

  for (auto &re : res) {
    result->addGeometry(new Polygon(re));
  }

  return result;
}

} // namespace SFCGAL::detail
