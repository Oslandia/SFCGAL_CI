// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/Geometry.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/algorithm/intersects.h>
#include <SFCGAL/detail/GeometrySet.h>
#include <SFCGAL/detail/GetPointsVisitor.h>
#include <SFCGAL/detail/algorithm/coversPoints.h>

namespace SFCGAL::detail::algorithm {
template <int Dim>
auto
_coversPoints(const Geometry &ga, const Geometry &gb) -> bool
{
  if (ga.isEmpty() || gb.isEmpty()) {
    return false;
  }

  GeometrySet<Dim> const gsa(ga);

  // get all points of gb;
  detail::GetPointsVisitor visitor;
  gb.accept(visitor);

  for (auto it = visitor.points.begin();
       it != visitor.points.end(); ++it) {
    const Point *ppt = *it;

    // a geometry set of one point
    GeometrySet<Dim> const gsp(*ppt);

    if (!SFCGAL::algorithm::intersects(gsp, gsa)) {
      return false;
    }
  }

  return true;
}

auto
coversPoints(const Geometry &ga, const Geometry &gb) -> bool
{
  return _coversPoints<2>(ga, gb);
}

auto
coversPoints3D(const Geometry &ga, const Geometry &gb) -> bool
{
  return _coversPoints<3>(ga, gb);
}
} // namespace SFCGAL
