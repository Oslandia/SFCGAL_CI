// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/Geometry.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/algorithm/intersects.h>
#include <SFCGAL/detail/GeometrySet.h>
#include <SFCGAL/detail/GetPointsVisitor.h>
#include <SFCGAL/detail/algorithm/coversPoints.h>

namespace SFCGAL {
namespace detail {
namespace algorithm {
template <int Dim>
bool
_coversPoints(const Geometry &ga, const Geometry &gb)
{
  if (ga.isEmpty() || gb.isEmpty()) {
    return false;
  }

  GeometrySet<Dim> gsa(ga);

  // get all points of gb;
  detail::GetPointsVisitor visitor;
  gb.accept(visitor);

  for (detail::GetPointsVisitor::const_iterator it = visitor.points.begin();
       it != visitor.points.end(); ++it) {
    const Point *ppt = *it;

    // a geometry set of one point
    GeometrySet<Dim> gsp(*ppt);

    if (!SFCGAL::algorithm::intersects(gsp, gsa)) {
      return false;
    }
  }

  return true;
}

bool
coversPoints(const Geometry &ga, const Geometry &gb)
{
  return _coversPoints<2>(ga, gb);
}

bool
coversPoints3D(const Geometry &ga, const Geometry &gb)
{
  return _coversPoints<3>(ga, gb);
}
} // namespace algorithm
} // namespace detail
} // namespace SFCGAL
