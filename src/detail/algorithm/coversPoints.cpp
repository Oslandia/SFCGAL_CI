// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/algorithm/coversPoints.h"
#include "SFCGAL/Geometry.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/algorithm/intersects.h"
#include "SFCGAL/detail/GeometrySet.h"
#include "SFCGAL/detail/GetPointsVisitor.h"

namespace SFCGAL::detail::algorithm {
/// @private
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

  for (const auto *ppt : visitor.points) {
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
} // namespace SFCGAL::detail::algorithm
