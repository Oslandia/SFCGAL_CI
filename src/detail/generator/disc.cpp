// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifdef _MSC_VER
  #define _USE_MATH_DEFINES
#endif

#include "SFCGAL/detail/generator/disc.h"

#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"

#include <cmath>
#include <cstddef>
#include <memory>

namespace SFCGAL::generator {

auto
disc(const Point &center, const double &radius,
     const unsigned int &nQuadrantSegments) -> std::unique_ptr<Polygon>
{
  BOOST_ASSERT(nQuadrantSegments > 1);

  std::unique_ptr<LineString> exteriorRing(new LineString());

  double const dTheta = M_PI_4 / nQuadrantSegments;

  for (size_t i = 0; i < static_cast<unsigned long>(nQuadrantSegments) * 4;
       i++) {
    Kernel::Vector_2 const p =
        center.toVector_2() +
        radius * Kernel::Vector_2(cos(i * dTheta), sin(i * dTheta));
    exteriorRing->addPoint(new Point(p.x(), p.y()));
  }

  exteriorRing->addPoint(exteriorRing->startPoint());

  return std::make_unique<Polygon>(exteriorRing.release());
}

} // namespace SFCGAL::generator
