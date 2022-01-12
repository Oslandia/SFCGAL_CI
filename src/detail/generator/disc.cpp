// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include <SFCGAL/detail/generator/disc.h>

#include <SFCGAL/LineString.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>

#include <cmath>
#include <memory>

namespace SFCGAL {
namespace generator {

///
///
///
auto
disc(const Point &center, const double &radius,
     const unsigned int &nQuadrantSegments) -> std::unique_ptr<Polygon>
{
  BOOST_ASSERT(nQuadrantSegments > 1);

  std::unique_ptr<LineString> exteriorRing(new LineString());

  double dTheta = M_PI_4 / nQuadrantSegments;

  for (size_t i = 0; i < nQuadrantSegments * 4; i++) {
    Kernel::Vector_2 p =
        center.toVector_2() +
        radius * Kernel::Vector_2(cos(i * dTheta), sin(i * dTheta));
    exteriorRing->addPoint(new Point(p.x(), p.y()));
  }

  exteriorRing->addPoint(exteriorRing->startPoint());

  return std::make_unique<Polygon>(exteriorRing.release());
}

} // namespace generator
} // namespace SFCGAL
