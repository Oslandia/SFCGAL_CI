// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/generator/hoch.h"

#include "SFCGAL/LineString.h"
#include "SFCGAL/Polygon.h"

namespace SFCGAL::generator {

auto
_hoch(const std::vector<Kernel::Vector_2> &points)
    -> std::vector<Kernel::Vector_2>
{
  std::vector<Kernel::Vector_2> result;
  result.reserve(points.size() * 2);
  size_t const numPoints = points.size();

  for (size_t i = 0; i < numPoints; i++) {
    const Kernel::Vector_2 &a = points[i];
    const Kernel::Vector_2 &b = points[(i + 1) % numPoints];

    Kernel::Vector_2 const ab = b - a;
    Kernel::Vector_2 const normal(-ab.y(), ab.x());

    result.push_back(a);
    result.push_back(a + ab / 3);
    result.push_back(a + ab / 2 + normal * sqrt(3.0) / 6.0);
    result.push_back(a + (ab * 2) / 3);
  }

  return result;
}

auto
hoch(const unsigned int &order) -> std::unique_ptr<Polygon>
{
  std::vector<Kernel::Vector_2> points;
  points.emplace_back(1.0, sqrt(3.0));
  points.emplace_back(2.0, 0.0);
  points.emplace_back(0.0, 0.0);

  for (unsigned int k = 0; k < order; k++) {
    points = _hoch(points);
  }

  std::unique_ptr<Polygon>    result(new Polygon());
  std::unique_ptr<LineString> ring(new LineString());

  for (auto &point : points) {
    ring->addPoint(new Point(point.x(), point.y()));
  }

  ring->addPoint(new Point(points.front().x(), points.front().y()));

  result->setExteriorRing(ring.release());

  return result;
}

} // namespace SFCGAL::generator
