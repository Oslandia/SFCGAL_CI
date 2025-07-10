// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/generator/sierpinski.h"

#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/Triangle.h"

namespace SFCGAL::generator {

auto
_sierpinski(const std::vector<Kernel::Triangle_2> &triangles)
    -> std::vector<Kernel::Triangle_2>
{
  std::vector<Kernel::Triangle_2> result;
  result.reserve(triangles.size() * 3);

  for (const auto &triangle : triangles) {
    const Kernel::Point_2 &a = triangle.vertex(0);
    const Kernel::Point_2 &b = triangle.vertex(1);
    const Kernel::Point_2 &c = triangle.vertex(2);

    Kernel::Point_2 const iAB = a + (b - a) / 2;
    Kernel::Point_2 const iBC = b + (c - b) / 2;
    Kernel::Point_2 const iCA = c + (a - c) / 2;

    result.emplace_back(a, iAB, iCA);
    result.emplace_back(b, iBC, iAB);
    result.emplace_back(c, iCA, iBC);
  }

  return result;
}

auto
sierpinski(const unsigned int &order) -> std::unique_ptr<MultiPolygon>
{
  std::vector<Kernel::Triangle_2> triangles;
  triangles.emplace_back(Kernel::Point_2(1.0, sqrt(3.0)),
                         Kernel::Point_2(2.0, 0.0), Kernel::Point_2(0.0, 0.0));

  for (unsigned int k = 0; k < order; k++) {
    triangles = _sierpinski(triangles);
  }

  std::unique_ptr<MultiPolygon> result(new MultiPolygon);

  for (auto &triangle : triangles) {
    result->addGeometry(Triangle(triangle).toPolygon());
  }

  return result;
}

} // namespace SFCGAL::generator
