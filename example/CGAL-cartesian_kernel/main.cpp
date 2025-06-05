// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <CGAL/Cartesian.h>
#include <CGAL/squared_distance_2.h>

// declaration Kernel simple
typedef CGAL::Cartesian<double> K;
typedef K::Point_2              Point_2;

int
main()
{
  Point_2                       a(0.0, 0.0);
  Point_2                       b(3.0, 4.0);
  K::Compute_squared_distance_2 squared_distance;
  std::cout << squared_distance(a, b) << std::endl;
  return 0;
}
