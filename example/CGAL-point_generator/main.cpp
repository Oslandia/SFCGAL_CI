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
#include <algorithm>
#include <iostream>

#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/point_generators_3.h>

typedef CGAL::Cartesian<double>                  Kernel;
typedef Kernel::Point_2                          Point_2;
typedef Kernel::Point_3                          Point_3;
typedef CGAL::Creator_uniform_2<double, Point_2> Creator_2;
typedef CGAL::Creator_uniform_3<double, Point_3> Creator_3;

int
main()
{

  //-- generate points in disc
  {
    CGAL::Random_points_in_disc_2<Point_2, Creator_2> g(150.0);

    size_t n = 100U;
    for (size_t i = 0; i < n; i++) {
      Point_2 p = *(++g);
      std::cout << p << std::endl;
    }
  }

  //-- generate points in sphere
  {
    CGAL::Random_points_in_sphere_3<Point_3, Creator_3> g(150.0);

    size_t n = 100U;
    for (size_t i = 0; i < n; i++) {
      Point_3 p = *(++g);
      std::cout << p << std::endl;
    }
  }

  return 0;
}
