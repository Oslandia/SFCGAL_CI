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
#include <iostream>

#include <boost/test/unit_test.hpp>

#include <SFCGAL/Coordinate.h>
#include <SFCGAL/Kernel.h>
#include <SFCGAL/LineString.h>

using namespace SFCGAL;

using Vector_2  = Kernel::Vector_2;
using Vector_3  = Kernel::Vector_3;
using Point_2   = Kernel::Point_2;
using Point_3   = Kernel::Point_3;
using Segment_2 = Kernel::Segment_2;
using Segment_3 = Kernel::Segment_3;

// always after CGAL
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_KernelTest)

/**
 * 1 - 1/3 - 1/3 - 1/3 = 0
 */
BOOST_AUTO_TEST_CASE(testRobustArithmetric)
{
  Kernel::FT v = 1;
  v -= Kernel::FT(1) / Kernel::FT(3);
  v -= Kernel::FT(1) / Kernel::FT(3);
  v -= Kernel::FT(1) / Kernel::FT(3);
  BOOST_CHECK_EQUAL(v, 0);
}

/**
 * Serialize/Deserialize 1/3
 * @todo check with hugo's code
 */
BOOST_AUTO_TEST_CASE(testSerializeDeserialize)
{
  Kernel::FT a = 1;
  a /= 3;

  std::stringstream ss;
  ss << CGAL::exact(a);

  Kernel::FT b;
  ss >> b;
  BOOST_CHECK_EQUAL(a, b);
}

/**
 * 3 lines intersecting on POINT(1/3 1)
 */
BOOST_AUTO_TEST_CASE(testIntersectsRobutness)
{
  LineString ab(Point(0.0, 0.0), Point(1.0, 3.0));
  LineString cd(Point(0.0, 1.0), Point(1.0, 1.0));
  LineString ef(Point(-1.0, 3.0), Point(1.0, 0.0));

  // ab, cd
  CGAL::Object const abIcd_ = CGAL::intersection(
      Segment_2(ab.startPoint().toPoint_2(), ab.endPoint().toPoint_2()),
      Segment_2(cd.startPoint().toPoint_2(), cd.endPoint().toPoint_2()));
  const auto *abIcd = CGAL::object_cast<Point_2>(&abIcd_);
  BOOST_REQUIRE(abIcd != NULL);

  // would break robustness if construction history is lost
  Point const intersectionA(*abIcd);

  CGAL::Object const abIef_ = CGAL::intersection(
      intersectionA.toPoint_2(),
      Segment_2(ef.startPoint().toPoint_2(), ef.endPoint().toPoint_2()));
  const auto *abIef = CGAL::object_cast<Point_2>(&abIef_);
  BOOST_REQUIRE(abIef != NULL);

  Point const intersectionB(*abIef);

  BOOST_CHECK_EQUAL(intersectionA.toPoint_2(), intersectionB.toPoint_2());
}

BOOST_AUTO_TEST_SUITE_END()
