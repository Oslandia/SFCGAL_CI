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
#include <boost/test/unit_test.hpp>

#include <SFCGAL/Kernel.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/detail/io/WktWriter.h>
#include <SFCGAL/detail/transform/ForceZOrderPoints.h>
#include <SFCGAL/io/wkt.h>

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_transform_ForceZOrderPointsTest)

BOOST_AUTO_TEST_CASE(simple)
{
  std::unique_ptr<Geometry> g1 = io::readWkt("POLYGON((0 0,0 1,1 1,1 0,0 0))");

  const Polygon &p = g1->as<Polygon>();
  BOOST_CHECK(!p.isCounterClockWiseOriented());

  transform::ForceZOrderPoints forceZ;
  g1->accept(forceZ);

  BOOST_CHECK(g1->is3D());
  BOOST_CHECK(g1->as<Polygon>().isCounterClockWiseOriented());
}

BOOST_AUTO_TEST_SUITE_END()
