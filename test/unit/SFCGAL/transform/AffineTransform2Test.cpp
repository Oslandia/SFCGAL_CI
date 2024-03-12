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

#include "SFCGAL/Kernel.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/detail/transform/AffineTransform2.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_transform_AffineTransform2Test)

BOOST_AUTO_TEST_CASE(simpleTranslate)
{
  LineString g(Point(0.0, 0.0), Point(1.0, 1.0));

  CGAL::Aff_transformation_2<Kernel> const affine(
      CGAL::TRANSLATION, CGAL::Vector_2<Kernel>(1.0, 2.0));

  transform::AffineTransform2 transform(affine);
  g.accept(transform);

  BOOST_CHECK_EQUAL(g.pointN(0).asText(5), "POINT(1.00000 2.00000)");
  BOOST_CHECK_EQUAL(g.pointN(1).asText(5), "POINT(2.00000 3.00000)");
}

BOOST_AUTO_TEST_SUITE_END()
