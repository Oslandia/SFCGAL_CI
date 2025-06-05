// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/Geometry.h"
#include "SFCGAL/io/wkt.h"

#include "../../../test_config.h"

#include <boost/test/unit_test.hpp>

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_IntersectionTest)

//
// https://trac.osgeo.org/postgis/ticket/4157
BOOST_AUTO_TEST_CASE(test_postgis_4157)
{
  std::unique_ptr<Geometry> const g1(
      io::readWkt("POLYGON Z (("
                  "122395.299 489126.697 8.61546664325712,"
                  "122389.298 489128.73 8.55588025324629,"
                  "122391.489 489135.198 8.5526708028059,"
                  "122397.49 489133.165 8.61225719281685,"
                  "122395.299 489126.697 8.61546664325712))"));
  std::unique_ptr<Geometry> const g2(
      io::readWkt("POLYHEDRALSURFACE Z ((("
                  "122390.998245685 489133.068537491 0,"
                  "122391.003145022 489133.066423547 0,"
                  "122391.003145022 489133.066423547 10,"
                  "122390.998245685 489133.068537491 10,"
                  "122390.998245685 489133.068537491 0"
                  ")),(("
                  "122391.003145022 489133.066423547 0,"
                  "122383.269575402 489114.842869866 0,"
                  "122383.269575402 489114.842869866 10,"
                  "122391.003145022 489133.066423547 10,"
                  "122391.003145022 489133.066423547 0"
                  ")))"));

  algorithm::intersection3D(*g1, *g2);
}

BOOST_AUTO_TEST_SUITE_END()
