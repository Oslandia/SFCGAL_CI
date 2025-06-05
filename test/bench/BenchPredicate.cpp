// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "../test_config.h"
#include "Bench.h"

#include <boost/test/unit_test.hpp>

#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/detail/GetPointsVisitor.h"
#include "SFCGAL/detail/generator/sierpinski.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_BenchPredicate)

double
randf()
{
  return double(rand()) / RAND_MAX;
}

BOOST_AUTO_TEST_CASE(testPointCompare)
{
  std::vector<CGAL::Point_3<Kernel>> points;
  for (size_t i = 0; i < 10000; i++) {
    points.push_back(CGAL::Point_3<Kernel>(randf(), randf(), randf()));
  }

  bench().start("xyz");

  for (size_t i = 0; i < points.size(); i++) {
    for (size_t j = 0; j < points.size(); j++) {
      points[i].x() < points[j].x();
    }
  }

  bench().stop();

  bench().start("compare");

  for (size_t i = 0; i < points.size(); i++) {
    for (size_t j = 0; j < points.size(); j++) {
      CGAL::compare_x<Kernel>(points[i], points[j]);
    }
  }

  bench().stop();
}

BOOST_AUTO_TEST_SUITE_END()
