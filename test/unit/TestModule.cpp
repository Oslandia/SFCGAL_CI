// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#define BOOST_TEST_MODULE UnitTestSFCGAL

#define BOOST_TEST_ALTERNATIVE_INIT_API

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

#include "SFCGAL/detail/tools/Log.h"

auto
init_unit_test_suite(int /*unused*/, char **const /*unused*/) -> test_suite *
{
  //	std::cerr << "init test suite" << std::endl;
  SFCGAL::Logger::get()->setLogLevel(SFCGAL::Logger::Info);
  return nullptr;
}
