// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/algorithm/Grid.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;
using namespace SFCGAL::algorithm;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_GridTest)

BOOST_AUTO_TEST_CASE(testMakeSquareGrid)
{
  std::unique_ptr<Geometry> extent =
      io::readWkt("POLYGON((0 0,10 0,10 10,0 10,0 0))");
  BOOST_REQUIRE(extent != nullptr);

  auto grid = makeSquareGrid(*extent, 5.0);
  BOOST_CHECK_EQUAL(grid.size(), 4U);

  auto gc = gridToGeometryCollection(grid);
  BOOST_CHECK_EQUAL(gc->numGeometries(), 4U);
}

BOOST_AUTO_TEST_CASE(testMakeRectangleGrid)
{
  std::unique_ptr<Geometry> extent =
      io::readWkt("POLYGON((0 0,12 0,12 8,0 8,0 0))");
  BOOST_REQUIRE(extent != nullptr);

  auto grid = makeRectangleGrid(*extent, 4.0, 4.0);
  BOOST_CHECK_EQUAL(grid.size(), 6U);

  auto gc = gridToGeometryCollection(grid);
  BOOST_CHECK_EQUAL(gc->numGeometries(), 6U);
}

BOOST_AUTO_TEST_CASE(testMakeHexagonGridFlatTop)
{
  std::unique_ptr<Geometry> extent =
      io::readWkt("POLYGON((0 0,50 0,50 50,0 50,0 0))");
  BOOST_REQUIRE(extent != nullptr);

  auto grid = makeHexagonGrid(*extent, 5.0, true);
  BOOST_CHECK(!grid.empty());

  auto gc = gridToGeometryCollection(grid);
  BOOST_CHECK_GT(gc->numGeometries(), 0U);
}

BOOST_AUTO_TEST_CASE(testMakeHexagonGridPointyTop)
{
  std::unique_ptr<Geometry> extent =
      io::readWkt("POLYGON((0 0,50 0,50 50,0 50,0 0))");
  BOOST_REQUIRE(extent != nullptr);

  auto grid = makeHexagonGrid(*extent, 5.0, false);
  BOOST_CHECK(!grid.empty());

  auto gc = gridToGeometryCollection(grid);
  BOOST_CHECK_GT(gc->numGeometries(), 0U);
}

BOOST_AUTO_TEST_CASE(testMakeTriangleGrid)
{
  std::unique_ptr<Geometry> extent =
      io::readWkt("POLYGON((0 0,20 0,20 20,0 20,0 0))");
  BOOST_REQUIRE(extent != nullptr);

  auto grid = makeTriangleGrid(*extent, 5.0);
  BOOST_CHECK(!grid.empty());

  auto gc = gridToGeometryCollection(grid);
  BOOST_CHECK_GT(gc->numGeometries(), 0U);
}

BOOST_AUTO_TEST_CASE(testMakeDiamondGrid)
{
  std::unique_ptr<Geometry> extent =
      io::readWkt("POLYGON((0 0,20 0,20 20,0 20,0 0))");
  BOOST_REQUIRE(extent != nullptr);

  auto grid = makeDiamondGrid(*extent, 5.0);
  BOOST_CHECK(!grid.empty());

  auto gc = gridToGeometryCollection(grid);
  BOOST_CHECK_GT(gc->numGeometries(), 0U);
}

BOOST_AUTO_TEST_CASE(testClipToExtent)
{
  std::unique_ptr<Geometry> extent =
      io::readWkt("POLYGON((0 0,10 0,10 10,0 10,0 0))");
  BOOST_REQUIRE(extent != nullptr);

  auto gridCover = makeSquareGrid(*extent, 5.0, GridClipMode::COVER_EXTENT);
  auto gridClip  = makeSquareGrid(*extent, 5.0, GridClipMode::CLIP_TO_EXTENT);

  // Both should have same count for exact fit
  BOOST_CHECK_EQUAL(gridCover.size(), gridClip.size());
}

BOOST_AUTO_TEST_CASE(testInvalidCellSize)
{
  std::unique_ptr<Geometry> extent =
      io::readWkt("POLYGON((0 0,10 0,10 10,0 10,0 0))");
  BOOST_REQUIRE(extent != nullptr);

  BOOST_CHECK_THROW(makeSquareGrid(*extent, 0.0), Exception);
  BOOST_CHECK_THROW(makeSquareGrid(*extent, -5.0), Exception);
}

BOOST_AUTO_TEST_SUITE_END()
