// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/algorithm/roofGeneration.h"

using namespace SFCGAL;
using namespace SFCGAL::algorithm;

BOOST_AUTO_TEST_SUITE(RoofGenerationTests)

BOOST_AUTO_TEST_CASE(testCalculateRidgeHeight)
{
  // Test basic height calculation
  double height = calculateRidgeHeight(10.0, 30.0); // 10m distance, 30° slope
  BOOST_CHECK_CLOSE(height, 5.773, 0.1); // tan(30°) ≈ 0.577

  height = calculateRidgeHeight(10.0, 45.0); // 10m distance, 45° slope
  BOOST_CHECK_CLOSE(height, 10.0, 0.1); // tan(45°) = 1.0

  // Test edge cases
  BOOST_CHECK_THROW(calculateRidgeHeight(10.0, 0.0), Exception);
  BOOST_CHECK_THROW(calculateRidgeHeight(10.0, 90.0), Exception);
}

BOOST_AUTO_TEST_CASE(testCalculateHorizontalDistance)
{
  // Test basic distance calculation
  double distance = calculateHorizontalDistance(10.0, 45.0); // 10m height, 45° slope
  BOOST_CHECK_CLOSE(distance, 10.0, 0.1); // cot(45°) = 1.0

  distance = calculateHorizontalDistance(5.773, 30.0); // height for 30° slope
  BOOST_CHECK_CLOSE(distance, 10.0, 0.1); // cot(30°) ≈ 1.732

  // Test edge cases
  BOOST_CHECK_THROW(calculateHorizontalDistance(10.0, 0.0), Exception);
  BOOST_CHECK_THROW(calculateHorizontalDistance(10.0, 90.0), Exception);
}

BOOST_AUTO_TEST_CASE(testGeneratePitchedRoofBasic)
{
  // Create a simple rectangular footprint
  std::vector<Point> points = {
      Point(0, 0, 0), Point(10, 0, 0), Point(10, 6, 0), Point(0, 6, 0), Point(0, 0, 0)
  };
  LineString ring(points);
  Polygon footprint(ring);

  // Create a ridge line through the center
  LineString ridgeLine(Point(2, 3, 0), Point(8, 3, 0));

  // Generate pitched roof
  auto roof = generatePitchedRoof(footprint, ridgeLine, 30.0);

  BOOST_CHECK(roof != nullptr);
  BOOST_CHECK(!roof->isEmpty());
  BOOST_CHECK(roof->numPatches() > 0);
}

BOOST_AUTO_TEST_CASE(testGenerateGableRoof)
{
  // Create a simple rectangular footprint
  std::vector<Point> points = {
      Point(0, 0, 0), Point(10, 0, 0), Point(10, 6, 0), Point(0, 6, 0), Point(0, 0, 0)
  };
  LineString ring(points);
  Polygon footprint(ring);

  // Create a ridge line along the center
  LineString ridgeLine(Point(0, 3, 0), Point(10, 3, 0));

  // Generate gable roof
  auto roof = generateGableRoof(footprint, ridgeLine, 30.0);

  BOOST_CHECK(roof != nullptr);
  BOOST_CHECK(!roof->isEmpty());
  BOOST_CHECK(roof->numPatches() > 0);
}

BOOST_AUTO_TEST_CASE(testGenerateSkillionRoof)
{
  // Create a simple rectangular footprint
  std::vector<Point> points = {
      Point(0, 0, 0), Point(10, 0, 0), Point(10, 6, 0), Point(0, 6, 0), Point(0, 0, 0)
  };
  LineString ring(points);
  Polygon footprint(ring);

  // Create a ridge line on one edge
  LineString ridgeLine(Point(0, 0, 0), Point(10, 0, 0));

  // Generate skillion roof
  auto roof = generateSkillionRoof(footprint, ridgeLine, 20.0);

  BOOST_CHECK(roof != nullptr);
  BOOST_CHECK(!roof->isEmpty());
  BOOST_CHECK(roof->numPatches() > 0);
}

BOOST_AUTO_TEST_CASE(testGenerateRoofWithParameters)
{
  // Create a simple rectangular footprint
  std::vector<Point> points = {
      Point(0, 0, 0), Point(10, 0, 0), Point(10, 6, 0), Point(0, 6, 0), Point(0, 0, 0)
  };
  LineString ring(points);
  Polygon footprint(ring);

  // Create a ridge line
  LineString ridgeLine(Point(2, 3, 0), Point(8, 3, 0));

  // Test different roof types
  RoofParameters params;

  // Test pitched roof
  params.type = RoofType::PITCHED;
  params.slopeAngle = 30.0;
  auto pitchedRoof = generateRoof(footprint, ridgeLine, params);
  BOOST_CHECK(pitchedRoof != nullptr);
  BOOST_CHECK(!pitchedRoof->isEmpty());

  // Test gable roof
  params.type = RoofType::GABLE;
  params.slopeAngle = 35.0;
  auto gableRoof = generateRoof(footprint, ridgeLine, params);
  BOOST_CHECK(gableRoof != nullptr);
  BOOST_CHECK(!gableRoof->isEmpty());

  // Test hipped roof (should use existing implementation)
  params.type = RoofType::HIPPED;
  params.roofHeight = 5.0;
  auto hippedRoof = generateRoof(footprint, ridgeLine, params);
  BOOST_CHECK(hippedRoof != nullptr);
  BOOST_CHECK(!hippedRoof->isEmpty());
}

BOOST_AUTO_TEST_CASE(testInvalidSlopeAngles)
{
  std::vector<Point> points = {
      Point(0, 0, 0), Point(10, 0, 0), Point(10, 6, 0), Point(0, 6, 0), Point(0, 0, 0)
  };
  LineString ring(points);
  Polygon footprint(ring);

  LineString ridgeLine(Point(2, 3, 0), Point(8, 3, 0));

  // Test invalid slope angles
  BOOST_CHECK_THROW(generatePitchedRoof(footprint, ridgeLine, 0.0), Exception);
  BOOST_CHECK_THROW(generatePitchedRoof(footprint, ridgeLine, 90.0), Exception);
  BOOST_CHECK_THROW(generatePitchedRoof(footprint, ridgeLine, -10.0), Exception);
  BOOST_CHECK_THROW(generatePitchedRoof(footprint, ridgeLine, 100.0), Exception);
}

BOOST_AUTO_TEST_CASE(testEmptyGeometry)
{
  Polygon emptyFootprint;
  LineString ridgeLine(Point(0, 0, 0), Point(1, 0, 0));

  // Should handle empty geometries gracefully
  auto roof = generatePitchedRoof(emptyFootprint, ridgeLine, 30.0);
  BOOST_CHECK(roof != nullptr);
  // The result behavior for empty input depends on implementation
}

BOOST_AUTO_TEST_CASE(testRidgePositions)
{
  std::vector<Point> points = {
      Point(0, 0, 0), Point(10, 0, 0), Point(10, 6, 0), Point(0, 6, 0), Point(0, 0, 0)
  };
  LineString ring(points);
  Polygon footprint(ring);

  // Interior ridge
  LineString interiorRidge(Point(2, 3, 0), Point(8, 3, 0));
  auto roofInterior = generatePitchedRoof(footprint, interiorRidge, 30.0,
                                         RidgePosition::INTERIOR);
  BOOST_CHECK(roofInterior != nullptr);
  BOOST_CHECK(!roofInterior->isEmpty());

  // Edge ridge
  LineString edgeRidge(Point(0, 0, 0), Point(10, 0, 0));
  auto roofEdge = generatePitchedRoof(footprint, edgeRidge, 30.0,
                                     RidgePosition::EDGE);
  BOOST_CHECK(roofEdge != nullptr);
  BOOST_CHECK(!roofEdge->isEmpty());

  // Exterior ridge
  LineString exteriorRidge(Point(-1, 3, 0), Point(11, 3, 0));
  auto roofExterior = generatePitchedRoof(footprint, exteriorRidge, 30.0,
                                         RidgePosition::EXTERIOR);
  BOOST_CHECK(roofExterior != nullptr);
  BOOST_CHECK(!roofExterior->isEmpty());
}

BOOST_AUTO_TEST_CASE(testComplexPolygon)
{
  // Test with a more complex polygon (L-shaped)
  std::vector<Point> points = {
      Point(0, 0, 0), Point(6, 0, 0), Point(6, 4, 0), Point(3, 4, 0),
      Point(3, 6, 0), Point(0, 6, 0), Point(0, 0, 0)
  };
  LineString ring(points);
  Polygon footprint(ring);

  LineString ridgeLine(Point(1, 3, 0), Point(5, 3, 0));

  // Should handle complex polygons
  auto roof = generatePitchedRoof(footprint, ridgeLine, 25.0);
  BOOST_CHECK(roof != nullptr);
  // Complex polygons may result in more patches
}

BOOST_AUTO_TEST_CASE(testBuildingHeights)
{
  // Test new building height and roof height parameters
  std::vector<Point> points = {
      Point(0, 0, 0), Point(10, 0, 0), Point(10, 6, 0), Point(0, 6, 0), Point(0, 0, 0)
  };
  LineString ring(points);
  Polygon footprint(ring);

  LineString ridgeLine(Point(0, 3, 0), Point(10, 3, 0));

  // Test pitched roof with building height
  auto building1 = generatePitchedRoof(footprint, ridgeLine, 3.0, 2.0, 30.0);
  BOOST_CHECK(building1 != nullptr);

  // Test gable roof with building height
  auto building2 = generateGableRoof(footprint, ridgeLine, 4.0, 3.0, 25.0);
  BOOST_CHECK(building2 != nullptr);

  // Test skillion roof with building height
  auto building3 = generateSkillionRoof(footprint, ridgeLine, 2.5, 1.5, 20.0);
  BOOST_CHECK(building3 != nullptr);

  // Test parameter validation
  BOOST_CHECK_THROW(generatePitchedRoof(footprint, ridgeLine, -1.0, 2.0, 30.0), Exception);
  BOOST_CHECK_THROW(generatePitchedRoof(footprint, ridgeLine, 3.0, 0.0, 30.0), Exception);
}

BOOST_AUTO_TEST_CASE(testRoofParameters)
{
  std::vector<Point> points = {
      Point(0, 0, 0), Point(10, 0, 0), Point(10, 6, 0), Point(0, 6, 0), Point(0, 0, 0)
  };
  LineString ring(points);
  Polygon footprint(ring);

  LineString ridgeLine(Point(0, 3, 0), Point(10, 3, 0));

  // Test new RoofParameters structure
  RoofParameters params;
  params.type = RoofType::PITCHED;
  params.buildingHeight = 5.0;
  params.roofHeight = 3.0;
  params.slopeAngle = 35.0;
  params.generateSolid = true;

  auto building = generateBuildingWithRoof(footprint, ridgeLine, params);
  BOOST_CHECK(building != nullptr);

  // Test flat roof with new parameters
  params.type = RoofType::FLAT;
  params.buildingHeight = 4.0;
  params.roofHeight = 1.0;
  auto flatBuilding = generateBuildingWithRoof(footprint, ridgeLine, params);
  BOOST_CHECK(flatBuilding != nullptr);

  // Test hipped roof with building height
  params.type = RoofType::HIPPED;
  params.buildingHeight = 6.0;
  params.roofHeight = 4.0;
  auto hippedBuilding = generateBuildingWithRoof(footprint, ridgeLine, params);
  BOOST_CHECK(hippedBuilding != nullptr);
}

BOOST_AUTO_TEST_CASE(testSolidGeneration)
{
  std::vector<Point> points = {
      Point(0, 0, 0), Point(5, 0, 0), Point(5, 5, 0), Point(0, 5, 0), Point(0, 0, 0)
  };
  LineString ring(points);
  Polygon footprint(ring);

  LineString ridgeLine(Point(0, 2.5, 0), Point(5, 2.5, 0));

  // Test that we can generate Solid geometry
  RoofParameters params;
  params.type = RoofType::GABLE;
  params.buildingHeight = 3.0;
  params.roofHeight = 2.0;
  params.slopeAngle = 30.0;
  params.generateSolid = true;

  auto building = generateBuildingWithRoof(footprint, ridgeLine, params);
  BOOST_CHECK(building != nullptr);

  // The result should be a valid geometry (Solid or PolyhedralSurface)
  BOOST_CHECK(!building->isEmpty());
}

BOOST_AUTO_TEST_CASE(testBackwardCompatibility)
{
  // Test that existing code still works with new RoofParameters
  std::vector<Point> points = {
      Point(0, 0, 0), Point(10, 0, 0), Point(10, 6, 0), Point(0, 6, 0), Point(0, 0, 0)
  };
  LineString ring(points);
  Polygon footprint(ring);

  LineString ridgeLine(Point(0, 3, 0), Point(10, 3, 0));

  // Old-style parameters should still work
  RoofParameters params;
  params.type = RoofType::PITCHED;
  params.roofHeight = 3.0;  // This should be used for backward compatibility
  params.slopeAngle = 30.0;

  auto roof = generateRoof(footprint, ridgeLine, params);
  BOOST_CHECK(roof != nullptr);
  BOOST_CHECK(!roof->isEmpty());
}

BOOST_AUTO_TEST_SUITE_END()