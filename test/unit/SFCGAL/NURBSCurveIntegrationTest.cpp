// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/Envelope.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/GeometryVisitor.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/algorithm/BoundaryVisitor.h"
#include "SFCGAL/algorithm/distance.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/detail/transform/AffineTransform2.h"
#include "SFCGAL/detail/transform/AffineTransform3.h"
#include "SFCGAL/io/wkt.h"
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Aff_transformation_3.h>
#include <cmath>
#include <memory>
#include <sstream>

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_NURBSCurveIntegrationTest)

//-- Helper functions for integration tests

std::vector<NURBSCurve::FT>
convertWeights(const std::vector<double> &doubleWeights)
{
  std::vector<NURBSCurve::FT> weights;
  weights.reserve(doubleWeights.size());
  for (const auto &weight : doubleWeights) {
    weights.emplace_back(weight);
  }
  return weights;
}

bool
isNearlyEqual(const Point &firstPoint, const Point &secondPoint,
              double tolerance = 1e-10)
{
  return algorithm::distance(firstPoint, secondPoint) <
         NURBSCurve::FT(tolerance);
}

std::vector<Point>
createTestPoints(bool is3D = false, bool isMeasured = false)
{
  std::vector<Point> points;

  if (is3D && isMeasured) {
    points.emplace_back(0.0, 0.0, 0.0, 10.0);
    points.emplace_back(2.0, 4.0, 2.0, 20.0);
    points.emplace_back(6.0, 4.0, -1.0, 30.0);
    points.emplace_back(8.0, 0.0, 0.0, 40.0);
  } else if (is3D) {
    points.emplace_back(0.0, 0.0, 0.0);
    points.emplace_back(2.0, 4.0, 2.0);
    points.emplace_back(6.0, 4.0, -1.0);
    points.emplace_back(8.0, 0.0, 0.0);
  } else if (isMeasured) {
    Point firstPoint(0.0, 0.0);
    firstPoint.setM(10.0);
    Point secondPoint(2.0, 4.0);
    secondPoint.setM(20.0);
    Point thirdPoint(6.0, 4.0);
    thirdPoint.setM(30.0);
    Point fourthPoint(8.0, 0.0);
    fourthPoint.setM(40.0);
    points.push_back(firstPoint);
    points.push_back(secondPoint);
    points.push_back(thirdPoint);
    points.push_back(fourthPoint);
  } else {
    points.emplace_back(0.0, 0.0);
    points.emplace_back(2.0, 4.0);
    points.emplace_back(6.0, 4.0);
    points.emplace_back(8.0, 0.0);
  }

  return points;
}

//-- Visitor pattern integration tests

class TestGeometryVisitor : public GeometryVisitor {
public:
  std::string lastVisited;
  size_t      visitCount = 0;

  void
  visit(Point &geom) override
  {
    lastVisited = "Point";
    visitCount++;
  }
  void
  visit(LineString &geom) override
  {
    lastVisited = "LineString";
    visitCount++;
  }
  void
  visit(Polygon &geom) override
  {
    lastVisited = "Polygon";
    visitCount++;
  }
  void
  visit(Triangle &geom) override
  {
    lastVisited = "Triangle";
    visitCount++;
  }
  void
  visit(Solid &geom) override
  {
    lastVisited = "Solid";
    visitCount++;
  }
  void
  visit(MultiPoint &geom) override
  {
    lastVisited = "MultiPoint";
    visitCount++;
  }
  void
  visit(MultiLineString &geom) override
  {
    lastVisited = "MultiLineString";
    visitCount++;
  }
  void
  visit(MultiPolygon &geom) override
  {
    lastVisited = "MultiPolygon";
    visitCount++;
  }
  void
  visit(MultiSolid &geom) override
  {
    lastVisited = "MultiSolid";
    visitCount++;
  }
  void
  visit(GeometryCollection &geom) override
  {
    lastVisited = "GeometryCollection";
    visitCount++;
  }
  void
  visit(PolyhedralSurface &geom) override
  {
    lastVisited = "PolyhedralSurface";
    visitCount++;
  }
  void
  visit(TriangulatedSurface &geom) override
  {
    lastVisited = "TriangulatedSurface";
    visitCount++;
  }
  void
  visit(NURBSCurve &geom) override
  {
    lastVisited = "NURBSCurve";
    visitCount++;
  }
};

class TestConstGeometryVisitor : public ConstGeometryVisitor {
public:
  std::string lastVisited;
  size_t      visitCount = 0;

  void
  visit(const Point &geom) override
  {
    lastVisited = "Point";
    visitCount++;
  }
  void
  visit(const LineString &geom) override
  {
    lastVisited = "LineString";
    visitCount++;
  }
  void
  visit(const Polygon &geom) override
  {
    lastVisited = "Polygon";
    visitCount++;
  }
  void
  visit(const Triangle &geom) override
  {
    lastVisited = "Triangle";
    visitCount++;
  }
  void
  visit(const Solid &geom) override
  {
    lastVisited = "Solid";
    visitCount++;
  }
  void
  visit(const MultiPoint &geom) override
  {
    lastVisited = "MultiPoint";
    visitCount++;
  }
  void
  visit(const MultiLineString &geom) override
  {
    lastVisited = "MultiLineString";
    visitCount++;
  }
  void
  visit(const MultiPolygon &geom) override
  {
    lastVisited = "MultiPolygon";
    visitCount++;
  }
  void
  visit(const MultiSolid &geom) override
  {
    lastVisited = "MultiSolid";
    visitCount++;
  }
  void
  visit(const GeometryCollection &geom) override
  {
    lastVisited = "GeometryCollection";
    visitCount++;
  }
  void
  visit(const PolyhedralSurface &geom) override
  {
    lastVisited = "PolyhedralSurface";
    visitCount++;
  }
  void
  visit(const TriangulatedSurface &geom) override
  {
    lastVisited = "TriangulatedSurface";
    visitCount++;
  }
  void
  visit(const NURBSCurve &geom) override
  {
    lastVisited = "NURBSCurve";
    visitCount++;
  }
};

BOOST_AUTO_TEST_CASE(testGeometryVisitorIntegration)
{
  auto       controlPoints = createTestPoints();
  NURBSCurve curve(controlPoints, 3);

  TestGeometryVisitor visitor;
  curve.accept(visitor);

  BOOST_CHECK_EQUAL(visitor.lastVisited, "NURBSCurve");
  BOOST_CHECK_EQUAL(visitor.visitCount, 1U);

  // Test const visitor
  TestConstGeometryVisitor constVisitor;
  const NURBSCurve        &constCurve = curve;
  constCurve.accept(constVisitor);

  BOOST_CHECK_EQUAL(constVisitor.lastVisited, "NURBSCurve");
  BOOST_CHECK_EQUAL(constVisitor.visitCount, 1U);
}

BOOST_AUTO_TEST_CASE(testGeometryCollectionWithNURBS)
{
  auto       controlPoints = createTestPoints();
  NURBSCurve curve(controlPoints, 2);

  GeometryCollection collection;
  collection.addGeometry(curve.clone());
  collection.addGeometry(new Point(10.0, 10.0));

  BOOST_CHECK_EQUAL(collection.numGeometries(), 2U);
  BOOST_CHECK(collection.geometryN(0).is<NURBSCurve>());
  BOOST_CHECK(collection.geometryN(1).is<Point>());

  // Test visitor on collection
  TestGeometryVisitor visitor;
  collection.accept(visitor);

  // Should visit the collection itself
  BOOST_CHECK_EQUAL(visitor.lastVisited, "GeometryCollection");
}

//-- Transform integration tests

BOOST_AUTO_TEST_CASE(testTransformNURBSCurve2D)
{
  auto       controlPoints = createTestPoints();
  NURBSCurve originalCurve(controlPoints, 2);

  // Test translation using AffineTransform2
  Kernel::Vector_2                   translationVector(5.0, 3.0);
  CGAL::Aff_transformation_2<Kernel> translationMatrix(CGAL::TRANSLATION,
                                                       translationVector);
  transform::AffineTransform2        translation(translationMatrix);

  NURBSCurve transformedCurve = originalCurve;
  transformedCurve.accept(translation);

  // Check that all control points were translated
  for (size_t idx = 0; idx < originalCurve.numControlPoints(); ++idx) {
    Point original    = originalCurve.controlPointN(idx);
    Point transformed = transformedCurve.controlPointN(idx);

    BOOST_CHECK_CLOSE(CGAL::to_double(transformed.x()),
                      CGAL::to_double(original.x()) + 5.0, 1e-10);
    BOOST_CHECK_CLOSE(CGAL::to_double(transformed.y()),
                      CGAL::to_double(original.y()) + 3.0, 1e-10);
  }

  // Curve properties should be preserved
  BOOST_CHECK_EQUAL(transformedCurve.degree(), originalCurve.degree());
  BOOST_CHECK_EQUAL(transformedCurve.numControlPoints(),
                    originalCurve.numControlPoints());
  BOOST_CHECK_EQUAL(transformedCurve.isRational(), originalCurve.isRational());
}

BOOST_AUTO_TEST_CASE(testTransformNURBSCurve3D)
{
  auto       controlPoints = createTestPoints(true, false);
  auto       weights       = convertWeights({1.0, 2.0, 1.5, 1.0});
  NURBSCurve originalCurve(controlPoints, weights, 2);

  // Test scaling and rotation using AffineTransform3
  double                             angle = M_PI / 4.0; // 45 degrees
  CGAL::Aff_transformation_3<Kernel> rotationMatrix(
      std::cos(angle), -std::sin(angle), 0.0, 0.0, std::sin(angle),
      std::cos(angle), 0.0, 0.0, 0.0, 0.0, 2.0, 0.0); // Scale Z by 2

  transform::AffineTransform3 rotation(rotationMatrix);

  NURBSCurve transformedCurve = originalCurve;
  transformedCurve.accept(rotation);

  BOOST_CHECK(transformedCurve.is3D());
  BOOST_CHECK(transformedCurve.isRational());

  // Check first control point transformation
  Point original    = originalCurve.controlPointN(0);
  Point transformed = transformedCurve.controlPointN(0);

  // For point (0,0,0), rotation should still give (0,0,0)
  BOOST_CHECK_SMALL(CGAL::to_double(transformed.x()), 1e-10);
  BOOST_CHECK_SMALL(CGAL::to_double(transformed.y()), 1e-10);
  BOOST_CHECK_SMALL(CGAL::to_double(transformed.z()), 1e-10);

  // Check that weights are preserved
  for (size_t idx = 0; idx < originalCurve.numControlPoints(); ++idx) {
    BOOST_CHECK_EQUAL(originalCurve.weight(idx), transformedCurve.weight(idx));
  }
}

BOOST_AUTO_TEST_CASE(testTransformDirectApplication)
{
  auto       controlPoints = createTestPoints();
  NURBSCurve originalCurve(controlPoints, 2);

  // Test direct application of transformation
  Kernel::Vector_2                   translationVector(10.0, 20.0);
  CGAL::Aff_transformation_2<Kernel> translationMatrix(CGAL::TRANSLATION,
                                                       translationVector);
  transform::AffineTransform2        translation(translationMatrix);

  NURBSCurve transformedCurve = originalCurve;
  transformedCurve.accept(translation);

  // Check translation was applied
  Point original    = originalCurve.controlPointN(0);
  Point transformed = transformedCurve.controlPointN(0);

  BOOST_CHECK_CLOSE(CGAL::to_double(transformed.x()),
                    CGAL::to_double(original.x()) + 10.0, 1e-10);
  BOOST_CHECK_CLOSE(CGAL::to_double(transformed.y()),
                    CGAL::to_double(original.y()) + 20.0, 1e-10);
}

//-- Algorithm integration tests

BOOST_AUTO_TEST_CASE(testEnvelopeAlgorithm)
{
  auto       controlPoints = createTestPoints();
  NURBSCurve curve(controlPoints, 3);

  Envelope envelope = curve.envelope();

  BOOST_CHECK(!envelope.isEmpty());
  BOOST_CHECK(!envelope.is3D()); // 2D curve

  // Envelope should contain all control points
  for (const auto &point : controlPoints) {
    BOOST_CHECK(envelope.xMin() <= point.x());
    BOOST_CHECK(envelope.xMax() >= point.x());
    BOOST_CHECK(envelope.yMin() <= point.y());
    BOOST_CHECK(envelope.yMax() >= point.y());
  }

  // Test 3D envelope
  auto       controlPoints3D = createTestPoints(true, false);
  NURBSCurve curve3D(controlPoints3D, 2);

  Envelope envelope3D = curve3D.envelope();
  BOOST_CHECK(envelope3D.is3D());

  for (const auto &point : controlPoints3D) {
    BOOST_CHECK(envelope3D.zMin() <= point.z());
    BOOST_CHECK(envelope3D.zMax() >= point.z());
  }
}

BOOST_AUTO_TEST_CASE(testBoundaryAlgorithm)
{
  auto       controlPoints = createTestPoints();
  NURBSCurve openCurve(controlPoints, 3);

  // Use BoundaryVisitor directly following SFCGAL patterns
  algorithm::BoundaryVisitor visitor;
  openCurve.accept(visitor);
  std::unique_ptr<Geometry> boundary(visitor.releaseBoundary());

  BOOST_REQUIRE(boundary != nullptr);
  BOOST_CHECK(boundary->is<MultiPoint>());

  const auto &multiPoint = boundary->as<MultiPoint>();
  BOOST_CHECK_EQUAL(multiPoint.numGeometries(), 2U); // Start and end points

  // Check boundary points are curve endpoints
  auto  bounds     = openCurve.parameterBounds();
  Point curveStart = openCurve.evaluate(bounds.first);
  Point curveEnd   = openCurve.evaluate(bounds.second);

  BOOST_CHECK(
      isNearlyEqual(multiPoint.geometryN(0).as<Point>(), curveStart, 1e-8) ||
      isNearlyEqual(multiPoint.geometryN(0).as<Point>(), curveEnd, 1e-8));
  BOOST_CHECK(
      isNearlyEqual(multiPoint.geometryN(1).as<Point>(), curveStart, 1e-8) ||
      isNearlyEqual(multiPoint.geometryN(1).as<Point>(), curveEnd, 1e-8));

  // Test closed curve (should have empty boundary)
  std::vector<Point> closedPoints = controlPoints;
  closedPoints.push_back(controlPoints[0]); // Close the curve
  NURBSCurve closedCurve(closedPoints, 3);

  algorithm::BoundaryVisitor closedVisitor;
  closedCurve.accept(closedVisitor);
  std::unique_ptr<Geometry> closedBoundary(closedVisitor.releaseBoundary());
  BOOST_REQUIRE(closedBoundary != nullptr);
  BOOST_CHECK(closedBoundary->isEmpty() ||
              (closedBoundary->is<MultiPoint>() &&
               closedBoundary->as<MultiPoint>().numGeometries() == 0));
}

BOOST_AUTO_TEST_CASE(testIsValidAlgorithm)
{
  // Valid curve
  auto       controlPoints = createTestPoints();
  auto       weights       = convertWeights({1.0, 2.0, 1.5, 1.0});
  NURBSCurve validCurve(controlPoints, weights, 2);

  auto validity = algorithm::isValid(validCurve);
  BOOST_CHECK(validity.valid());
  BOOST_CHECK(validity.reason().empty());

  // Test with different tolerance
  auto validityTol = algorithm::isValid(validCurve, 1e-6);
  BOOST_CHECK(validityTol.valid());

  // Empty curve should be valid
  NURBSCurve emptyCurve;
  auto       emptyValidity = algorithm::isValid(emptyCurve);
  BOOST_CHECK(emptyValidity.valid());
}

// TODO
// BOOST_AUTO_TEST_CASE(testDistanceAlgorithm)
// {
//   auto       controlPoints1 = createTestPoints();
//   NURBSCurve curve1(controlPoints1, 2);
//
//   // Create second curve offset in Y direction
//   std::vector<Point> controlPoints2;
//   for (const auto &point : controlPoints1) {
//     controlPoints2.emplace_back(point.x(), point.y() + 10.0);
//   }
//   NURBSCurve curve2(controlPoints2, 2);
//
//   // Distance should be approximately 10 (minimum separation)
//   NURBSCurve::FT distance = algorithm::distance(curve1, curve2);
//   std::cout << "DISTANCE : " << distance << "\n";
//   BOOST_CHECK(distance > NURBSCurve::FT(9.0));
//   BOOST_CHECK(distance < NURBSCurve::FT(11.0));
//
//   // Distance to point
//   Point          testPoint(100.0, 100.0);
//   NURBSCurve::FT pointDistance = algorithm::distance(curve1, testPoint);
//   std::cout << "DISTANCE : " << distance << "\n";
//   BOOST_CHECK(pointDistance > NURBSCurve::FT(90.0)); // Far away
//
//   // Distance to self should be 0
//   NURBSCurve::FT selfDistance = algorithm::distance(curve1, curve1);
//   BOOST_CHECK_SMALL(CGAL::to_double(selfDistance), 1e-10);
// }

//-- IO integration tests

BOOST_AUTO_TEST_CASE(testWKTIntegrationRoundTrip)
{
  // Test various coordinate types through WKT roundtrip
  std::vector<std::tuple<bool, bool, std::string>> testCases = {
      {false, false, "2D"},
      {true, false, "3D"},
      {false, true, "2DM"},
      {true, true, "3DM"}};

  for (const auto &[is3D, isMeasured, desc] : testCases) {
    auto originalPoints = createTestPoints(is3D, isMeasured);
    auto weights        = convertWeights({1.0, 2.5, 0.8, 1.2});

    NURBSCurve originalCurve(originalPoints, weights, 2);

    // Write to WKT
    std::string wktText = originalCurve.asText(6);
    BOOST_CHECK(!wktText.empty());

    // Read back
    auto readGeometry = io::readWkt(wktText);
    BOOST_REQUIRE(readGeometry->is<NURBSCurve>());

    const auto &readCurve = readGeometry->as<NURBSCurve>();

    // Check properties preserved
    BOOST_CHECK_EQUAL(originalCurve.is3D(), readCurve.is3D());
    BOOST_CHECK_EQUAL(originalCurve.isMeasured(), readCurve.isMeasured());
    BOOST_CHECK_EQUAL(originalCurve.isRational(), readCurve.isRational());
    BOOST_CHECK_EQUAL(originalCurve.degree(), readCurve.degree());
    BOOST_CHECK_EQUAL(originalCurve.numControlPoints(),
                      readCurve.numControlPoints());

    // Check curve evaluation equivalence
    auto                  bounds = originalCurve.parameterBounds();
    NURBSCurve::Parameter testParam =
        bounds.first + (bounds.second - bounds.first) * NURBSCurve::FT(0.5);

    Point originalPoint = originalCurve.evaluate(testParam);
    Point readPoint     = readCurve.evaluate(testParam);

    BOOST_CHECK(isNearlyEqual(originalPoint, readPoint, 1e-6));
  }
}

BOOST_AUTO_TEST_CASE(testComplexGeometryCollectionWithNURBS)
{
  // Create complex geometry collection with NURBS and other types
  auto       nurbsPoints = createTestPoints();
  NURBSCurve nurbsCurve(nurbsPoints, 3);

  GeometryCollection collection;
  collection.addGeometry(nurbsCurve.clone());
  collection.addGeometry(new Point(0.0, 0.0));

  std::vector<Point> linePoints = {Point(0.0, 0.0), Point(5.0, 5.0)};
  LineString         lineString;
  for (const auto &point : linePoints) {
    lineString.addPoint(point);
  }
  collection.addGeometry(lineString.clone());

  // Test WKT roundtrip of collection
  std::string collectionWkt = collection.asText(3);
  BOOST_CHECK(collectionWkt.find("NURBSCURVE") != std::string::npos);
  BOOST_CHECK(collectionWkt.find("POINT") != std::string::npos);
  BOOST_CHECK(collectionWkt.find("LINESTRING") != std::string::npos);

  auto readCollection = io::readWkt(collectionWkt);
  BOOST_REQUIRE(readCollection->is<GeometryCollection>());

  const auto &readGC = readCollection->as<GeometryCollection>();
  BOOST_CHECK_EQUAL(readGC.numGeometries(), 3U);
  BOOST_CHECK(readGC.geometryN(0).is<NURBSCurve>());
  BOOST_CHECK(readGC.geometryN(1).is<Point>());
  BOOST_CHECK(readGC.geometryN(2).is<LineString>());
}

//-- Multi-dimensional coordinate handling

BOOST_AUTO_TEST_CASE(testMixedDimensionHandling)
{
  // Test error handling when mixing dimensions
  std::vector<Point> mixedPoints;
  mixedPoints.emplace_back(0.0, 0.0);      // 2D
  mixedPoints.emplace_back(1.0, 1.0, 1.0); // 3D - should cause error

  // Constructor should detect inconsistent dimensions
  BOOST_CHECK_THROW(NURBSCurve curve(mixedPoints, 1), Exception);

  // Test setting inconsistent control point
  auto       validPoints = createTestPoints();
  NURBSCurve validCurve(validPoints, 2);

  Point inconsistent3D(10.0, 10.0, 10.0);
  BOOST_CHECK_THROW(validCurve.setControlPoint(0, inconsistent3D), Exception);
}

BOOST_AUTO_TEST_CASE(testMeasuredCoordinatePreservation)
{
  auto       measuredPoints = createTestPoints(false, true);
  NURBSCurve measuredCurve(measuredPoints, 3);

  BOOST_CHECK(measuredCurve.isMeasured());

  // Test that M coordinates are preserved through operations
  auto  bounds    = measuredCurve.parameterBounds();
  Point evalPoint = measuredCurve.evaluate(bounds.first);

  BOOST_CHECK(evalPoint.isMeasured());
  BOOST_CHECK(!std::isnan(evalPoint.m()));

  // Test derivatives preserve M coordinate structure
  Point derivative = measuredCurve.derivative(bounds.first, 1);
  BOOST_CHECK(derivative.isMeasured());

  // Test transformation preserves M coordinates
  CGAL::Aff_transformation_2<Kernel> identityMatrix(CGAL::IDENTITY);
  transform::AffineTransform2        identity(identityMatrix);

  NURBSCurve transformedCurve = measuredCurve;
  transformedCurve.accept(identity);

  BOOST_CHECK(transformedCurve.isMeasured());
  for (size_t idx = 0; idx < transformedCurve.numControlPoints(); ++idx) {
    BOOST_CHECK(transformedCurve.controlPointN(idx).isMeasured());
    BOOST_CHECK_EQUAL(transformedCurve.controlPointN(idx).m(),
                      measuredCurve.controlPointN(idx).m());
  }
}

BOOST_AUTO_TEST_CASE(test3DZMCoordinateHandling)
{
  auto       points3DM = createTestPoints(true, true);
  auto       weights   = convertWeights({1.0, 1.5, 2.0, 0.5});
  NURBSCurve curve3DM(points3DM, weights, 3);

  BOOST_CHECK(curve3DM.is3D());
  BOOST_CHECK(curve3DM.isMeasured());
  BOOST_CHECK(curve3DM.isRational());
  BOOST_CHECK_EQUAL(curve3DM.coordinateDimension(),
                    3); // XYZ, M handled separately

  // Test evaluation preserves all coordinate types
  auto  bounds = curve3DM.parameterBounds();
  Point evalPoint =
      curve3DM.evaluate((bounds.first + bounds.second) / NURBSCurve::FT(2));

  BOOST_CHECK(evalPoint.is3D());
  BOOST_CHECK(evalPoint.isMeasured());
  BOOST_CHECK(CGAL::is_finite(evalPoint.x()));
  BOOST_CHECK(CGAL::is_finite(evalPoint.y()));
  BOOST_CHECK(CGAL::is_finite(evalPoint.z()));
  BOOST_CHECK(!std::isnan(evalPoint.m()));

  // Test coordinate dropping
  curve3DM.dropZ();
  BOOST_CHECK(!curve3DM.is3D());
  BOOST_CHECK(curve3DM.isMeasured());

  curve3DM.dropM();
  BOOST_CHECK(!curve3DM.is3D());
  BOOST_CHECK(!curve3DM.isMeasured());
}

//-- Performance and stress integration tests

BOOST_AUTO_TEST_CASE(testLargeNURBSCurveIntegration)
{
  // Create large NURBS curve
  std::vector<Point>  manyPoints;
  std::vector<double> manyWeights;

  for (int idx = 0; idx < 200; ++idx) {
    double angle  = 4.0 * M_PI * idx / 199.0; // Two full spirals
    double radius = 1.0 + 0.5 * idx / 199.0;  // Expanding spiral
    manyPoints.emplace_back(radius * std::cos(angle), radius * std::sin(angle),
                            0.1 * idx); // 3D spiral
    manyWeights.push_back(1.0 + 0.2 * std::sin(6.0 * angle));
  }

  auto       weights = convertWeights(manyWeights);
  NURBSCurve largeCurve(manyPoints, weights, 5);

  // All operations should handle large curve without issues
  BOOST_CHECK_NO_THROW(auto validity = algorithm::isValid(largeCurve));
  BOOST_CHECK_NO_THROW(auto envelope = largeCurve.envelope());
  BOOST_CHECK_NO_THROW(std::string wktText = largeCurve.asText(2));

  // Test transformation
  CGAL::Aff_transformation_3<Kernel> scalingMatrix(
      2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0);
  transform::AffineTransform3 scaling(scalingMatrix);

  BOOST_CHECK_NO_THROW(largeCurve.accept(scaling));

  // Test evaluation at multiple points
  auto bounds = largeCurve.parameterBounds();
  for (int testIdx = 0; testIdx <= 20; ++testIdx) {
    NURBSCurve::Parameter param =
        bounds.first + (bounds.second - bounds.first) *
                           NURBSCurve::FT(testIdx) / NURBSCurve::FT(20);

    BOOST_CHECK_NO_THROW(Point point = largeCurve.evaluate(param));
    Point point = largeCurve.evaluate(param);
    BOOST_CHECK(CGAL::is_finite(point.x()));
    BOOST_CHECK(CGAL::is_finite(point.y()));
    BOOST_CHECK(CGAL::is_finite(point.z()));
  }
}

BOOST_AUTO_TEST_CASE(testNURBSCurveConversionCompatibility)
{
  auto       controlPoints = createTestPoints();
  auto       weights       = convertWeights({1.0, 3.0, 1.5, 1.0});
  NURBSCurve originalCurve(controlPoints, weights, 3);

  // Test LineString conversion compatibility
  auto lineString = originalCurve.toLineString(50);
  BOOST_REQUIRE(lineString != nullptr);
  BOOST_CHECK_EQUAL(lineString->numPoints(), 51U);

  // LineString should preserve coordinate dimensions
  BOOST_CHECK_EQUAL(lineString->is3D(), originalCurve.is3D());
  BOOST_CHECK_EQUAL(lineString->isMeasured(), originalCurve.isMeasured());

  // Test that LineString can be used in algorithms
  auto lineEnvelope  = lineString->envelope();
  auto curveEnvelope = originalCurve.envelope();

  // Envelopes should be similar (LineString is approximation)
  BOOST_CHECK(std::abs(CGAL::to_double(lineEnvelope.xMin() -
                                       curveEnvelope.xMin())) < 1.0);
  BOOST_CHECK(std::abs(CGAL::to_double(lineEnvelope.yMin() -
                                       curveEnvelope.yMin())) < 1.0);

  // Test adaptive LineString conversion
  auto adaptiveLineString =
      originalCurve.toLineStringAdaptive(NURBSCurve::FT(0.1), 10, 100);
  BOOST_REQUIRE(adaptiveLineString != nullptr);
  BOOST_CHECK(adaptiveLineString->numPoints() >=
              11U); // At least minSegments + 1
  BOOST_CHECK(adaptiveLineString->numPoints() <=
              101U); // At most maxSegments + 1
}

//-- Error propagation through system

BOOST_AUTO_TEST_CASE(testErrorPropagationIntegration)
{
  NURBSCurve emptyCurve;

  // Operations on empty curves should propagate errors appropriately
  BOOST_CHECK_NO_THROW(auto envelope = emptyCurve.envelope());
  auto envelope = emptyCurve.envelope();
  BOOST_CHECK(envelope.isEmpty());

  algorithm::BoundaryVisitor visitor;
  BOOST_CHECK_NO_THROW(emptyCurve.accept(visitor));
  std::unique_ptr<Geometry> boundary(visitor.releaseBoundary());
  BOOST_CHECK(boundary->isEmpty());

  BOOST_CHECK_NO_THROW(auto validity = algorithm::isValid(emptyCurve));
  auto validity = algorithm::isValid(emptyCurve);
  BOOST_CHECK(validity.valid()); // Empty should be valid

  // Test transform on empty curve
  CGAL::Aff_transformation_3<Kernel> someTransformMatrix(
      2.0, 0.0, 0.0, 5.0, 0.0, 2.0, 0.0, 3.0, 0.0, 0.0, 2.0, 1.0);
  transform::AffineTransform3 someTransform(someTransformMatrix);

  BOOST_CHECK_NO_THROW(emptyCurve.accept(someTransform));
  BOOST_CHECK(emptyCurve.isEmpty()); // Should remain empty
}

BOOST_AUTO_TEST_SUITE_END()
