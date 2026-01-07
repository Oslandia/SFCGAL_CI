// Copyright (c) 2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/insertPointsWithinTolerance.h"

#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/algorithm/distance.h"
#include "SFCGAL/detail/GetPointsVisitor.h"
#include "SFCGAL/numeric.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

namespace SFCGAL::algorithm {

/**
 * @brief Extract all points from a geometry using GetPointsVisitor
 * @param geometry The input geometry to extract points from
 * @param points The output vector to store extracted points
 */
void
extractPoints(const Geometry &geometry, std::vector<Point> &points)
{
  SFCGAL::detail::GetPointsVisitor visitor;
  geometry.accept(visitor);

  // Copy the collected points to the output vector
  for (const Point *point : visitor.points) {
    points.push_back(*point);
  }
}

/**
 * @brief Find the insertion position for a point on a linestring
 * @param lineString The linestring to find insertion position on
 * @param point The point to find insertion position for
 * @return The index where the point should be inserted
 */
auto
findInsertionPosition(const LineString &lineString, const Point &point)
    -> size_t
{
  if (lineString.numPoints() < 2) {
    return 0; // If there's only one point or none, insert at the beginning
  }

  // Find the segment closest to the point
  double minDistance = std::numeric_limits<double>::max();
  size_t bestPos     = 0;

  for (size_t i = 0; i < lineString.numSegments(); ++i) {
    const Point &segStart = lineString.pointN(i);
    const Point &segEnd   = lineString.pointN(i + 1);

    double dist = distancePointSegment(point, segStart, segEnd);

    if (dist < minDistance) {
      minDistance = dist;
      bestPos     = i + 1; // Insert after the start of the segment
    }
  }

  return bestPos;
}

/**
 * @brief Insert a point into a linestring at a specific position
 * @param lineString The original linestring
 * @param point The point to insert
 * @param position The position to insert the point at
 * @return A new linestring with the point inserted
 */
auto
insertPointIntoLineString(const LineString &lineString, const Point &point,
                          size_t position) -> std::unique_ptr<LineString>
{
  auto newLineString = std::make_unique<LineString>();

  // Add points up to the insertion position
  for (size_t i = 0; i < position && i < lineString.numPoints(); ++i) {
    newLineString->addPoint(lineString.pointN(i));
  }

  // Add the new point
  newLineString->addPoint(point);

  // Add remaining points
  for (size_t i = position; i < lineString.numPoints(); ++i) {
    newLineString->addPoint(lineString.pointN(i));
  }

  return newLineString;
}

/**
 * @brief Insert points into a linestring based on distance tolerance
 * @param lineString The original linestring to insert points into
 * @param points The points to insert
 * @param tolerance The maximum distance for a point to be considered for
 * insertion
 * @return A new linestring with the points inserted
 */
auto
insertPointsIntoLineString(const LineString         &lineString,
                           const std::vector<Point> &points, double tolerance)
    -> std::unique_ptr<LineString>
{
  // Start with the original lineString
  auto result = lineString.clone();

  // Create a list of points to insert with their distances to the closest
  // segments
  std::vector<std::pair<Point, double>> pointsToInsert;

  for (const Point &point : points) {
    // Find the closest segment and its distance
    double minDist = std::numeric_limits<double>::max();

    for (size_t i = 0; i < result->numSegments(); ++i) {
      const Point &segStart = result->pointN(i);
      const Point &segEnd   = result->pointN(i + 1);

      double dist = distancePointSegment(point, segStart, segEnd);

      minDist = std::min(minDist, dist);
    }

    // If the distance is within tolerance, mark for insertion
    if (minDist <= tolerance) {
      pointsToInsert.emplace_back(point, minDist);
    }
  }

  // Sort points by distance to ensure closer points are inserted first
  std::sort(pointsToInsert.begin(), pointsToInsert.end(),
            [](const std::pair<Point, double> &a,
               const std::pair<Point, double> &b) -> bool {
              return a.second < b.second;
            });

  // Insert points one by one
  for (const auto &pointPair : pointsToInsert) {
    const Point &point = pointPair.first;

    // Recalculate the best position after previous insertions
    size_t pos = findInsertionPosition(*result, point);

    // Make sure the position is valid
    if (pos <= result->numPoints()) {
      auto temp = insertPointIntoLineString(*result, point, pos);
      result    = std::move(temp);
    }
  }

  return result;
}

namespace {
/**
 * @brief Process a single polygon by inserting points into its rings
 * @param polygon The polygon to process
 * @param sourcePoints The points to insert
 * @param tolerance The maximum distance for point insertion
 * @return A new polygon with points inserted
 */
auto
processPolygon(const Polygon &polygon, const std::vector<Point> &sourcePoints,
               double tolerance) -> std::unique_ptr<Polygon>
{
  auto resultPolygon = std::make_unique<Polygon>();

  // Process exterior ring
  auto exteriorRing = insertPointsIntoLineString(polygon.exteriorRing(),
                                                 sourcePoints, tolerance);
  resultPolygon->setExteriorRing(std::move(exteriorRing));

  // Process interior rings
  for (size_t i = 0; i < polygon.numInteriorRings(); ++i) {
    auto newHole = insertPointsIntoLineString(polygon.interiorRingN(i),
                                              sourcePoints, tolerance);
    resultPolygon->addInteriorRing(std::move(newHole));
  }

  return resultPolygon;
}
} // namespace

auto
insertPointsWithinTolerance(const Geometry &baseGeometry,
                            const Geometry &sourceGeometry, double tolerance)
    -> std::unique_ptr<Geometry>
{
  // If the base or source geometry is empty, return base geometry directly
  if ((baseGeometry.isEmpty()) || (sourceGeometry.isEmpty())) {
    return baseGeometry.clone();
  }

  // Extract all points from the source geometry
  std::vector<Point> sourcePoints;
  extractPoints(sourceGeometry, sourcePoints);

  // Handle different base geometry types
  switch (baseGeometry.geometryTypeId()) {
  case TYPE_LINESTRING: {
    const auto &baseLineString = baseGeometry.as<LineString>();
    auto        result =
        insertPointsIntoLineString(baseLineString, sourcePoints, tolerance);
    return std::unique_ptr<Geometry>(result.release());
  }

  case TYPE_POLYGON: {
    const auto &basePolygon = baseGeometry.as<Polygon>();
    return processPolygon(basePolygon, sourcePoints, tolerance);
  }

  case TYPE_MULTIPOINT: {
    // For MultiPoint, just return the original since we're adding points to it
    return baseGeometry.clone();
  }

  case TYPE_MULTILINESTRING: {
    const auto &baseMultiLineString   = baseGeometry.as<MultiLineString>();
    auto        resultMultiLineString = std::make_unique<MultiLineString>();

    for (size_t i = 0; i < baseMultiLineString.numGeometries(); ++i) {
      const auto &lineString =
          baseMultiLineString.geometryN(i).as<LineString>();
      auto newLineString =
          insertPointsIntoLineString(lineString, sourcePoints, tolerance);
      resultMultiLineString->addGeometry(std::move(newLineString));
    }

    return std::unique_ptr<Geometry>(resultMultiLineString.release());
  }

  case TYPE_MULTIPOLYGON: {
    const auto &baseMultiPolygon   = baseGeometry.as<MultiPolygon>();
    auto        resultMultiPolygon = std::make_unique<MultiPolygon>();

    for (size_t i = 0; i < baseMultiPolygon.numGeometries(); ++i) {
      resultMultiPolygon->addGeometry(processPolygon(
          baseMultiPolygon.polygonN(i), sourcePoints, tolerance));
    }

    return resultMultiPolygon;
  }

  default:
    // For unsupported types, return a clone of the base geometry
    return baseGeometry.clone();
  }
}

} // namespace SFCGAL::algorithm
