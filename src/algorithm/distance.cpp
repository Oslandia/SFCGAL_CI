// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/distance.h"

#include "SFCGAL/Curve.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/algorithm/isValid.h"

#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_with_holes_2.h>

#include "SFCGAL/algorithm/intersects.h"
#include "SFCGAL/detail/GetPointsVisitor.h"
#include "SFCGAL/detail/transform/AffineTransform3.h"
#include <algorithm>
#include <limits>

namespace SFCGAL::algorithm {

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

/**
 * @brief Compute 2D distance between two geometries
 */
auto
distance(const Geometry &geometry1, const Geometry &geometry2) -> double
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(geometry1);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(geometry2);
  return distance(geometry1, geometry2, NoValidityCheck());
}

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection
#ifndef DOXYGEN_SHOULD_SKIP_THIS

auto
distance(const Geometry &geometry1, const Geometry &geometry2,
         [[maybe_unused]] NoValidityCheck noCheck) -> double
{
  switch (geometry1.geometryTypeId()) {
  case TYPE_POINT:
    return distancePointGeometry(geometry1.as<Point>(), geometry2);

  case TYPE_LINESTRING:
    return distanceLineStringGeometry(geometry1.as<LineString>(), geometry2);

  case TYPE_NURBSCURVE: {
    auto lineString =
        geometry1.as<NURBSCurve>().toLineString(); // default parameters
    if (!lineString || lineString->isEmpty()) {
      return std::numeric_limits<double>::infinity();
    }
    return distanceLineStringGeometry(*lineString, geometry2);
  }

  case TYPE_POLYGON:
    return distancePolygonGeometry(geometry1.as<Polygon>(), geometry2);

  case TYPE_TRIANGLE:
    return distanceTriangleGeometry(geometry1.as<Triangle>(), geometry2);

  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
    return distanceGeometryCollectionToGeometry(geometry1, geometry2);

  case TYPE_SOLID:
    BOOST_THROW_EXCEPTION(NotImplementedException(
        (boost::format("distance(%s,%s) is not implemented") %
         geometry1.geometryType() % geometry2.geometryType())
            .str()));
  }

  BOOST_ASSERT(false);
  return 0;
}

auto
distancePointGeometry(const Point &point, const Geometry &geometry) -> double
{
  switch (geometry.geometryTypeId()) {
  case TYPE_POINT:
    return distancePointPoint(point, geometry.as<Point>());

  case TYPE_LINESTRING:
    return distancePointLineString(point, geometry.as<LineString>());

  case TYPE_NURBSCURVE: {
    auto lineString =
        geometry.as<NURBSCurve>().toLineString(); // default parameters
    if (!lineString || lineString->isEmpty()) {
      return std::numeric_limits<double>::infinity();
    }
    return distancePointLineString(point, *lineString);
  }

  case TYPE_POLYGON:
    return distancePointPolygon(point, geometry.as<Polygon>());

  case TYPE_TRIANGLE:
    return distancePointTriangle(point, geometry.as<Triangle>());

    // collection dispatch
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
    return distanceGeometryCollectionToGeometry(geometry, point);

  case TYPE_SOLID:
    BOOST_THROW_EXCEPTION(NotImplementedException(
        (boost::format("distance(%s,%s) is not implemented") %
         point.geometryType() % geometry.geometryType())
            .str()));
  }

  BOOST_ASSERT(false);
  return 0;
}

auto
distancePointPoint(const Point &point1, const Point &point2) -> double
{
  if (point1.isEmpty() || point2.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  return CGAL::sqrt(CGAL::to_double(
      CGAL::squared_distance(point1.toPoint_2(), point2.toPoint_2())));
}

auto
distancePointLineString(const Point &point, const LineString &lineString)
    -> double
{
  if (point.isEmpty() || lineString.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  size_t const numSegments = lineString.numSegments();

  double dMin = std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < numSegments; i++) {
    double const segmentDistance = distancePointSegment(
        point, lineString.pointN(i), lineString.pointN(i + 1));

    if (i == 0 || segmentDistance < dMin) {
      dMin = segmentDistance;
    }
  }

  return dMin;
}

auto
distancePointPolygon(const Point &point, const Polygon &polygon) -> double
{
  if (point.isEmpty() || polygon.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  if (intersects(point, polygon, NoValidityCheck())) {
    return 0.0;
  }

  double dMin = 0.0;

  // check if the point is in the polygon
  for (size_t i = 0; i < polygon.numRings(); i++) {
    double const ringDistance =
        distancePointLineString(point, polygon.ringN(i));

    if (i == 0 || ringDistance < dMin) {
      dMin = ringDistance;
    }
  }

  return dMin;
}

auto
distancePointTriangle(const Point &point, const Triangle &triangle) -> double
{
  return distancePointPolygon(point, triangle.toPolygon());
}

auto
distanceLineStringGeometry(const LineString &lineString,
                           const Geometry   &geometry) -> double
{
  if (lineString.isEmpty() || geometry.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  switch (geometry.geometryTypeId()) {
  case TYPE_POINT:
    return distancePointLineString(geometry.as<Point>(),
                                   lineString); // symetric

  case TYPE_LINESTRING:
    return distanceLineStringLineString(lineString, geometry.as<LineString>());

  case TYPE_NURBSCURVE: {
    auto lineStringFromCurve =
        geometry.as<NURBSCurve>().toLineString(); // default parameters
    if (!lineStringFromCurve || lineStringFromCurve->isEmpty()) {
      return std::numeric_limits<double>::infinity();
    }
    return distanceLineStringLineString(lineString, *lineStringFromCurve);
  }

  case TYPE_POLYGON:
    return distanceLineStringPolygon(lineString, geometry.as<Polygon>());

  case TYPE_TRIANGLE:
    return distanceLineStringTriangle(lineString, geometry.as<Triangle>());

    // collection dispatch
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
    return distanceGeometryCollectionToGeometry(geometry, lineString);

  case TYPE_SOLID:
    BOOST_THROW_EXCEPTION(NotImplementedException(
        (boost::format("distance(%s,%s) is not implemented") %
         lineString.geometryType() % geometry.geometryType())
            .str()));
  }

  BOOST_ASSERT(false);
  return 0;
}

auto
distanceLineStringLineString(const LineString &lineString1,
                             const LineString &lineString2) -> double
{
  if (lineString1.isEmpty() || lineString2.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  size_t const ns1 = lineString1.numSegments();
  size_t const ns2 = lineString2.numSegments();

  double dMin = std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < ns1; i++) {
    for (size_t j = 0; j < ns2; j++) {
      dMin = std::min(dMin, distanceSegmentSegment(lineString1.pointN(i),
                                                   lineString1.pointN(i + 1),
                                                   lineString2.pointN(j),
                                                   lineString2.pointN(j + 1)));
    }
  }

  return dMin;
}

auto
distanceLineStringPolygon(const LineString &lineString, const Polygon &polygon)
    -> double
{
  if (lineString.isEmpty() || polygon.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  if (intersects(lineString, polygon, NoValidityCheck())) {
    return 0.0;
  }

  double dMin = std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < polygon.numRings(); i++) {
    double const lineStringDistance =
        distanceLineStringLineString(lineString, polygon.ringN(i));

    dMin = std::min(lineStringDistance, dMin);
  }

  return dMin;
}

auto
distanceLineStringTriangle(const LineString &lineString,
                           const Triangle   &triangle) -> double
{
  return distanceLineStringPolygon(lineString, triangle.toPolygon());
}

auto
distancePolygonGeometry(const Polygon &polygon, const Geometry &geometry)
    -> double
{
  if (polygon.isEmpty() || geometry.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  switch (geometry.geometryTypeId()) {
  case TYPE_POINT:
    return distancePointPolygon(geometry.as<Point>(), polygon); // symetric

  case TYPE_LINESTRING:
    return distanceLineStringPolygon(geometry.as<LineString>(),
                                     polygon); // symetric

  case TYPE_NURBSCURVE: {
    auto lineString =
        geometry.as<NURBSCurve>().toLineString(); // default parameters
    if (!lineString || lineString->isEmpty()) {
      return std::numeric_limits<double>::infinity();
    }
    return distanceLineStringPolygon(*lineString, polygon); // symetric
  }

  case TYPE_POLYGON:
    return distancePolygonPolygon(polygon, geometry.as<Polygon>());

  case TYPE_TRIANGLE:
    return distancePolygonTriangle(polygon, geometry.as<Triangle>());

    // collection dispatch
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
    return distanceGeometryCollectionToGeometry(geometry, polygon);

  case TYPE_SOLID:
    BOOST_THROW_EXCEPTION(NotImplementedException(
        (boost::format("distance(%s,%s) is not implemented") %
         polygon.geometryType() % geometry.geometryType())
            .str()));
  }

  BOOST_ASSERT(false);
  return 0;
}

auto
distancePolygonPolygon(const Polygon &polygon1, const Polygon &polygon2)
    -> double
{
  if (polygon1.isEmpty() || polygon2.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  if (intersects(polygon1, polygon2, NoValidityCheck())) {
    return 0.0;
  }

  double dMin = std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < polygon1.numRings(); i++) {
    for (size_t j = 0; j < polygon2.numRings(); j++) {
      double const d =
          distanceLineStringLineString(polygon1.ringN(i), polygon2.ringN(j));

      dMin = std::min(d, dMin);
    }
  }

  return dMin;
}

auto
distancePolygonTriangle(const Polygon &polygon, const Triangle &triangle)
    -> double
{
  return distancePolygonPolygon(polygon, triangle.toPolygon());
}

auto
distanceTriangleGeometry(const Triangle &triangle, const Geometry &geometry)
    -> double
{
  return distancePolygonGeometry(triangle.toPolygon(), geometry);
}

/**
 * @struct Circle
 * @brief Represents a 2D circle with a center and a radius.
 */
struct Circle {
  /**
   * @brief Main constructor.
   *
   */
  Circle(const double &radius, CGAL::Vector_2<Kernel> &center)
      : _radius(radius), _center(center), _empty(false)
  {
  }

  /**
   * @brief Default constructor.
   *
   * Initializes an empty circle.
   */
  Circle() = default;

  /**
   * @brief Checks if the circle is empty.
   * @return true if the circle is empty, false otherwise.
   */
  [[nodiscard]] auto
  isEmpty() const -> bool
  {
    return _empty;
  }

  /**
   * @brief Returns the radius of the circle.
   * @return The radius of the circle.
   *
   * @note Assertion: the circle must not be empty.
   */
  [[nodiscard]] auto
  radius() const -> double
  {
    BOOST_ASSERT(!_empty);
    return _radius;
  }

  /**
   * @brief Returns the center of the circle.
   * @return A constant reference to the circle's center.
   *
   * @note Assertion: the circle must not be empty.
   */
  [[nodiscard]] auto
  center() const -> const CGAL::Vector_2<Kernel> &
  {
    BOOST_ASSERT(!_empty);
    return _center;
  }

private:
  double                 _radius{}; ///< Radius of the circle
  CGAL::Vector_2<Kernel> _center;   ///< Center of the circle
  bool _empty{true};                ///< Indicates whether the circle is empty
};

/// \cond IGNORE
auto
boundingCircle(const Geometry &geom) -> const Circle
{
  if (geom.isEmpty()) {
    return Circle();
  }

  using namespace SFCGAL::detail;
  GetPointsVisitor v;
  const_cast<Geometry &>(geom).accept(v);

  if (v.points.empty()) {
    return Circle();
  }

  const auto end = v.points.end();

  // centroid
  Vector_2 c(0, 0);
  int      numPoint = 0;

  for (auto x = v.points.begin(); x != end; ++x) {
    c = c + (*x)->toVector_2();
    ++numPoint;
  }

  BOOST_ASSERT(numPoint);
  c = c / numPoint;

  // farest point from centroid
  Vector_2   f             = c;
  Kernel::FT maxDistanceSq = 0;

  for (auto x = v.points.begin(); x != end; ++x) {
    const Vector_2   cx  = (*x)->toVector_2() - c;
    const Kernel::FT dSq = cx * cx;

    if (dSq > maxDistanceSq) {
      f             = (*x)->toVector_2();
      maxDistanceSq = dSq;
    }
  }

  return Circle(std::sqrt(CGAL::to_double(maxDistanceSq)), c);
}
/// \endcond

auto
distanceGeometryCollectionToGeometry(const Geometry &geometryCollection,
                                     const Geometry &geometry) -> double
{
  if (geometryCollection.isEmpty() || geometry.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  // if bounding spheres (BS) of geometry and geometryCollectioni don't
  // intersect and if the closest point of BS(geometryCollectionj) is further
  // than the farest point of BS(geometryCollectioni) there is no need to
  // compute the distance(geometryCollectionj, geometry) since it will be
  // greater than distance(geometryCollectioni, geometry)
  //
  // The aim is not to find the minimal bounding sphere, but a good enough
  // sphere than encloses all points
  std::set<size_t> noTest;

  if (true) {
    std::vector<Circle> bcCollection;

    bcCollection.reserve(geometryCollection.numGeometries());
    for (size_t i = 0; i < geometryCollection.numGeometries(); i++) {
      bcCollection.push_back(boundingCircle(geometryCollection.geometryN(i)));
    }

    Circle const bcGeometry(boundingCircle(geometry));

    if (bcGeometry.isEmpty()) {
      return std::numeric_limits<double>::infinity();
    }

    std::vector<size_t> noIntersect;

    for (size_t i = 0; i < geometryCollection.numGeometries(); i++) {
      if (bcCollection[i].isEmpty()) {
        continue;
      }

      const double l2 = CGAL::to_double(
          (bcGeometry.center() - bcCollection[i].center()).squared_length());

      if (std::pow(bcGeometry.radius() + bcCollection[i].radius(), 2) < l2) {
        noIntersect.push_back(i);
      }
    }

    for (size_t i = 0; i < noIntersect.size(); i++) {
      const double li = std::sqrt(CGAL::to_double(
          (bcCollection[noIntersect[i]].center() - bcGeometry.center())
              .squared_length()));

      for (size_t j = i; j < noIntersect.size(); j++) {
        const double lj = std::sqrt(CGAL::to_double(
            (bcCollection[noIntersect[j]].center() - bcGeometry.center())
                .squared_length()));

        if (li + bcCollection[noIntersect[i]].radius() <
            lj - bcCollection[noIntersect[j]].radius()) {
          noTest.insert(noIntersect[j]);
        } else if (lj + bcCollection[noIntersect[j]].radius() <
                   li - bcCollection[noIntersect[i]].radius()) {
          noTest.insert(noIntersect[i]);
        }
      }
    }

    // if (!noTest.empty()) std::cout << "pruning " << noTest.size() << "/" <<
    // geometryCollection.numGeometries() << "\n";
  }

  double dMin = std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < geometryCollection.numGeometries(); i++) {
    if (noTest.end() != noTest.find(i)) {
      continue;
    }

    dMin = std::min(dMin, distance(geometryCollection.geometryN(i), geometry));
  }

  return dMin;
}

auto
distancePointSegment(const Point &point, const Point &segmentStart,
                     const Point &segmentEnd) -> double
{
  // empty already checked
  BOOST_ASSERT(!point.isEmpty());
  BOOST_ASSERT(!segmentStart.isEmpty());
  BOOST_ASSERT(!segmentEnd.isEmpty());

  return CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(
      point.toPoint_2(),
      Segment_2(segmentStart.toPoint_2(), segmentEnd.toPoint_2()))));
}

auto
distanceSegmentSegment(const Point &point1, const Point &point2,
                       const Point &point3, const Point &point4) -> double
{
  // empty already checked
  BOOST_ASSERT(!point1.isEmpty());
  BOOST_ASSERT(!point2.isEmpty());
  BOOST_ASSERT(!point3.isEmpty());
  BOOST_ASSERT(!point4.isEmpty());

  return CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(
      CGAL::Segment_2<Kernel>(point1.toPoint_2(), point2.toPoint_2()),
      CGAL::Segment_2<Kernel>(point3.toPoint_2(), point4.toPoint_2()))));
}

#endif // ifndef DOXYGEN_SHOULD_SKIP_THIS
/// @} end of private section

} // namespace SFCGAL::algorithm
