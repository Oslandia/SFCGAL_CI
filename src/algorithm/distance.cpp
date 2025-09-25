// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/distance.h"

#include "SFCGAL/Curve.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
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
distance(const Geometry &gA, const Geometry &gB) -> double
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(gA);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(gB);
  return distance(gA, gB, NoValidityCheck());
}

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection
#ifndef DOXYGEN_SHOULD_SKIP_THIS

auto
distance(const Geometry &gA, const Geometry &gB,
         [[maybe_unused]] NoValidityCheck noCheck) -> double
{
  switch (gA.geometryTypeId()) {
  case TYPE_POINT:
    return distancePointGeometry(gA.as<Point>(), gB);

  case TYPE_LINESTRING:
    return distanceLineStringGeometry(gA.as<LineString>(), gB);

  case TYPE_NURBSCURVE: {
    auto lineString =
        gA.as<Curve>().toLineStringAdaptive(); // default tolerance FT(1e-3)
    if (!lineString || lineString->isEmpty()) {
      lineString = gA.as<Curve>().toLineString(
          256); // fallback to denser uniform sampling
    }
    if (!lineString || lineString->isEmpty()) {
      return std::numeric_limits<double>::infinity();
    }
    return distanceLineStringGeometry(*lineString, gB);
  }

  case TYPE_POLYGON:
    return distancePolygonGeometry(gA.as<Polygon>(), gB);

  case TYPE_TRIANGLE:
    return distanceTriangleGeometry(gA.as<Triangle>(), gB);

  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
    return distanceGeometryCollectionToGeometry(gA, gB);

  case TYPE_SOLID:
    BOOST_THROW_EXCEPTION(NotImplementedException(
        (boost::format("distance(%s,%s) is not implemented") %
         gA.geometryType() % gB.geometryType())
            .str()));
  }

  BOOST_ASSERT(false);
  return 0;
}

auto
distancePointGeometry(const Point &gA, const Geometry &gB) -> double
{
  switch (gB.geometryTypeId()) {
  case TYPE_POINT:
    return distancePointPoint(gA, gB.as<Point>());

  case TYPE_LINESTRING:
    return distancePointLineString(gA, gB.as<LineString>());

  case TYPE_NURBSCURVE: {
    auto lineString =
        gB.as<Curve>().toLineStringAdaptive(); // default tolerance FT(1e-3)
    if (!lineString || lineString->isEmpty()) {
      lineString = gB.as<Curve>().toLineString(
          256); // fallback to denser uniform sampling
    }
    if (!lineString || lineString->isEmpty()) {
      return std::numeric_limits<double>::infinity();
    }
    return distancePointLineString(gA, *lineString);
  }

  case TYPE_POLYGON:
    return distancePointPolygon(gA, gB.as<Polygon>());

  case TYPE_TRIANGLE:
    return distancePointTriangle(gA, gB.as<Triangle>());

    // collection dispatch
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
    return distanceGeometryCollectionToGeometry(gB, gA);

  case TYPE_SOLID:
    BOOST_THROW_EXCEPTION(NotImplementedException(
        (boost::format("distance(%s,%s) is not implemented") %
         gA.geometryType() % gB.geometryType())
            .str()));
  }

  BOOST_ASSERT(false);
  return 0;
}

auto
distancePointPoint(const Point &gA, const Point &gB) -> double
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  return CGAL::sqrt(
      CGAL::to_double(CGAL::squared_distance(gA.toPoint_2(), gB.toPoint_2())));
}

auto
distancePointLineString(const Point &gA, const LineString &gB) -> double
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  size_t const numSegments = gB.numSegments();

  double dMin = std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < numSegments; i++) {
    double const segmentDistance =
        distancePointSegment(gA, gB.pointN(i), gB.pointN(i + 1));

    if (i == 0 || segmentDistance < dMin) {
      dMin = segmentDistance;
    }
  }

  return dMin;
}

auto
distancePointPolygon(const Point &gA, const Polygon &gB) -> double
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  if (intersects(gA, gB, NoValidityCheck())) {
    return 0.0;
  }

  double dMin = 0.0;

  // check if the point is in the polygon
  for (size_t i = 0; i < gB.numRings(); i++) {
    double const ringDistance = distancePointLineString(gA, gB.ringN(i));

    if (i == 0 || ringDistance < dMin) {
      dMin = ringDistance;
    }
  }

  return dMin;
}

auto
distancePointTriangle(const Point &gA, const Triangle &gB) -> double
{
  return distancePointPolygon(gA, gB.toPolygon());
}

auto
distanceLineStringGeometry(const LineString &gA, const Geometry &gB) -> double
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  switch (gB.geometryTypeId()) {
  case TYPE_POINT:
    return distancePointLineString(gB.as<Point>(), gA); // symetric

  case TYPE_LINESTRING:
    return distanceLineStringLineString(gA, gB.as<LineString>());

  case TYPE_NURBSCURVE: {
    auto lineString =
        gB.as<Curve>().toLineStringAdaptive(); // default tolerance FT(1e-3)
    if (!lineString || lineString->isEmpty()) {
      lineString = gB.as<Curve>().toLineString(
          256); // fallback to denser uniform sampling
    }
    if (!lineString || lineString->isEmpty()) {
      return std::numeric_limits<double>::infinity();
    }
    return distanceLineStringLineString(gA, *lineString);
  }

  case TYPE_POLYGON:
    return distanceLineStringPolygon(gA, gB.as<Polygon>());

  case TYPE_TRIANGLE:
    return distanceLineStringTriangle(gA, gB.as<Triangle>());

    // collection dispatch
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
    return distanceGeometryCollectionToGeometry(gB, gA);

  case TYPE_SOLID:
    BOOST_THROW_EXCEPTION(NotImplementedException(
        (boost::format("distance(%s,%s) is not implemented") %
         gA.geometryType() % gB.geometryType())
            .str()));
  }

  BOOST_ASSERT(false);
  return 0;
}

auto
distanceLineStringLineString(const LineString &gA, const LineString &gB)
    -> double
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  size_t const nsA = gA.numSegments();
  size_t const nsB = gB.numSegments();

  double dMin = std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < nsA; i++) {
    for (size_t j = 0; j < nsB; j++) {
      dMin = std::min(dMin,
                      distanceSegmentSegment(gA.pointN(i), gA.pointN(i + 1),
                                             gB.pointN(j), gB.pointN(j + 1)));
    }
  }

  return dMin;
}

auto
distanceLineStringPolygon(const LineString &gA, const Polygon &gB) -> double
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  if (intersects(gA, gB, NoValidityCheck())) {
    return 0.0;
  }

  double dMin = std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < gB.numRings(); i++) {
    double const lineStringDistance =
        distanceLineStringLineString(gA, gB.ringN(i));

    dMin = std::min(lineStringDistance, dMin);
  }

  return dMin;
}

auto
distanceLineStringTriangle(const LineString &gA, const Triangle &gB) -> double
{
  return distanceLineStringPolygon(gA, gB.toPolygon());
}

auto
distancePolygonGeometry(const Polygon &gA, const Geometry &gB) -> double
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  switch (gB.geometryTypeId()) {
  case TYPE_POINT:
    return distancePointPolygon(gB.as<Point>(), gA); // symetric

  case TYPE_LINESTRING:
    return distanceLineStringPolygon(gB.as<LineString>(), gA); // symetric

  case TYPE_NURBSCURVE: {
    auto lineString =
        gB.as<Curve>().toLineStringAdaptive(); // default tolerance FT(1e-3)
    if (!lineString || lineString->isEmpty()) {
      lineString = gB.as<Curve>().toLineString(
          256); // fallback to denser uniform sampling
    }
    if (!lineString || lineString->isEmpty()) {
      return std::numeric_limits<double>::infinity();
    }
    return distanceLineStringPolygon(*lineString, gA); // symetric
  }

  case TYPE_POLYGON:
    return distancePolygonPolygon(gA, gB.as<Polygon>());

  case TYPE_TRIANGLE:
    return distancePolygonTriangle(gA, gB.as<Triangle>());

    // collection dispatch
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
    return distanceGeometryCollectionToGeometry(gB, gA);

  case TYPE_SOLID:
    BOOST_THROW_EXCEPTION(NotImplementedException(
        (boost::format("distance(%s,%s) is not implemented") %
         gA.geometryType() % gB.geometryType())
            .str()));
  }

  BOOST_ASSERT(false);
  return 0;
}

auto
distancePolygonPolygon(const Polygon &gA, const Polygon &gB) -> double
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  if (intersects(gA, gB, NoValidityCheck())) {
    return 0.0;
  }

  double dMin = std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < gA.numRings(); i++) {
    for (size_t j = 0; j < gB.numRings(); j++) {
      double const d = distanceLineStringLineString(gA.ringN(i), gB.ringN(j));

      if (d < dMin) {
        dMin = d;
      }
    }
  }

  return dMin;
}

auto
distancePolygonTriangle(const Polygon &gA, const Triangle &gB) -> double
{
  return distancePolygonPolygon(gA, gB.toPolygon());
}

auto
distanceTriangleGeometry(const Triangle &gA, const Geometry &gB) -> double
{
  return distancePolygonGeometry(gA.toPolygon(), gB);
}

/**
 * @struct Circle
 * @brief Represents a 2D circle with a center and a radius.
 */
struct Circle {
  /**
   * @brief Main constructor.
   * @param radius The radius of the circle.
   * @param center The center of the circle.
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
distanceGeometryCollectionToGeometry(const Geometry &gA, const Geometry &gB)
    -> double
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  // if bounding spheres (BS) of gB and gAi don't intersect and
  // if the closest point of BS(gAj) is further than the farest
  // point of BS(gAi) there is no need to compute the distance(gAj, gB)
  // since it will be greater than distance(gAi, gB)
  //
  // The aim is not to find the minimal bounding sphere, but a good enough
  // sphere than encloses all points
  std::set<size_t> noTest;

  if (true) {
    std::vector<Circle> bcA;

    for (size_t i = 0; i < gA.numGeometries(); i++) {
      bcA.push_back(boundingCircle(gA.geometryN(i)));
    }

    Circle const bcB(boundingCircle(gB));

    if (bcB.isEmpty()) {
      return std::numeric_limits<double>::infinity();
    }

    std::vector<size_t> noIntersect;

    for (size_t i = 0; i < gA.numGeometries(); i++) {
      if (bcA[i].isEmpty()) {
        continue;
      }

      const double l2 =
          CGAL::to_double((bcB.center() - bcA[i].center()).squared_length());

      if (std::pow(bcB.radius() + bcA[i].radius(), 2) < l2) {
        noIntersect.push_back(i);
      }
    }

    for (size_t i = 0; i < noIntersect.size(); i++) {
      const double li = std::sqrt(CGAL::to_double(
          (bcA[noIntersect[i]].center() - bcB.center()).squared_length()));

      for (size_t j = i; j < noIntersect.size(); j++) {
        const double lj = std::sqrt(CGAL::to_double(
            (bcA[noIntersect[j]].center() - bcB.center()).squared_length()));

        if (li + bcA[noIntersect[i]].radius() <
            lj - bcA[noIntersect[j]].radius()) {
          noTest.insert(noIntersect[j]);
        } else if (lj + bcA[noIntersect[j]].radius() <
                   li - bcA[noIntersect[i]].radius()) {
          noTest.insert(noIntersect[i]);
        }
      }
    }

    // if (!noTest.empty()) std::cout << "pruning " << noTest.size() << "/" <<
    // gA.numGeometries() << "\n";
  }

  double dMin = std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < gA.numGeometries(); i++) {
    if (noTest.end() != noTest.find(i)) {
      continue;
    }

    dMin = std::min(dMin, distance(gA.geometryN(i), gB));
  }

  return dMin;
}

auto
distancePointSegment(const Point &p, const Point &a, const Point &b) -> double
{
  // empty already checked
  BOOST_ASSERT(!p.isEmpty());
  BOOST_ASSERT(!a.isEmpty());
  BOOST_ASSERT(!b.isEmpty());

  return CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(
      p.toPoint_2(), Segment_2(a.toPoint_2(), b.toPoint_2()))));
}

auto
distanceSegmentSegment(const Point &a, const Point &b, const Point &c,
                       const Point &d) -> double
{
  // empty already checked
  BOOST_ASSERT(!a.isEmpty());
  BOOST_ASSERT(!b.isEmpty());
  BOOST_ASSERT(!c.isEmpty());
  BOOST_ASSERT(!d.isEmpty());

  return CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(
      CGAL::Segment_2<Kernel>(a.toPoint_2(), b.toPoint_2()),
      CGAL::Segment_2<Kernel>(c.toPoint_2(), d.toPoint_2()))));
}

#endif // ifndef DOXYGEN_SHOULD_SKIP_THIS
/// @} end of private section

} // namespace SFCGAL::algorithm
