// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/distance3d.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/detail/tools/Log.h"

#include "SFCGAL/algorithm/intersects.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/detail/GetPointsVisitor.h"
#include "SFCGAL/detail/transform/AffineTransform3.h"
#include "SFCGAL/triangulate/triangulatePolygon.h"

using squared_distance_t = SFCGAL::Kernel::FT;

using Point_3    = SFCGAL::Kernel::Point_3;
using Vector_3   = SFCGAL::Kernel::Vector_3;
using Segment_3  = SFCGAL::Kernel::Segment_3;
using Triangle_3 = SFCGAL::Kernel::Triangle_3;
using Plane_3    = SFCGAL::Kernel::Plane_3;

namespace SFCGAL::algorithm {

///
///
///
auto
distance3D(const Geometry &gA, const Geometry &gB, NoValidityCheck /*unused*/)
    -> double
{
  // SFCGAL_DEBUG( boost::format("dispatch distance3D(%s,%s)") % gA.asText() %
  // gB.asText() );

  switch (gA.geometryTypeId()) {
  case TYPE_POINT:
    return distancePointGeometry3D(gA.as<Point>(), gB);

  case TYPE_LINESTRING:
    return distanceLineStringGeometry3D(gA.as<LineString>(), gB);

  case TYPE_POLYGON:
    // SFCGAL_DEBUG( boost::format("gA is a Polygon (%s)") % gA.geometryTypeId()
    // );
    return distancePolygonGeometry3D(gA.as<Polygon>(), gB);

  case TYPE_TRIANGLE:
    // SFCGAL_DEBUG( boost::format("gA is a Triangle (%s)") %
    // gA.geometryTypeId() );
    return distanceTriangleGeometry3D(gA.as<Triangle>(), gB);

  case TYPE_SOLID:
    return distanceSolidGeometry3D(gA.as<Solid>(), gB);

    // collection dispatch
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
    return distanceGeometryCollectionToGeometry3D(gB, gA);
  }

  BOOST_THROW_EXCEPTION(
      Exception((boost::format("distance3D(%s,%s) is not implemented") %
                 gA.geometryType() % gB.geometryType())
                    .str()));
}

auto
distance3D(const Geometry &gA, const Geometry &gB) -> double
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(gA);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(gB);

  return distance3D(gA, gB, NoValidityCheck());
}
///
///
///
auto
distancePointGeometry3D(const Point &gA, const Geometry &gB) -> double
{
  // SFCGAL_DEBUG( boost::format("dispatch distancePointGeometry3D(%s,%s)") %
  // gA.asText() % gB.asText() );

  switch (gB.geometryTypeId()) {
  case TYPE_POINT:
    return distancePointPoint3D(gA, gB.as<Point>());

  case TYPE_LINESTRING:
    return distancePointLineString3D(gA, gB.as<LineString>());

  case TYPE_TRIANGLE:
    return distancePointTriangle3D(gA, gB.as<Triangle>());

  case TYPE_POLYGON:
    return distancePointPolygon3D(gA, gB.as<Polygon>());

  case TYPE_SOLID:
    return distancePointSolid3D(gA, gB.as<Solid>());

  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
    return distanceGeometryCollectionToGeometry3D(gB, gA);
  }

  BOOST_THROW_EXCEPTION(
      Exception((boost::format("distance3D(%s,%s) is not implemented") %
                 gA.geometryType() % gB.geometryType())
                    .str()));
}

///
///
///
auto
distancePointPoint3D(const Point &gA, const Point &gB) -> double
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  return CGAL::sqrt(
      CGAL::to_double(CGAL::squared_distance(gA.toPoint_3(), gB.toPoint_3())));
}

///
///
///
auto
distancePointLineString3D(const Point &gA, const LineString &gB) -> double
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  double dMin = std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < gB.numSegments(); i++) {
    dMin = std::min(dMin,
                    distancePointSegment3D(gA, gB.pointN(i), gB.pointN(i + 1)));
  }

  return dMin;
}

///
///
///
auto
distancePointTriangle3D(const Point &gA, const Triangle &gB) -> double
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  return distancePointTriangle3D(gA, gB.vertex(0), gB.vertex(1), gB.vertex(2));
}

///
///
///
auto
distancePointPolygon3D(const Point &gA, const Polygon &gB) -> double
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  TriangulatedSurface triangulateSurfaceB;
  triangulate::triangulatePolygon3D(gB, triangulateSurfaceB);
  return distanceGeometryCollectionToGeometry3D(triangulateSurfaceB, gA);
}

///
///
///
auto
distancePointSolid3D(const Point &gA, const Solid &gB) -> double
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  if (intersects3D(gA, gB, NoValidityCheck())) {
    return 0.0;
  }

  double dMin = std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < gB.numShells(); i++) {
    dMin = std::min(dMin,
                    distanceGeometryCollectionToGeometry3D(gB.shellN(i), gA));
  }

  return dMin;
}

///
///
///
auto
distanceLineStringGeometry3D(const LineString &gA, const Geometry &gB) -> double
{
  // SFCGAL_DEBUG( boost::format("dispatch distanceLineStringGeometry3D(%s,%s)")
  // % gA.asText() % gB.asText() );

  switch (gB.geometryTypeId()) {
  case TYPE_POINT:
    return distancePointLineString3D(gB.as<Point>(), gA); // symetric

  case TYPE_LINESTRING:
    return distanceLineStringLineString3D(gA, gB.as<LineString>());

  case TYPE_TRIANGLE:
    return distanceLineStringTriangle3D(gA, gB.as<Triangle>());

  case TYPE_POLYGON:
    return distanceLineStringPolygon3D(gA, gB.as<Polygon>());

  case TYPE_SOLID:
    return distanceLineStringSolid3D(gA, gB.as<Solid>());

  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
    return distanceGeometryCollectionToGeometry3D(gB, gA);
  }

  BOOST_THROW_EXCEPTION(
      Exception((boost::format("distance3D(%s,%s) is not implemented") %
                 gA.geometryType() % gB.geometryType())
                    .str()));
}

///
///
///
auto
distanceLineStringLineString3D(const LineString &gA, const LineString &gB)
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
                      distanceSegmentSegment3D(gA.pointN(i), gA.pointN(i + 1),
                                               gB.pointN(j), gB.pointN(j + 1)));
    }
  }

  return dMin;
}

///
///
///
auto
distanceLineStringTriangle3D(const LineString &gA, const Triangle &gB) -> double
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  double dMin = std::numeric_limits<double>::infinity();

  const Point &tA = gB.vertex(0);
  const Point &tB = gB.vertex(1);
  const Point &tC = gB.vertex(2);

  for (size_t i = 0; i < gA.numSegments(); i++) {
    dMin = std::min(dMin, distanceSegmentTriangle3D(
                              gA.pointN(i), gA.pointN(i + 1), tA, tB, tC));
  }

  return dMin;
}

///
///
///
auto
distanceLineStringPolygon3D(const LineString &gA, const Polygon &gB) -> double
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  TriangulatedSurface triangulateSurfaceB;
  triangulate::triangulatePolygon3D(gB, triangulateSurfaceB);
  return distanceGeometryCollectionToGeometry3D(triangulateSurfaceB, gA);
}

///
///
///
auto
distanceLineStringSolid3D(const LineString &gA, const Solid &gB) -> double
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  if (intersects(gA, gB, NoValidityCheck())) {
    return 0.0;
  }

  double dMin = std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < gB.numShells(); i++) {
    dMin = std::min(dMin, gB.shellN(i).distance3D(gA));
  }

  return dMin;
}

///
///
///
auto
distanceTriangleGeometry3D(const Triangle &gA, const Geometry &gB) -> double
{
  // SFCGAL_DEBUG( boost::format("dispatch distanceTriangleGeometry3D(%s,%s)") %
  // gA.asText() % gB.asText() );

  switch (gB.geometryTypeId()) {
  case TYPE_POINT:
    return distancePointTriangle3D(gB.as<Point>(), gA); // symetric

  case TYPE_LINESTRING:
    return distanceLineStringTriangle3D(gB.as<LineString>(), gA); // symetric

  case TYPE_TRIANGLE:
    return distanceTriangleTriangle3D(gA, gB.as<Triangle>());

  case TYPE_POLYGON:
    return distancePolygonGeometry3D(gB.as<Polygon>(), gA);

  case TYPE_SOLID:
    return distanceTriangleSolid3D(gA, gB.as<Solid>());

  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
    return distanceGeometryCollectionToGeometry3D(gB, gA);
  }

  BOOST_THROW_EXCEPTION(
      Exception((boost::format("distance3D(%s,%s) is not implemented") %
                 gA.geometryType() % gB.geometryType())
                    .str()));
}

///
///
///
auto
distanceTriangleSolid3D(const Triangle &gA, const Solid &gB) -> double
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  if (intersects3D(gA, gB, NoValidityCheck())) {
    return 0.0;
  }

  double dMin = std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < gB.numShells(); i++) {
    dMin = std::min(dMin, gB.shellN(i).distance3D(gA));
  }

  return dMin;
}

///
///
///
auto
distancePolygonGeometry3D(const Polygon &gA, const Geometry &gB) -> double
{
  // SFCGAL_DEBUG( boost::format("dispatch distancePolygonGeometry3D(%s,%s)") %
  // gA.asText() % gB.asText() );

  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  TriangulatedSurface triangulateSurfaceA;
  triangulate::triangulatePolygon3D(gA, triangulateSurfaceA);
  return distanceGeometryCollectionToGeometry3D(triangulateSurfaceA, gB);
}

///
///
///
auto
distanceSolidGeometry3D(const Solid &gA, const Geometry &gB) -> double
{
  // SFCGAL_DEBUG( boost::format("dispatch distanceSolidGeometry3D(%s,%s)") %
  // gA.asText() % gB.asText() );

  switch (gB.geometryTypeId()) {
  case TYPE_POINT:
    return distancePointSolid3D(gB.as<Point>(), gA); // symetric

  case TYPE_LINESTRING:
    return distanceLineStringSolid3D(gB.as<LineString>(), gA); // symetric

  case TYPE_TRIANGLE:
    return distanceTriangleSolid3D(gB.as<Triangle>(), gA); // symetric

  case TYPE_POLYGON:
    return distancePolygonGeometry3D(gB.as<Polygon>(), gA); // symetric

  case TYPE_SOLID:
    return distanceSolidSolid3D(gA, gB.as<Solid>());

  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
    return distanceGeometryCollectionToGeometry3D(gB, gA);
  }

  BOOST_THROW_EXCEPTION(
      Exception((boost::format("distance3D(%s,%s) is not implemented") %
                 gA.geometryType() % gB.geometryType())
                    .str()));
}

///
///
///
auto
distanceSolidSolid3D(const Solid &gA, const Solid &gB) -> double
{
  // SFCGAL_DEBUG( boost::format("dispatch distancePolygonGeometry3D(%s,%s)") %
  // gA.asText() % gB.asText() );

  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  if (intersects(gA, gB, NoValidityCheck())) {
    return 0.0;
  }

  double dMin = std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < gA.numShells(); i++) {
    for (size_t j = 0; j < gB.numShells(); j++) {
      dMin = std::min(dMin, gA.shellN(i).distance3D(gB.shellN(j)));
    }
  }

  return dMin;
}

struct Sphere {
  Sphere(const double &r, CGAL::Vector_3<Kernel> &c)
      : _radius(r), _center(c), _empty(false)
  {
  }
  Sphere() = default;
  [[nodiscard]] auto
  isEmpty() const -> bool
  {
    return _empty;
  }
  [[nodiscard]] auto
  radius() const -> double
  {
    BOOST_ASSERT(!_empty);
    return _radius;
  }
  [[nodiscard]] auto
  center() const -> const CGAL::Vector_3<Kernel> &
  {
    BOOST_ASSERT(!_empty);
    return _center;
  }

private:
  double                 _radius{};
  CGAL::Vector_3<Kernel> _center;
  bool                   _empty{true};
};

auto
boundingSphere(const Geometry &geom) -> const Sphere
{
  if (geom.isEmpty()) {
    return Sphere();
  }

  using namespace SFCGAL::detail;
  GetPointsVisitor v;
  const_cast<Geometry &>(geom).accept(v);

  if (v.points.empty()) {
    return Sphere();
  }

  using Vector_3 = CGAL::Vector_3<Kernel>;

  const auto end = v.points.end();

  // centroid
  Vector_3 c(0, 0, 0);
  int      numPoint = 0;

  for (auto x = v.points.begin(); x != end; ++x) {
    c = c + (*x)->toVector_3();
    ++numPoint;
  }

  BOOST_ASSERT(numPoint);
  c = c / numPoint;

  // farest point from centroid
  Vector_3   f             = c;
  Kernel::FT maxDistanceSq = 0;

  for (auto x = v.points.begin(); x != end; ++x) {
    const Vector_3   cx  = (*x)->toVector_3() - c;
    const Kernel::FT dSq = cx * cx;

    if (dSq > maxDistanceSq) {
      f             = (*x)->toVector_3();
      maxDistanceSq = dSq;
    }
  }

  return Sphere(std::sqrt(CGAL::to_double(maxDistanceSq)), c);
}

///
///
///
auto
distanceGeometryCollectionToGeometry3D(const Geometry &gA, const Geometry &gB)
    -> double
{
  // SFCGAL_DEBUG( boost::format("dispatch
  // distanceGeometryCollectionToGeometry3D(%s,%s)") % gA.asText() % gB.asText()
  // );

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
    std::vector<Sphere> bsA;

    for (size_t i = 0; i < gA.numGeometries(); i++) {
      bsA.push_back(boundingSphere(gA.geometryN(i)));
    }

    Sphere const bsB(boundingSphere(gB));

    if (bsB.isEmpty()) {
      return std::numeric_limits<double>::infinity();
    }

    std::vector<size_t> noIntersect;

    for (size_t i = 0; i < gA.numGeometries(); i++) {
      if (bsA[i].isEmpty()) {
        continue;
      }

      const double l2 =
          CGAL::to_double((bsB.center() - bsA[i].center()).squared_length());

      if (std::pow(bsB.radius() + bsA[i].radius(), 2) < l2) {
        noIntersect.push_back(i);
      }
    }

    for (size_t i = 0; i < noIntersect.size(); i++) {
      const double li = std::sqrt(CGAL::to_double(
          (bsA[noIntersect[i]].center() - bsB.center()).squared_length()));

      for (size_t j = i; j < noIntersect.size(); j++) {
        const double lj = std::sqrt(CGAL::to_double(
            (bsA[noIntersect[j]].center() - bsB.center()).squared_length()));

        if (li + bsA[noIntersect[i]].radius() <
            lj - bsA[noIntersect[j]].radius()) {
          noTest.insert(noIntersect[j]);
        } else if (lj + bsA[noIntersect[j]].radius() <
                   li - bsA[noIntersect[i]].radius()) {
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

    dMin = std::min(dMin, distance3D(gA.geometryN(i), gB));
  }

  return dMin;
}

///
///
///
auto
distancePointSegment3D(const Point &p, const Point &a, const Point &b) -> double
{
  // empty already checked
  BOOST_ASSERT(!p.isEmpty());
  BOOST_ASSERT(!a.isEmpty());
  BOOST_ASSERT(!b.isEmpty());

  return CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(
      p.toPoint_3(), Segment_3(a.toPoint_3(), b.toPoint_3()))));
}

/*
 * missing in CGAL?
 */
auto
squaredDistancePointTriangle3D(const Point_3 &p, const Triangle_3 &abc)
    -> squared_distance_t
{
#if CGAL_VERSION_NR >= 1041001000 // >= 4.10
  return CGAL::squared_distance(p, abc);
#else
  Point_3 a = abc.vertex(0);
  Point_3 b = abc.vertex(1);
  Point_3 c = abc.vertex(2);

  /*
   * project P on ABC plane as projP.
   */
  Point_3 projP = Plane_3(a, b, c).projection(p);

  squared_distance_t dMin;

  if (abc.has_on(projP)) {
    // Is projP is in the triangle, return distance from P to its projection
    // on the plane
    dMin = CGAL::squared_distance(p, projP);
  } else {
    // Else, the distance is the minimum from P to triangle sides
    dMin = CGAL::squared_distance(p, Segment_3(a, b));
    dMin = std::min(dMin, CGAL::squared_distance(p, Segment_3(b, c)));
    dMin = std::min(dMin, CGAL::squared_distance(p, Segment_3(c, a)));
  }

  return dMin;
#endif
}

///
///
///
auto
distancePointTriangle3D(const Point &p_, const Point &a_, const Point &b_,
                        const Point &c_) -> double
{
  // empty already checked
  BOOST_ASSERT(!p_.isEmpty());
  BOOST_ASSERT(!a_.isEmpty());
  BOOST_ASSERT(!b_.isEmpty());
  BOOST_ASSERT(!c_.isEmpty());

  Point_3 const    p = p_.toPoint_3();
  Triangle_3 const abc(a_.toPoint_3(), b_.toPoint_3(), c_.toPoint_3());

  squared_distance_t const dMin = squaredDistancePointTriangle3D(p, abc);
  return CGAL::sqrt(CGAL::to_double(dMin));
}

///
///
///
auto
distanceSegmentSegment3D(const Point &a, const Point &b, const Point &c,
                         const Point &d) -> double
{
  // empty already checked
  BOOST_ASSERT(!a.isEmpty());
  BOOST_ASSERT(!b.isEmpty());
  BOOST_ASSERT(!c.isEmpty());
  BOOST_ASSERT(!d.isEmpty());

  return CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(
      CGAL::Segment_3<Kernel>(a.toPoint_3(), b.toPoint_3()),
      CGAL::Segment_3<Kernel>(c.toPoint_3(), d.toPoint_3()))));
}

auto
squaredDistanceSegmentTriangle3D(const Segment_3 &sAB, const Triangle_3 &tABC)
    -> squared_distance_t
{
  using squared_distance_t = Kernel::FT;

  /*
   * If [sAsB] intersects the triangle (tA,tB,tC), distance is 0.0
   */
  if (CGAL::do_intersect(sAB, tABC)) {
    return 0.0;
  }

  /*
   * else, distance is the min of the following values :
   * - distance from sA to the Triangle
   * - distance from sB to the Triangle
   * - distance from sAB to the side of the Triangles
   */
  squared_distance_t dMin = squaredDistancePointTriangle3D(sAB.vertex(0), tABC);
  dMin = std::min(dMin, squaredDistancePointTriangle3D(sAB.vertex(1), tABC));

  for (int i = 0; i < 3; i++) {
    dMin =
        std::min(dMin, CGAL::squared_distance(
                           sAB, Segment_3(tABC.vertex(i), tABC.vertex(i + 1))));
  }

  return dMin;
}

///
///
///
auto
distanceSegmentTriangle3D(const Point &sA_, const Point &sB_, const Point &tA_,
                          const Point &tB_, const Point &tC_) -> double
{
  // empty already checked
  BOOST_ASSERT(!sA_.isEmpty());
  BOOST_ASSERT(!sB_.isEmpty());
  BOOST_ASSERT(!tA_.isEmpty());
  BOOST_ASSERT(!tB_.isEmpty());
  BOOST_ASSERT(!tC_.isEmpty());

  Point_3 const   sA = sA_.toPoint_3();
  Point_3 const   sB = sB_.toPoint_3();
  Segment_3 const sAB(sA, sB);

  Point_3 const    tA = tA_.toPoint_3();
  Point_3 const    tB = tB_.toPoint_3();
  Point_3 const    tC = tC_.toPoint_3();
  Triangle_3 const tABC(tA, tB, tC);

  squared_distance_t const dMin = squaredDistanceSegmentTriangle3D(sAB, tABC);
  return CGAL::sqrt(CGAL::to_double(dMin));
}

/*
 * missing in CGAL?
 */
auto
squaredDistanceTriangleTriangle3D(const Triangle_3 &triangleA,
                                  const Triangle_3 &triangleB)
    -> squared_distance_t
{
  if (CGAL::do_intersect(triangleA, triangleB)) {
    return {0};
  }

  /*
   * min of distance from A segments to B triangle and B segments to A triangle
   */

  squared_distance_t dMin = squaredDistanceSegmentTriangle3D(
      Segment_3(triangleA.vertex(0), triangleA.vertex(1)), triangleB);
  dMin = std::min(dMin, squaredDistanceSegmentTriangle3D(
                            Segment_3(triangleA.vertex(1), triangleA.vertex(2)),
                            triangleB));
  dMin = std::min(dMin, squaredDistanceSegmentTriangle3D(
                            Segment_3(triangleA.vertex(2), triangleA.vertex(0)),
                            triangleB));

  dMin = std::min(dMin, squaredDistanceSegmentTriangle3D(
                            Segment_3(triangleB.vertex(0), triangleB.vertex(1)),
                            triangleA));
  dMin = std::min(dMin, squaredDistanceSegmentTriangle3D(
                            Segment_3(triangleB.vertex(1), triangleB.vertex(2)),
                            triangleA));
  dMin = std::min(dMin, squaredDistanceSegmentTriangle3D(
                            Segment_3(triangleB.vertex(2), triangleB.vertex(0)),
                            triangleA));

  return dMin;
}

///
///
///
auto
distanceTriangleTriangle3D(const Triangle &gA, const Triangle &gB) -> double
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return std::numeric_limits<double>::infinity();
  }

  Triangle_3 const triangleA = gA.toTriangle_3();
  Triangle_3 const triangleB = gB.toTriangle_3();

  squared_distance_t const dMin =
      squaredDistanceTriangleTriangle3D(triangleA, triangleB);
  return CGAL::sqrt(CGAL::to_double(dMin));
}

} // namespace SFCGAL::algorithm
