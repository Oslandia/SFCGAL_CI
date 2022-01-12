// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _SFCGAL_ALGORITHM_PLANE_H_
#define _SFCGAL_ALGORITHM_PLANE_H_

#include <boost/format.hpp>

//#include <SFCGAL/detail/ublas.h>

#include <SFCGAL/Exception.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/algorithm/normal.h>
#include <SFCGAL/detail/GetPointsVisitor.h>

namespace SFCGAL {
namespace algorithm {

/**
 * @brief Test if a 3D plane can be extracted from a Polygon
 * @ingroup public_api
 */
template <typename Kernel>
bool
hasPlane3D(const Polygon &polygon, CGAL::Point_3<Kernel> &a,
           CGAL::Point_3<Kernel> &b, CGAL::Point_3<Kernel> &c)
{
  typedef CGAL::Point_3<Kernel> Point_3;

  const LineString &exteriorRing = polygon.exteriorRing();

  /*
   * look for 3 non collinear points
   */
  size_t n = 0;

  for (size_t i = 0; i < exteriorRing.numPoints(); i++) {
    Point_3 p = exteriorRing.pointN(i).toPoint_3();

    if (n == 0) {
      a = p;
      n++;
    } else if (n == 1 && a != p) {
      b = p;
      n++;
    } else if (n == 2 && !CGAL::collinear(a, b, p)) {
      c = p;
      n++;
      return true;
    }
  }

  BOOST_ASSERT(n < 3);
  return false;
}

/**
 * Test if a 3D plane can be extracted from a Polygon
 */
template <typename Kernel>
bool
hasPlane3D(const Polygon &polygon)
{
  // temporary arguments
  CGAL::Point_3<Kernel> a, b, c;
  return hasPlane3D(polygon, a, b, c);
}

/**
 * Get 3 non collinear points from a Polygon
 */
template <typename Kernel>
void
plane3D(const Polygon &polygon, CGAL::Point_3<Kernel> &a,
        CGAL::Point_3<Kernel> &b, CGAL::Point_3<Kernel> &c)
{
  if (!hasPlane3D(polygon, a, b, c)) {
    BOOST_THROW_EXCEPTION(
        Exception((boost::format("can't find plane for Polygon '%1%'") %
                   polygon.asText(3))
                      .str()));
  }
}

/**
 * Returns the oriented 3D plane of a polygon (supposed to be planar).
 * May return degenerate plane.
 */
template <typename Kernel>
CGAL::Plane_3<Kernel>
plane3D(const Polygon &polygon)
{
  CGAL::Vector_3<Kernel> nrml = normal3D<Kernel>(polygon, true);

  return CGAL::Plane_3<Kernel>(polygon.exteriorRing().pointN(0).toPoint_3(),
                               nrml);
}

struct Plane3DInexactUnsafe {
};

/**
 * Returns the oriented 3D plane of a polygon (supposed to be planar) - inexact
 * version.
 * @warning Will divide by zero if polygon is degenerate.
 * @warning result is rounded to double (avoid huge expression tree).
 */
template <typename Kernel>
CGAL::Plane_3<Kernel>
plane3D(const Polygon &polygon, const Plane3DInexactUnsafe &)
{
  CGAL::Vector_3<Kernel> nrml = normal3D<Kernel>(polygon, false);

  const double nrm = std::sqrt(CGAL::to_double(nrml.squared_length()));
  nrml = CGAL::Vector_3<Kernel>(nrml.x() / nrm, nrml.y() / nrm, nrml.z() / nrm);

  return CGAL::Plane_3<Kernel>(polygon.exteriorRing().pointN(0).toPoint_3(),
                               nrml);
}

/**
 * Returns the oriented 3D plane of a polygon (supposed to be planar).
 * This is legacy code for SFCGAL users and should be deprecated.
 * @warning result is rounded to double if exact is false (avoid huge expression
 * tree).
 * @warning Will divide by zero if polygon is degenerate. This maintains the
 * previous behaviour.
 */
template <typename Kernel>
CGAL::Plane_3<Kernel>
plane3D(const Polygon &polygon, bool exact)
{
  if (exact)
    return plane3D<Kernel>(polygon);
  else
    return plane3D<Kernel>(polygon, Plane3DInexactUnsafe());
}

/**
 * Test if all points of a geometry lie in the same plane
 * @ingroup detail
 */
template <typename Kernel>
bool
isPlane3D(const Geometry &geom, const double &toleranceAbs)
{
  if (geom.isEmpty()) {
    return true;
  }

  using namespace SFCGAL::detail;
  GetPointsVisitor v;
  const_cast<Geometry &>(geom).accept(v);

  if (v.points.size() == 0) {
    return true;
  }

  // the present approach is to find a good plane by:
  // - computing the centroid C of the point set
  // - finding the farest point F from C
  // - finding the farest point G from (CF)
  // - we define the unit normal N to the plane from CFxCG
  // - we check that points Xi are in the plane CXi.N < tolerance
  //
  // note that we could compute the covarence matrix of the points and use SVD
  // but we would need a lib for that, and it may be overkill

  typedef CGAL::Vector_3<Kernel> Vector_3;

  const GetPointsVisitor::const_iterator end = v.points.end();

  // centroid
  Vector_3 c(0, 0, 0);
  int      numPoint = 0;

  for (GetPointsVisitor::const_iterator x = v.points.begin(); x != end; ++x) {
    c = c + (*x)->toVector_3();
    ++numPoint;
  }

  BOOST_ASSERT(numPoint);
  c = c / numPoint;

  // farest point from centroid
  Vector_3            f             = c;
  typename Kernel::FT maxDistanceSq = 0;

  for (GetPointsVisitor::const_iterator x = v.points.begin(); x != end; ++x) {
    const Vector_3            cx  = (*x)->toVector_3() - c;
    const typename Kernel::FT dSq = cx * cx;

    if (dSq > maxDistanceSq) {
      f             = (*x)->toVector_3();
      maxDistanceSq = dSq;
    }
  }

  if (std::sqrt(CGAL::to_double(maxDistanceSq)) < toleranceAbs) {
    // std::cout << "all points in the same location\n";
    return true;
  }

  // farest point from line
  Vector_3       g  = c;
  const Vector_3 cf = f - c; // direction of (CF)
  maxDistanceSq     = 0;     // watch out, we reuse the variable

  for (GetPointsVisitor::const_iterator x = v.points.begin(); x != end; ++x) {
    const Vector_3 cx = (*x)->toVector_3() - c;
    const Vector_3 cp =
        (cx * cf) * cf / cf.squared_length(); // projection of x on line (CF)
    const typename Kernel::FT dSq = (cx - cp).squared_length();

    if (dSq > maxDistanceSq) {
      g             = (*x)->toVector_3();
      maxDistanceSq = dSq;
    }
  }

  if (std::sqrt(CGAL::to_double(maxDistanceSq)) < toleranceAbs) {
    // std::cout << "all points aligned\n";
    return true;
  }

  const Vector_3 n = CGAL::cross_product(cf, g - c);

  const Vector_3 nNormed = n / std::sqrt(CGAL::to_double(n.squared_length()));

  for (GetPointsVisitor::const_iterator x = v.points.begin(); x != end; ++x) {
    const Vector_3 cx = (*x)->toVector_3() - c;

    if (std::abs(CGAL::to_double(cx * nNormed)) > toleranceAbs) {
      // std::cout << "point out of plane\n";
      return false;
    }
  }

  // std::cout << "plane general case\n";
  return true;
}

} // namespace algorithm
} // namespace SFCGAL

#endif
