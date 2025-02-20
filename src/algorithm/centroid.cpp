// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/centroid.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"

#include "SFCGAL/algorithm/area.h"
#include "SFCGAL/algorithm/isValid.h"

#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangle_2.h>

#include <CGAL/Plane_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Triangle_3.h>

#include "SFCGAL/Exception.h"
#include <boost/format.hpp>

namespace SFCGAL::algorithm {

using Point_2    = CGAL::Point_2<SFCGAL::Kernel>;
using Triangle_2 = CGAL::Triangle_2<SFCGAL::Kernel>;
using Polygon_2  = CGAL::Polygon_2<SFCGAL::Kernel>;
using Vector_2   = CGAL::Vector_2<SFCGAL::Kernel>;

using Point_3    = CGAL::Point_3<SFCGAL::Kernel>;
using Triangle_3 = CGAL::Triangle_3<SFCGAL::Kernel>;
using Plane_3    = CGAL::Plane_3<SFCGAL::Kernel>;
using Vector_3   = CGAL::Vector_3<SFCGAL::Kernel>;

///
///
///
template <typename T>
std::function<WeightedCentroid(const T &)> weightedCentroidLambda =
    [](const T &g) -> WeightedCentroid {
  SFCGAL::Kernel::FT totalArea = 0.0;
  SFCGAL::Kernel::FT totalM    = 0.0;
  Vector_3           totalWeightedCentroid;

  for (typename T::const_iterator ite = g.begin(); ite != g.end(); ite++) {
    // compute geometry area+centroid
    WeightedCentroid wc = weightedCentroid(*ite);

    // update totals
    if (ite == g.begin())
      totalWeightedCentroid = wc.area * wc.centroid;
    else
      totalWeightedCentroid += wc.area * wc.centroid;
    totalM += wc.area * wc.m;
    totalArea += wc.area;
  }

  totalWeightedCentroid /= totalArea;
  totalM /= totalArea;

  return {totalArea, totalWeightedCentroid, totalM};
};

///
///
///
auto
centroid(const Geometry &g) -> std::unique_ptr<Point>
{
  if (g.isEmpty()) {
    BOOST_THROW_EXCEPTION(
        InappropriateGeometryException("No point in geometry."));
  }

  WeightedCentroid wCent = weightedCentroid(g);

  Point out;
  if (g.is3D())
    out = Point(wCent.centroid.x(), wCent.centroid.y(), wCent.centroid.z());
  else
    out = Point(wCent.centroid.x(), wCent.centroid.y());

  if (g.isMeasured())
    out.setM(CGAL::to_double(wCent.m));

  return std::make_unique<Point>(out);
}

auto
weightedCentroid(const Geometry &g) -> WeightedCentroid
{
  WeightedCentroid wCent;
  switch (g.geometryTypeId()) {
  case TYPE_POINT:
    wCent = WeightedCentroid(0.0, g.as<Point>().toVector_3());
    break;
  case TYPE_LINESTRING:
    wCent = weightedCentroid(g.as<LineString>());
    break;

  case TYPE_POLYGON:
    wCent = weightedCentroid(g.as<Polygon>());
    break;

  case TYPE_TRIANGLE:
    wCent = weightedCentroid(g.as<Triangle>());
    break;

  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_GEOMETRYCOLLECTION:
    wCent = weightedCentroid(g.as<GeometryCollection>());
    break;

  case TYPE_TRIANGULATEDSURFACE:
    wCent = weightedCentroid(g.as<TriangulatedSurface>());
    break;

  case TYPE_POLYHEDRALSURFACE:
    // wCent = weightedCentroid(g.as<PolyhedralSurface>());

  case TYPE_SOLID:
  case TYPE_MULTISOLID:
    BOOST_THROW_EXCEPTION(NotImplementedException(
        "Centroid of solid gemotry is not implemented."));
  }

  return wCent;
  BOOST_THROW_EXCEPTION(
      Exception((boost::format("Unexpected geometry type (%s) in "
                               "SFCGAL::algorithm::weightedCentroid") %
                 g.geometryType())
                    .str()));
}

///
///
///
auto
weightedCentroid(const Triangle &g) -> WeightedCentroid
{
  g.toTriangle_2().area();
  return weightedCentroid(g.vertex(0), g.vertex(1), g.vertex(2));
}

///
///
///
auto
weightedCentroid(const Point &a, const Point &b, const Point &c)
    -> WeightedCentroid
{
  // std::cout << "============== weightedCentroid(Triangle)\n";
  // compute triangle area
  SFCGAL::Kernel::FT area = SFCGAL::Kernel().compute_area_2_object()(
      a.toPoint_2(), b.toPoint_2(), c.toPoint_2());

  // compute triangle centroid
  Vector_3 out = c.toVector_3();
  out += b.toVector_3();
  out += a.toVector_3();
  out /= 3.0;

  SFCGAL::Kernel::FT m = 0.0;
  if (a.isMeasured() && b.isMeasured() && c.isMeasured()) {
    m = (a.m() + b.m() + c.m()) / 3.0;
  }

  // std::cout << "Triangle weightedCentroid: " << out << "\n";
  // std::cout << "Triangle area: " << area << "\n";

  return {area, out, m};
}

///
///
///
auto
weightedCentroid(const LineString &g) -> WeightedCentroid
{
  // std::cout << "============== weightedCentroid(LineString)\n";
  SFCGAL::Kernel::FT totalArea = 0.0;
  SFCGAL::Kernel::FT totalM    = 0.0;
  Vector_3           totalWeightedCentroid;

  if (g.isClosed()) { // ie. a polygon
    for (size_t i = 1; i < g.numPoints() - 2; i++) {
      // compute triangle area+centroid
      WeightedCentroid wc =
          weightedCentroid(g.pointN(0), g.pointN(i), g.pointN(i + 1));

      // update totals
      if (i == 1)
        totalWeightedCentroid = wc.area * wc.centroid;
      else
        totalWeightedCentroid += wc.area * wc.centroid;
      totalM += wc.area * wc.m;
      totalArea += wc.area;
      // std::cout << "Polygon weightedCentroid: " << wc.area * wc.centroid <<
      // "\n"; std::cout << "Polygon area: " << wc.area << "\n"; std::cout <<
      // "Polygon totalweightedCentroid: " << totalWeightedCentroid << "\n";
      // std::cout << "Polygon totalArea: " << totalArea << "\n";
    }

  } else { // ie. a linestring
    for (size_t i = 0; i < g.numPoints() - 1; i++) {
      Point    a  = g.pointN(i);
      Point    b  = g.pointN(i + 1);
      Vector_3 ba = b.toPoint_3() - a.toPoint_3();

      // TODO: 3D length! should be 2D?
      SFCGAL::Kernel::FT len = CGAL::sqrt(CGAL::to_double(ba.squared_length()));
      if (i == 0)
        totalWeightedCentroid = (a.toVector_3() + (ba / 2.0)) * len;
      else
        totalWeightedCentroid += (a.toVector_3() + (ba / 2.0)) * len;
      if (a.isMeasured() && b.isMeasured())
        totalM += (a.m() + b.m()) * 0.5 * len;
      totalArea += len;
      // std::cout << "LineString weightedCentroid: " << (a.toVector_3() + (ba
      // / 2.0)) * len << "\n"; std::cout << "LineString len: " << len << "\n";
      // std::cout << "LineString totalweightedCentroid: " <<
      // totalWeightedCentroid << "\n"; std::cout << "LineString totalArea: " <<
      // totalArea << "\n";
    }
  }

  if (totalArea == 0.0)
    BOOST_THROW_EXCEPTION(InappropriateGeometryException(
        "SFCGAL::algorithm::Centroid of LineString without 2D area is not "
        "valid."));

  totalWeightedCentroid /= totalArea;
  totalM /= totalArea;

  return {totalArea, totalWeightedCentroid, totalM};
}

///
///
///
auto
weightedCentroid(const Polygon &g) -> WeightedCentroid
{
  // std::cout << "============== weightedCentroid(Polygon)\n";
  SFCGAL::Kernel::FT totalArea = 0.0;
  SFCGAL::Kernel::FT totalM    = 0.0;
  Vector_3           totalWeightedCentroid;

  for (size_t i = 0; i < g.numRings(); i++) {
    try {
      WeightedCentroid ringCentroid = weightedCentroid(g.ringN(i));
      // Kernel::FT const ringCentroid = CGAL::abs( signedCentroid( g.ringN( i )
      // ) );

      if (i == 0) {
        // exterior ring
        totalArea             = ringCentroid.area;
        totalM                = ringCentroid.area * ringCentroid.m;
        totalWeightedCentroid = ringCentroid.area * ringCentroid.centroid;
      } else {
        // interior ring
        totalArea -= ringCentroid.area;
        totalM -= ringCentroid.area * ringCentroid.m;
        totalWeightedCentroid -= ringCentroid.area * ringCentroid.centroid;
      }
      // std::cout << "sub Polygon weightedCentroid: " << ringCentroid.area *
      // ringCentroid.centroid
      //           << "\n";
      // std::cout << "sub Polygon area: " << ringCentroid.area << "\n";
      // std::cout << "sub Polygon totalweightedCentroid: " <<
      // totalWeightedCentroid << "\n"; std::cout << "sub Polygon totalArea: "
      // << totalArea << "\n";

    } catch (InappropriateGeometryException &e) {
      BOOST_THROW_EXCEPTION(InappropriateGeometryException(
          (boost::format("SFCGAL::algorithm::Centroid of Polygon failed at "
                         "ring %d with cause: \n%s") %
           i % e.what())
              .str()));
    }
  }

  totalWeightedCentroid /= totalArea;
  totalM /= totalArea;

  return {totalArea, totalWeightedCentroid, totalM};
}

///
///
///
auto
weightedCentroid(const GeometryCollection &g) -> WeightedCentroid
{
  return weightedCentroidLambda<GeometryCollection>(g);
}

///
///
///
auto
weightedCentroid(const TriangulatedSurface &g) -> WeightedCentroid
{
  return weightedCentroidLambda<TriangulatedSurface>(g);
}

///
///
///
auto
weightedCentroid(const PolyhedralSurface &g) -> WeightedCentroid
{
  return weightedCentroidLambda<PolyhedralSurface>(g);
}

} // namespace SFCGAL::algorithm
