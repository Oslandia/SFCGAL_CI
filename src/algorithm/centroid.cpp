// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/centroid.h"

#include "SFCGAL/Curve.h"
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

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

template <typename T>
std::function<WeightedCentroid(const T &, bool)> weightedCentroidLambda =
    [](const T &g, bool enable3DComputation) -> WeightedCentroid {
  SFCGAL::Kernel::FT totalArea = 0.0;
  SFCGAL::Kernel::FT totalM    = 0.0;
  Vector_3           totalWeightedCentroid;

  for (typename T::const_iterator ite = g.begin(); ite != g.end(); ite++) {
    // compute geometry area+centroid
    WeightedCentroid wc = weightedCentroid(*ite, enable3DComputation);

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

/// @} end of private section

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

auto
centroid(const Geometry &geom) -> std::unique_ptr<Point>
{
  if (geom.isEmpty()) {
    BOOST_THROW_EXCEPTION(
        InappropriateGeometryException("No point in geometry."));
  }

  WeightedCentroid wCent = weightedCentroid(geom);

  Point out;
  if (geom.is3D())
    out = Point(wCent.centroid.x(), wCent.centroid.y(), wCent.centroid.z());
  else
    out = Point(wCent.centroid.x(), wCent.centroid.y());

  if (geom.isMeasured())
    out.setM(CGAL::to_double(wCent.m));

  return std::make_unique<Point>(out);
}

auto
centroid3D(const Geometry &geom) -> std::unique_ptr<Point>
{
  if (geom.isEmpty()) {
    BOOST_THROW_EXCEPTION(
        InappropriateGeometryException("No point in geometry."));
  }

  WeightedCentroid wCent = weightedCentroid(geom, true);

  Point out;
  if (geom.is3D())
    out = Point(wCent.centroid.x(), wCent.centroid.y(), wCent.centroid.z());
  else
    out = Point(wCent.centroid.x(), wCent.centroid.y());

  if (geom.isMeasured())
    out.setM(CGAL::to_double(wCent.m));

  return std::make_unique<Point>(out);
}
auto
weightedCentroid(const Geometry &geom, bool enable3DComputation)
    -> WeightedCentroid
{
  WeightedCentroid wCent;
  switch (geom.geometryTypeId()) {
  case TYPE_POINT:
    wCent = WeightedCentroid(0.0, geom.as<Point>().toVector_3(),
                             (geom.isMeasured() ? geom.as<Point>().m() : 0.0));
    break;
  case TYPE_LINESTRING:
    wCent = weightedCentroid(geom.as<LineString>(), enable3DComputation);
    break;

  case TYPE_NURBSCURVE:
    wCent = weightedCentroid(geom.as<Curve>(), enable3DComputation);
    break;

  case TYPE_POLYGON:
    wCent = weightedCentroid(geom.as<Polygon>(), enable3DComputation);
    break;

  case TYPE_TRIANGLE:
    wCent = weightedCentroid(geom.as<Triangle>(), enable3DComputation);
    break;

  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION:
    wCent =
        weightedCentroid(geom.as<GeometryCollection>(), enable3DComputation);
    break;

  case TYPE_TRIANGULATEDSURFACE:
    wCent =
        weightedCentroid(geom.as<TriangulatedSurface>(), enable3DComputation);
    break;

  case TYPE_POLYHEDRALSURFACE:
    wCent = weightedCentroid(geom.as<PolyhedralSurface>(), enable3DComputation);
    break;

  case TYPE_SOLID:
    wCent = weightedCentroid(geom.as<Solid>(), enable3DComputation);
    break;
  }

  return wCent;
  BOOST_THROW_EXCEPTION(
      Exception((boost::format("Unexpected geometry type (%s) in "
                               "SFCGAL::algorithm::weightedCentroid") %
                 geom.geometryType())
                    .str()));
}

auto
weightedCentroid(const Triangle &triangle, bool enable3DComputation)
    -> WeightedCentroid
{
  return weightedCentroid(triangle.vertex(0), triangle.vertex(1),
                          triangle.vertex(2), enable3DComputation);
}

auto
weightedCentroid(const Point &pta, const Point &ptb, const Point &ptc,
                 bool enable3DComputation) -> WeightedCentroid
{
  // compute triangle area
  SFCGAL::Kernel::FT area = 0.0;
  if (enable3DComputation) {
    // compute 3D area
    Triangle_3 t(pta.toPoint_3(), ptb.toPoint_3(), ptc.toPoint_3());
    area = CGAL::sqrt(CGAL::to_double(t.squared_area()));
  } else {
    // compute 2D signed area
    area = SFCGAL::Kernel().compute_area_2_object()(
        pta.toPoint_2(), ptb.toPoint_2(), ptc.toPoint_2());
  }

  // compute triangle centroid
  Vector_3 out = ptc.toVector_3();
  out += ptb.toVector_3();
  out += pta.toVector_3();
  out /= 3.0;

  SFCGAL::Kernel::FT m = 0.0;
  if (pta.isMeasured() && ptb.isMeasured() && ptc.isMeasured()) {
    m = (pta.m() + ptb.m() + ptc.m()) / 3.0;
  }

  return {area, out, m};
}

auto
weightedCentroid(const Curve &curve, bool enable3DComputation)
    -> WeightedCentroid
{
  // Convert curve to LineString approximation for centroid calculation
  auto lineString = curve.toLineString(); // default parameters
  if (!lineString || lineString->isEmpty()) {
    return {}; // Return empty centroid for invalid/empty curves
  }
  return weightedCentroid(*lineString, enable3DComputation);
}

auto
weightedCentroid(const LineString &lineString, bool enable3DComputation)
    -> WeightedCentroid
{
  SFCGAL::Kernel::FT totalArea = 0.0;
  SFCGAL::Kernel::FT totalM    = 0.0;
  Vector_3           totalWeightedCentroid;

  if (lineString.isClosed()) { // ie. a polygon
    for (size_t i = 1; i < lineString.numPoints() - 2; i++) {
      // compute triangle area+centroid
      WeightedCentroid wc =
          weightedCentroid(lineString.pointN(0), lineString.pointN(i),
                           lineString.pointN(i + 1), enable3DComputation);

      // update totals
      if (i == 1)
        totalWeightedCentroid = wc.area * wc.centroid;
      else
        totalWeightedCentroid += wc.area * wc.centroid;
      totalM += wc.area * wc.m;
      totalArea += wc.area;
    }

  } else { // ie. a linestring
    for (size_t i = 0; i < lineString.numPoints() - 1; i++) {
      Point    a  = lineString.pointN(i);
      Point    b  = lineString.pointN(i + 1);
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

auto
weightedCentroid(const Polygon &polygon, bool enable3DComputation)
    -> WeightedCentroid
{
  SFCGAL::Kernel::FT totalArea = 0.0;
  SFCGAL::Kernel::FT totalM    = 0.0;
  Vector_3           totalWeightedCentroid;

  for (size_t i = 0; i < polygon.numRings(); i++) {
    try {
      WeightedCentroid ringCentroid =
          weightedCentroid(polygon.ringN(i), enable3DComputation);

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

auto
weightedCentroid(const GeometryCollection &collection, bool enable3DComputation)
    -> WeightedCentroid
{
  return weightedCentroidLambda<GeometryCollection>(collection,
                                                    enable3DComputation);
}

auto
weightedCentroid(const TriangulatedSurface &tin, bool enable3DComputation)
    -> WeightedCentroid
{
  return weightedCentroidLambda<TriangulatedSurface>(tin, enable3DComputation);
}

auto
weightedCentroid(const PolyhedralSurface &surface, bool enable3DComputation)
    -> WeightedCentroid
{
  return weightedCentroidLambda<PolyhedralSurface>(surface,
                                                   enable3DComputation);
}

auto
weightedCentroid(const Solid &solid, bool enable3DComputation)
    -> WeightedCentroid
{
  SFCGAL::Kernel::FT totalArea = 0.0;
  SFCGAL::Kernel::FT totalM    = 0.0;
  Vector_3           totalWeightedCentroid;

  for (size_t i = 0; i < solid.numShells(); i++) {
    try {
      WeightedCentroid shellCentroid =
          weightedCentroid(solid.shellN(i), enable3DComputation);

      if (i == 0) {
        // exterior shell
        totalArea             = shellCentroid.area;
        totalM                = shellCentroid.area * shellCentroid.m;
        totalWeightedCentroid = shellCentroid.area * shellCentroid.centroid;
      } else {
        // interior shell
        totalArea -= shellCentroid.area;
        totalM -= shellCentroid.area * shellCentroid.m;
        totalWeightedCentroid -= shellCentroid.area * shellCentroid.centroid;
      }
    } catch (InappropriateGeometryException &e) {
      BOOST_THROW_EXCEPTION(InappropriateGeometryException(
          (boost::format("SFCGAL::algorithm::Centroid of Solid failed at "
                         "shell %d with cause: \n%s") %
           i % e.what())
              .str()));
    }
  }

  totalWeightedCentroid /= totalArea;
  totalM /= totalArea;

  return {totalArea, totalWeightedCentroid, totalM};
}

} // namespace SFCGAL::algorithm
