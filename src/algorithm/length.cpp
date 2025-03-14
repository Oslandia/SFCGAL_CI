// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/length.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"

#include "SFCGAL/Exception.h"

namespace SFCGAL::algorithm {

auto
length(const LineString &g) -> double
{
  double result = 0.0;

  for (size_t i = 0; i < g.numSegments(); i++) {
    Kernel::Segment_2 const segment(g.pointN(i).toPoint_2(),
                                    g.pointN(i + 1).toPoint_2());
    result += CGAL::sqrt(CGAL::to_double(segment.squared_length()));
  }

  return result;
}

auto
length(const GeometryCollection &g) -> double
{
  double result = 0.0;

  for (const auto &it : g) {
    result += length(it);
  }

  return result;
}

auto
length(const Geometry &g) -> double
{
  switch (g.geometryTypeId()) {
  case TYPE_POINT:
    return 0.0;

  case TYPE_LINESTRING:
    return length(g.as<LineString>());

  case TYPE_POLYGON:
    return 0.0;

  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_GEOMETRYCOLLECTION:
    return length(g.as<GeometryCollection>());

  case TYPE_POLYHEDRALSURFACE:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_TRIANGLE:
  case TYPE_SOLID:
  case TYPE_MULTISOLID:
    return 0.0;
  }

  BOOST_THROW_EXCEPTION(
      Exception((boost::format("undefined length for geometry type %s") %
                 g.geometryType())
                    .str()));
  return 0.0;
}

//------ 3D

auto
length3D(const LineString &g) -> double
{
  double result = 0.0;

  for (size_t i = 0; i < g.numSegments(); i++) {
    Kernel::Segment_3 const segment(g.pointN(i).toPoint_3(),
                                    g.pointN(i + 1).toPoint_3());
    result += CGAL::sqrt(CGAL::to_double(segment.squared_length()));
  }

  return result;
}

auto
length3D(const GeometryCollection &g) -> double
{
  double result = 0.0;

  for (const auto &it : g) {
    result += length3D(it);
  }

  return result;
}

auto
length3D(const Geometry &g) -> double
{
  switch (g.geometryTypeId()) {
  case TYPE_POINT:
    return 0.0;

  case TYPE_LINESTRING:
    return length3D(g.as<LineString>());

  case TYPE_POLYGON:
    return 0.0;

  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_GEOMETRYCOLLECTION:
    return length3D(g.as<GeometryCollection>());

  case TYPE_POLYHEDRALSURFACE:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_TRIANGLE:
  case TYPE_SOLID:
  case TYPE_MULTISOLID:
    return 0.0;
  }

  BOOST_THROW_EXCEPTION(
      Exception((boost::format("undefined length for geometry type %s") %
                 g.geometryType())
                    .str()));
  return 0.0;
}

} // namespace SFCGAL::algorithm
