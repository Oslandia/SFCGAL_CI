// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _SFCGAL_IO_VTK_H_
#define _SFCGAL_IO_VTK_H_

#include <fstream>
#include <iostream>
#include <sstream>

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/MultiSolid.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>

using namespace SFCGAL;

namespace SFCGAL {
namespace io {
// print each ring has a different polygon
inline void
vtk(const Polygon &poly, const std::string &file)
{
  std::stringstream pointStr;
  std::stringstream polyStr;
  size_t            numPoints = 0;
  size_t            numData   = 0;

  for (size_t r = 0; r != poly.numRings(); ++r) {
    polyStr << poly.ringN(r).numPoints();

    for (size_t p = 0; p != poly.ringN(r).numPoints(); ++p) {
      pointStr << poly.ringN(r).pointN(p).x() << " "
               << poly.ringN(r).pointN(p).y() << " "
               << poly.ringN(r).pointN(p).z() << "\n";
      polyStr << " " << numPoints;
      ++numPoints;
      ++numData;
    }

    ++numData;
    polyStr << "\n";
  }

  std::ofstream out(file.c_str());
  out << "# vtk DataFile Version 1.0\n"
      << "Polygon output\n"
      << "ASCII\n"
      << "\n"
      << "DATASET POLYDATA\n"
      << "POINTS " << numPoints << " float\n"
      << pointStr.str() << "POLYGONS " << poly.numRings() << " " << numData
      << "\n"
      << polyStr.str();
}

template <typename MultiPolygonOrPolyhedraSurface>
void
vtk(const MultiPolygonOrPolyhedraSurface &multiPoly, const std::string &file)
{
  std::stringstream pointStr;
  std::stringstream polyStr;
  size_t            numPoints = 0;
  size_t            numRings  = 0;
  size_t            numData   = 0;

  for (size_t i = 0; i != multiPoly.numGeometries(); ++i) {
    const Polygon &poly = multiPoly.polygonN(i);

    for (size_t r = 0; r != poly.numRings(); ++r) {
      polyStr << poly.ringN(r).numPoints();

      for (size_t p = 0; p != poly.ringN(r).numPoints(); ++p) {
        pointStr << poly.ringN(r).pointN(p).x() << " "
                 << poly.ringN(r).pointN(p).y() << " "
                 << poly.ringN(r).pointN(p).z() << "\n";
        polyStr << " " << numPoints;
        ++numPoints;
        ++numData;
      }

      ++numData;
      ++numRings;
      polyStr << "\n";
    }
  }

  std::ofstream out(file.c_str());
  out << "# vtk DataFile Version 1.0\n"
      << "Polygon output\n"
      << "ASCII\n"
      << "\n"
      << "DATASET POLYDATA\n"
      << "POINTS " << numPoints << " float\n"
      << pointStr.str() << "POLYGONS " << numRings << " " << numData << "\n"
      << polyStr.str();
}

inline void
vtk(const Triangle &tri, const std::string &file)
{
  vtk(tri.toPolygon(), file);
}

inline void
vtk(const TriangulatedSurface &s, const std::string &file)
{
  std::stringstream pointStr;
  std::stringstream polyStr;
  size_t            numPoints = 0;
  size_t            numTri    = 0;
  size_t            numData   = 0;

  for (size_t i = 0; i != s.numTriangles(); ++i) {
    const Triangle &tri = s.triangleN(i);
    polyStr << 3;

    for (size_t p = 0; p != 3; ++p) {
      pointStr << tri.vertex(p).x() << " " << tri.vertex(p).y() << " "
               << tri.vertex(p).z() << "\n";
      polyStr << " " << numPoints;
      ++numPoints;
      ++numData;
    }

    ++numData;
    ++numTri;
    polyStr << "\n";
  }

  std::ofstream out(file.c_str());
  out << "# vtk DataFile Version 1.0\n"
      << "Polygon output\n"
      << "ASCII\n"
      << "\n"
      << "DATASET POLYDATA\n"
      << "POINTS " << numPoints << " float\n"
      << pointStr.str() << "POLYGONS " << numTri << " " << numData << "\n"
      << polyStr.str();
}

inline void
vtk(const Geometry &g, const std::string &file)
{
  switch (g.geometryTypeId()) {
  case TYPE_POLYGON:
    vtk(g.as<Polygon>(), file);
    return;

  case TYPE_MULTIPOLYGON:
    vtk(g.as<MultiPolygon>(), file);
    return;

  case TYPE_TRIANGLE:
    vtk(g.as<Triangle>(), file);
    return;

  case TYPE_TRIANGULATEDSURFACE:
    vtk(g.as<TriangulatedSurface>(), file);
    return;

  case TYPE_POLYHEDRALSURFACE:
    vtk(g.as<PolyhedralSurface>(), file);
    return;

  case TYPE_SOLID:
    vtk(g.as<Solid>().exteriorShell(), file);
    return;

  case TYPE_POINT:
  case TYPE_LINESTRING:
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION:
    BOOST_THROW_EXCEPTION(NotImplementedException(
        (boost::format("vtk(%s, file) is not implemented") % g.geometryType())
            .str()));
  }
}
} // namespace io
} // namespace SFCGAL

#endif
