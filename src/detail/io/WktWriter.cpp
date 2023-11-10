// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/detail/io/WktWriter.h>

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

#include <boost/exception/all.hpp>
#include <exception>

namespace SFCGAL::detail::io {

namespace impl {
auto
writeFT(std::ostream &s, const CGAL::Gmpq &ft) -> std::ostream &
{
  s << ft;
  return s;
}

#ifdef CGAL_USE_GMPXX
auto
writeFT(std::ostream &s, const mpq_class &ft) -> std::ostream &
{
  s << ft.get_num() << "/" << ft.get_den();
  return s;
}
#endif
} // namespace impl

///
///
///
WktWriter::WktWriter(std::ostream &s) : _s(s), _exactWrite(false) {}

void
WktWriter::writeRec(const Geometry &g)
{
  switch (g.geometryTypeId()) {
  case TYPE_POINT:
    write(g.as<Point>());
    return;

  case TYPE_LINESTRING:
    write(g.as<LineString>());
    return;

  case TYPE_POLYGON:
    write(g.as<Polygon>());
    return;

  case TYPE_GEOMETRYCOLLECTION:
    write(g.as<GeometryCollection>());
    return;

  case TYPE_MULTIPOINT:
    write(g.as<MultiPoint>());
    return;

  case TYPE_MULTILINESTRING:
    write(g.as<MultiLineString>());
    return;

  case TYPE_MULTIPOLYGON:
    write(g.as<MultiPolygon>());
    return;

  case TYPE_TRIANGLE:
    write(g.as<Triangle>());
    return;

  case TYPE_TRIANGULATEDSURFACE:
    write(g.as<TriangulatedSurface>());
    return;

  case TYPE_POLYHEDRALSURFACE:
    write(g.as<PolyhedralSurface>());
    return;

  case TYPE_SOLID:
    write(g.as<Solid>());
    return;

  case TYPE_MULTISOLID:
    write(g.as<MultiSolid>());
    return;
  }

  std::ostringstream oss;
  oss << "WktWriter : '" << g.geometryType() << "' is not supported";
  BOOST_THROW_EXCEPTION(std::runtime_error(oss.str()));
}

///
///
///
void
WktWriter::write(const Geometry &g, bool exact)
{
  _exactWrite = exact;
  writeRec(g);
}

///
///
///
void
WktWriter::writeCoordinateType(const Geometry &g)
{
  if (g.is3D() && !g.isMeasured()) {
    _s << " Z";
  } else if (!g.is3D() && g.isMeasured()) {
    _s << " M";
  } else if (g.is3D() && g.isMeasured()) {
    _s << " ZM";
  }
}

static auto
fixZeroNeg(double val, int precision) -> double
{
  if (std::abs(val) < std::pow(10, -precision)) {
    return 0;
  }
  return val;
}

///
///
///
void
WktWriter::writeCoordinate(const Point &g)
{
  if (_exactWrite) {
    impl::writeFT(_s, CGAL::exact(g.x())) << " ";
    impl::writeFT(_s, CGAL::exact(g.y()));

    if (g.is3D()) {
      _s << " ";
      impl::writeFT(_s, CGAL::exact(g.z()));
    }
  } else {
    _s << fixZeroNeg(CGAL::to_double(g.x()), _s.precision()) << " "
       << fixZeroNeg(CGAL::to_double(g.y()), _s.precision());

    if (g.is3D()) {
      _s << " " << fixZeroNeg(CGAL::to_double(g.z()), _s.precision());
    }
  }

  // m coordinate
  if (g.isMeasured()) {
    _s << " " << fixZeroNeg(CGAL::to_double(g.m()), _s.precision());
  }
}

///
///
///
void
WktWriter::write(const Point &g)
{
  _s << "POINT";
  writeCoordinateType(g);

  if (g.isEmpty()) {
    _s << " EMPTY";
    return;
  }

  writeInner(g);
}

///
///
///
void
WktWriter::writeInner(const Point &g)
{
  if (g.isEmpty()) {
    _s << "EMPTY";
    return;
  }

  _s << "(";
  writeCoordinate(g);
  _s << ")";
}

///
///
///
void
WktWriter::write(const LineString &g)
{
  _s << "LINESTRING";
  writeCoordinateType(g);

  if (g.isEmpty()) {
    _s << " EMPTY";
    return;
  }

  writeInner(g);
}

///
///
///
void
WktWriter::writeInner(const LineString &g)
{
  _s << "(";

  for (size_t i = 0; i < g.numPoints(); i++) {
    if (i != 0) {
      _s << ",";
    }

    writeCoordinate(g.pointN(i));
  }

  _s << ")";
}

///
///
///
void
WktWriter::write(const Polygon &g)
{
  _s << "POLYGON";
  writeCoordinateType(g);

  if (g.isEmpty()) {
    _s << " EMPTY";
    return;
  }

  writeInner(g);
}

///
///
///
void
WktWriter::writeInner(const Polygon &g)
{
  _s << "(";
  writeInner(g.exteriorRing());

  for (size_t i = 0; i < g.numInteriorRings(); i++) {
    _s << ",";
    writeInner(g.interiorRingN(i));
  }

  _s << ")";
}

///
///
///
void
WktWriter::write(const GeometryCollection &g)
{
  _s << "GEOMETRYCOLLECTION";
  writeCoordinateType(g);

  if (g.isEmpty()) {
    _s << " EMPTY";
    return;
  }

  _s << "(";

  for (size_t i = 0; i < g.numGeometries(); i++) {
    if (i != 0) {
      _s << ",";
    }

    writeRec(g.geometryN(i));
  }

  _s << ")";
}

///
///
///
void
WktWriter::write(const MultiPoint &g)
{
  _s << "MULTIPOINT";
  writeCoordinateType(g);

  if (g.isEmpty()) {
    _s << " EMPTY";
    return;
  }

  _s << "(";

  for (size_t i = 0; i < g.numGeometries(); i++) {
    if (i != 0) {
      _s << ",";
    }

    writeInner(g.geometryN(i).as<Point>());
  }

  _s << ")";
}

///
///
///
void
WktWriter::write(const MultiLineString &g)
{
  _s << "MULTILINESTRING";
  writeCoordinateType(g);

  if (g.isEmpty()) {
    _s << " EMPTY";
    return;
  }

  _s << "(";

  for (size_t i = 0; i < g.numGeometries(); i++) {
    if (i != 0) {
      _s << ",";
    }

    writeInner(g.geometryN(i).as<LineString>());
  }

  _s << ")";
}

///
///
///
void
WktWriter::write(const MultiPolygon &g)
{
  _s << "MULTIPOLYGON";
  writeCoordinateType(g);

  if (g.isEmpty()) {
    _s << " EMPTY";
    return;
  }

  _s << "(";

  for (size_t i = 0; i < g.numGeometries(); i++) {
    if (i != 0) {
      _s << ",";
    }

    writeInner(g.geometryN(i).as<Polygon>());
  }

  _s << ")";
}

///
///
///
void
WktWriter::write(const MultiSolid &g)
{
  _s << "MULTISOLID";
  writeCoordinateType(g);

  if (g.isEmpty()) {
    _s << " EMPTY";
    return;
  }

  _s << "(";

  for (size_t i = 0; i < g.numGeometries(); i++) {
    if (i != 0) {
      _s << ",";
    }

    writeInner(g.geometryN(i).as<Solid>());
  }

  _s << ")";
}

///
///
///
void
WktWriter::write(const Triangle &g)
{
  _s << "TRIANGLE";
  writeCoordinateType(g);

  if (g.isEmpty()) {
    _s << " EMPTY";
    return;
  }

  writeInner(g);
}

///
///
///
void
WktWriter::writeInner(const Triangle &g)
{
  _s << "(";
  _s << "(";

  // close triangle
  for (size_t i = 0; i < 4; i++) {
    if (i != 0) {
      _s << ",";
    }

    writeCoordinate(g.vertex(i));
  }

  _s << ")";
  _s << ")";
}

///
///
///
void
WktWriter::write(const TriangulatedSurface &g)
{
  _s << "TIN";
  writeCoordinateType(g);

  if (g.isEmpty()) {
    _s << " EMPTY";
    return;
  }

  _s << "("; // begin TIN

  for (size_t i = 0; i < g.numGeometries(); i++) {
    if (i != 0) {
      _s << ",";
    }

    writeInner(g.geometryN(i));
  }

  _s << ")"; // end TIN
}

///
///
///
void
WktWriter::write(const PolyhedralSurface &g)
{
  _s << "POLYHEDRALSURFACE";
  writeCoordinateType(g);

  if (g.isEmpty()) {
    _s << " EMPTY";
    return;
  }

  writeInner(g);
}

///
///
///
void
WktWriter::writeInner(const PolyhedralSurface &g)
{
  _s << "("; // begin POLYHEDRALSURFACE

  for (size_t i = 0; i < g.numPolygons(); i++) {
    if (i != 0) {
      _s << ",";
    }

    writeInner(g.polygonN(i));
  }

  _s << ")"; // end POLYHEDRALSURFACE
}

///
///
///
void
WktWriter::write(const Solid &g)
{
  _s << "SOLID";
  writeCoordinateType(g);

  if (g.isEmpty()) {
    _s << " EMPTY";
    return;
  }

  writeInner(g);
}

///
///
///
void
WktWriter::writeInner(const Solid &g)
{
  _s << "("; // begin SOLID
  writeInner(g.exteriorShell());

  for (size_t i = 0; i < g.numInteriorShells(); i++) {
    _s << ",";
    writeInner(g.interiorShellN(i));
  }

  _s << ")"; // end SOLID
}

} // namespace SFCGAL
