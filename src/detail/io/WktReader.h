// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_IO_WKTREADER_H_
#define SFCGAL_IO_WKTREADER_H_

#include <sstream>

#include "SFCGAL/config.h"

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/PreparedGeometry.h"

#include "SFCGAL/detail/tools/InputStreamReader.h"

namespace SFCGAL {
namespace detail {
namespace io {

/**
 * Reader for Well-Known Text (WKT) geometry format.
 *
 * Parses WKT strings to create SFCGAL Geometry objects. Supports all
 * standard geometry types including Points, LineStrings, Polygons,
 * and collections, as well as 3D and measured coordinates.
 *
 * @warning M is ignored
 * @todo take M in account?
 */
class SFCGAL_API WktReader {
public:
  /**
   * Constructor: read WKT from input stream
   *
   * @param inputStream input stream from which WKT will be read
   */
  WktReader(std::istream &inputStream);

  /**
   * read an SRID, if present
   *
   * @return SRID value
   */
  auto
  readSRID() -> srid_t;

  /**
   * read a geometry from a string
   *
   * @warning returns new instance
   * @return pointer to newly created Geometry object
   */
  auto
  readGeometry() -> Geometry *;

  /**
   * read geometry type
   *
   * @return GeometryType of the next geometry in the stream
   */
  auto
  readGeometryType() -> GeometryType;

  /**
   * read coordinate type [Z][M]
   *
   * @return CoordinateType of the geometry
   */
  auto
  readCoordinateType() -> CoordinateType;

  /**
   * Read Point content from WKT
   *
   * ex : (1.0 2.0 14.0)
   *
   * @param point Point object to fill with coordinates
   */
  void
  readInnerPoint(Point &point);

  /**
   * Read LineString content from WKT
   *
   * ex : (1.0 2.0,1.0,6.0)
   *
   * @param lineString LineString object to fill with coordinates
   */
  void
  readInnerLineString(LineString &lineString);

  /**
   * Read Polygon content from WKT
   *
   * ex : ((30 10, 10 20, 20 40, 40 40, 30 10))
   *
   * @param polygon Polygon object to fill with coordinates
   */
  void
  readInnerPolygon(Polygon &polygon);

  /**
   * Read Triangle content from WKT
   *
   * @param triangle Triangle object to fill with coordinates
   */
  void
  readInnerTriangle(Triangle &triangle);

  /**
   * Read MultiPoint content from WKT
   *
   * ex : (0.0 1.0,5.0 6.0) or ((0.0 4.0),(5.0 6.0))
   *
   * @param multiPoint MultiPoint object to fill with coordinates
   */
  void
  readInnerMultiPoint(MultiPoint &multiPoint);

  /**
   * Read MultiLineString content from WKT
   *
   * @param multiLineString MultiLineString object to fill with coordinates
   */
  void
  readInnerMultiLineString(MultiLineString &multiLineString);

  /**
   * Read MultiPolygon content from WKT
   *
   * @param multiPolygon MultiPolygon object to fill with coordinates
   */
  void
  readInnerMultiPolygon(MultiPolygon &multiPolygon);

  /**
   * Read GeometryCollection content from WKT
   *
   * @param collection GeometryCollection object to fill
   */
  void
  readInnerGeometryCollection(GeometryCollection &collection);

  /**
   * Read TriangulatedSurface content from WKT
   *
   * @param tin TriangulatedSurface object to fill
   */
  void
  readInnerTriangulatedSurface(TriangulatedSurface &tin);

  /**
   * Read PolyhedralSurface content from WKT
   *
   * @param surface PolyhedralSurface object to fill
   */
  void
  readInnerPolyhedralSurface(PolyhedralSurface &surface);

  /**
   * Read Solid content from WKT
   *
   * @param solid Solid object to fill
   */
  void
  readInnerSolid(Solid &solid);

  /**
   * Read MultiSolid content from WKT
   *
   * @param multiSolid MultiSolid object to fill
   */
  void
  readInnerMultiSolid(MultiSolid &multiSolid);

  /**
   * Read NURBSCurve content from WKT (NEW ISO FORMAT).
   *
   * Supports syntax:
   * - NURBSCURVE(degree, (points))
   * - NURBSCURVE(degree, (points), (weights))
   * - NURBSCURVE(degree, (points), (weights), (knots))
   *
   * @param g The NURBSCurve geometry to populate with parsed data.
   */
  void
  readInnerNURBSCurve(NURBSCurve &g);

  /**
   * Read vector of weights from WKT format: (w1, w2, w3, ...)
   * @return Vector of weight values
   */
  auto
  readWeightsVector() -> std::vector<Kernel::FT>;

  /**
   * Read vector of knots from WKT format: (k1, k2, k3, ...)
   * @return Vector of knot values
   */
  auto
  readKnotsVector() -> std::vector<Kernel::FT>;

  /**
   * Read degree value (unsigned integer)
   * @return NURBS curve degree
   */
  auto
  readDegree() -> unsigned int;

  /**
   * Read coordinate from WKT
   *
   * @todo ZM management
   * @param point Point object to fill with coordinate
   * @return true if coordinate was read successfully, false otherwise
   */
  auto
  readPointCoordinate(Point &point) -> bool;

private:
  /**
   * input stream
   */
  tools::InputStreamReader _reader;

  /**
   * actually reading 3D ?
   */
  bool _is3D;
  /**
   * actually reading Measured ?
   */
  bool _isMeasured;

  /**
   * returns default parse error message
   *
   * @return string containing error message
   */
  auto
  parseErrorMessage() -> std::string;
};

} // namespace io
} // namespace detail
} // namespace SFCGAL

#endif
