// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_IO_WKTWRITER_H_
#define SFCGAL_IO_WKTWRITER_H_

#include <sstream>

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/config.h"

namespace SFCGAL {
namespace detail {
namespace io {

/**
 * Writer for Well-Known Text (WKT) geometry format.
 *
 * Converts SFCGAL Geometry objects to their WKT string representation.
 * Supports all standard geometry types including 3D coordinates and
 * measured values.
 *
 * @warning Triangles are transformed into polygons
 */
class SFCGAL_API WktWriter {
public:
  /**
   * Constructor
   *
   * @param outStream output stream where WKT will be written
   */
  WktWriter(std::ostream &outStream);

  /**
   * Write a geometry as WKT
   *
   * @todo replace with visitor dispatch
   * @param geometry geometry to serialize as WKT
   * @param exact if true, write exact representation (without simplifications)
   */
  void
  write(const Geometry &geometry, bool exact = false);

protected:
  /**
   * Write coordinate type ("" | " Z" | " ZM")
   *
   * @param geometry geometry whose coordinate type is written
   */
  void
  writeCoordinateType(const Geometry &geometry);

  /**
   * Write a single coordinate
   *
   * @param point point whose coordinates are written
   */
  void
  writeCoordinate(const Point &point);

  /**
   * Write a Point geometry
   *
   * @param point Point to serialize
   */
  void
  write(const Point &point);

  /**
   * Write the inner coordinates of a Point
   *
   * @param point Point whose coordinates are written
   */
  void
  writeInner(const Point &point);

  /**
   * Write a LineString geometry
   *
   * @param lineString LineString to serialize
   */
  void
  write(const LineString &lineString);

  /**
   * Write the inner coordinates of a LineString
   *
   * @param lineString LineString whose coordinates are written
   */
  void
  writeInner(const LineString &lineString);

  /**
   * Write a Polygon geometry
   *
   * @param polygon Polygon to serialize
   */
  void
  write(const Polygon &polygon);

  /**
   * Write the inner rings of a Polygon
   *
   * @param polygon Polygon whose rings are written
   */
  void
  writeInner(const Polygon &polygon);

  /**
   * Write a GeometryCollection
   *
   * @param collection GeometryCollection to serialize
   */
  void
  write(const GeometryCollection &collection);

  /**
   * Write a MultiPoint geometry
   *
   * @param multipoint MultiPoint to serialize
   */
  void
  write(const MultiPoint &multipoint);

  /**
   * Write a MultiLineString geometry
   *
   * @param multilinestring MultiLineString to serialize
   */
  void
  write(const MultiLineString &multilinestring);

  /**
   * Write a MultiPolygon geometry
   *
   * @param multipolygon MultiPolygon to serialize
   */
  void
  write(const MultiPolygon &multipolygon);

  /**
   * Write vector of weights in WKT format: (w1, w2, w3, ...)
   */
  void
  writeWeights(const std::vector<Kernel::FT> &weights);

  /**
   * Write vector of knots in WKT format: (k1, k2, k3, ...)
   */
  void
  writeKnots(const std::vector<Kernel::FT> &knots);

  /**
   * Write NURBSCurve to WKT format.
   *
   * Outputs the complete WKT representation including geometry type
   * and coordinate information.
   *
   * @param g The NURBSCurve geometry to write.
   */
  void
  write(const NURBSCurve &g);
  /**
   * Write inner content of NURBSCurve to WKT format.
   *
   * Outputs only the coordinate content without the geometry type prefix,
   * used for nested geometry structures.
   *
   * @param g The NURBSCurve geometry to write.
   */
  void
  writeInner(const NURBSCurve &g);

  // for recursive call use
  /**
   * Write a MultiSolid geometry
   *
   * @param multisolid MultiSolid to serialize
   */
  void
  write(const MultiSolid &multisolid);

  /**
   * Write a Triangle geometry (as Polygon)
   *
   * @param triangle Triangle to serialize
   */
  void
  write(const Triangle &triangle);

  /**
   * Write the inner coordinates of a Triangle
   *
   * @param triangle Triangle whose coordinates are written
   */
  void
  writeInner(const Triangle &triangle);

  /**
   * Write a TriangulatedSurface geometry
   *
   * @param tin TriangulatedSurface to serialize
   */
  void
  write(const TriangulatedSurface &tin);

  /**
   * Write a PolyhedralSurface geometry
   *
   * @param surface PolyhedralSurface to serialize
   */
  void
  write(const PolyhedralSurface &surface);

  /**
   * Write the inner content of a PolyhedralSurface
   *
   * @param surface PolyhedralSurface whose content is written
   */
  void
  writeInner(const PolyhedralSurface &surface);

  /**
   * Write a Solid geometry
   *
   * @param solid Solid to serialize
   */
  void
  write(const Solid &solid);

  /**
   * Write the inner shells of a Solid
   *
   * @param solid Solid whose shells are written
   */
  void
  writeInner(const Solid &solid);

  /**
   * Recursive call to write a geometry
   *
   * @param geometry geometry to serialize recursively
   */
  void
  writeRec(const Geometry &geometry);

private:
  std::ostream &_s;
  bool          _exactWrite = false;
};

} // namespace io
} // namespace detail
} // namespace SFCGAL

#endif
