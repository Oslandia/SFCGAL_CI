// Copyright (c) 2023-2023, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef _SFCGAL_IO_WKBWRITER_H_
#define _SFCGAL_IO_WKBWRITER_H_

#include <boost/endian/conversion.hpp>
#include <cstddef>
#include <vector>

#include <SFCGAL/Geometry.h>
#include <SFCGAL/config.h>

using srid_t = uint32_t;

namespace SFCGAL::detail::io {

/**
 * Writer for WKB
 *
 */
class SFCGAL_API WkbWriter {
public:
  /**
   * write WKB for a geometry
   * wkbOrder is the native endianness by default.
   */
  void
  write(const Geometry      &g,
        boost::endian::order wkbOrder = boost::endian::order::native);

  /**
   * write EWKB for a geometry
   * wkbOrder is the native endianness by default.
   */
  void
  write(const Geometry &g, const srid_t &srid,
        boost::endian::order wkbOrder = boost::endian::order::native);

  /**
   * Export (E)WKB as a string
   *
   * asHex will transform the string with \x suffix.
   */
  auto
  toString(bool asHex = false) -> std::string;

private:
  /**
   * Dedicated method to write the geometry type into _wkb data
   */
  void
  writeGeometryType(const Geometry &g, boost::endian::order wkbOrder =
                                           boost::endian::order::native);

  /**
   * Dedicated method to write Point into _wkb data
   */
  void
  writeInner(const Point         &g,
             boost::endian::order wkbOrder = boost::endian::order::native);

  /**
   * Dedicated method to write LineString into _wkb data
   */
  void
  writeInner(const LineString    &g,
             boost::endian::order wkbOrder = boost::endian::order::native);
  /**
   * Dedicated method to write Ring into _wkb data
   *
   * This method is shared by LineString and Polygon.
   */
  void
  writeInnerRing(const LineString    &g,
                 boost::endian::order wkbOrder = boost::endian::order::native);

  /**
   * Dedicated method to write Polygon into _wkb data
   */
  void
  writeInner(const Polygon       &g,
             boost::endian::order wkbOrder = boost::endian::order::native);

  /**
   * Dedicated method to write GeometryCollection into _wkb data
   */
  void
  writeInner(const GeometryCollection &g,
             boost::endian::order      wkbOrder = boost::endian::order::native);

  /**
   * Dedicated method to write Multi Geometries into _wkb data
   * Multi Geometries are: MultiPoint, MultiLineString, MultiPolygon,
   * PolyhedralSurface and TriangulatedSurface.
   */
  template <typename M, typename G>
  void
  writeInner(const M &g, boost::endian::order wkbOrder);

  /**
   * Dedicated method to write Polygon into _wkb data
   */
  void
  writeInner(const Triangle      &g,
             boost::endian::order wkbOrder = boost::endian::order::native);

  /**
   * Dedicated method to write Polygon into _wkb data
   */
  void
  writeCoordinate(const Point         &g,
                  boost::endian::order wkbOrder = boost::endian::order::native);

  /**
   * Method to write Geometry into _wkb data
   * Only for recursive call use
   */
  void
  writeRec(const Geometry      &g,
           boost::endian::order wkbOrder = boost::endian::order::native);

  std::vector<std::byte> _wkb;

  srid_t _srid;

  bool _useSrid = false;
  bool _isEWKB  = false;
};

} // namespace SFCGAL::detail::io

#endif
