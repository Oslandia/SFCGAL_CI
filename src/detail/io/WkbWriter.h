// Copyright (c) 2023-2023, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef _SFCGAL_IO_WKBWRITER_H_
#define _SFCGAL_IO_WKBWRITER_H_

#include <boost/endian/conversion.hpp>
#include <cstddef>
#include <vector>

#include <SFCGAL/Geometry.h>
#include <SFCGAL/config.h>

typedef uint32_t srid_t;

namespace SFCGAL::detail::io {

/**
 * Writer for WKB
 *
 */
class SFCGAL_API WkbWriter {
public:
  void
  write(const Geometry      &g,
        boost::endian::order wkbOrder = boost::endian::order::native);

  void
  write(const Geometry &g, const srid_t &srid,
        boost::endian::order wkbOrder = boost::endian::order::native);

  auto
  toString(const bool asHex = false) -> std::string;

protected:
  void
  writeInner(const Point         &g,
             boost::endian::order wkbOrder = boost::endian::order::native);

  void
  writeInner(const LineString    &g,
             boost::endian::order wkbOrder = boost::endian::order::native);
  void
  writeInnerRing(const LineString    &g,
                 boost::endian::order wkbOrder = boost::endian::order::native);

  void
  writeInner(const Polygon       &g,
             boost::endian::order wkbOrder = boost::endian::order::native);
  ;

  void
  writeInner(const GeometryCollection &g,
             boost::endian::order      wkbOrder = boost::endian::order::native);
  ;

  // for Multi Geometries (MultiPoint, MultiLineString, MultiPolygon, PolyhedralSurface and TriangulatedSurface)
template <typename M, typename G>
void
writeInner(const M &g, boost::endian::order wkbOrder);

  void
  writeInner(const MultiSolid    &g,
             boost::endian::order wkbOrder = boost::endian::order::native);
  ;

  void
  writeInner(const Triangle      &g,
             boost::endian::order wkbOrder = boost::endian::order::native);
  ;

  void
  writeInner(const Solid         &g,
             boost::endian::order wkbOrder = boost::endian::order::native);
  ;

private:
  void
  writeCoordinateType(const Geometry &g, boost::endian::order wkbOrder =
                                             boost::endian::order::native);
  ;

  void
  writeCoordinate(const Point         &g,
                  boost::endian::order wkbOrder = boost::endian::order::native);
  ;

  void
  writeGeometryType(const Geometry &g, boost::endian::order wkbOrder =
                                           boost::endian::order::native);
  ;
  // for recursive call use
  void
  writeRec(const Geometry      &g,
           boost::endian::order wkbOrder = boost::endian::order::native);
  ;

  std::vector<std::byte> _wkb;

  srid_t _srid;

  bool _useSrid = false;
  bool _isEWKB  = false;
};

} // namespace SFCGAL::detail::io

#endif
