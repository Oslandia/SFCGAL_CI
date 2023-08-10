// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef _SFCGAL_IO_WKBWRITER_H_
#define _SFCGAL_IO_WKBWRITER_H_

#include <boost/endian/conversion.hpp>
#include <cstddef>
#include <vector>

#include <SFCGAL/Geometry.h>
#include <SFCGAL/config.h>

typedef uint32_t srid_t;

namespace SFCGAL {
namespace detail {
namespace io {

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
  write(const Point         &g,
        boost::endian::order wkbOrder = boost::endian::order::native);
  void
  writeInner(const Point         &g,
             boost::endian::order wkbOrder = boost::endian::order::native);

  void
  write(const LineString    &g,
        boost::endian::order wkbOrder = boost::endian::order::native);
  void
  writeInner(const LineString    &g,
             boost::endian::order wkbOrder = boost::endian::order::native);
  void
  writeInnerRing(const LineString    &g,
                 boost::endian::order wkbOrder = boost::endian::order::native);

  void
  write(const Polygon       &g,
        boost::endian::order wkbOrder = boost::endian::order::native);
  ;
  void
  writeInner(const Polygon       &g,
             boost::endian::order wkbOrder = boost::endian::order::native);
  ;

  void
  write(const GeometryCollection &g,
        boost::endian::order      wkbOrder = boost::endian::order::native);
  ;
  void
  writeInner(const GeometryCollection &g,
             boost::endian::order      wkbOrder = boost::endian::order::native);
  ;

  void
  write(const MultiPoint    &g,
        boost::endian::order wkbOrder = boost::endian::order::native);
  ;
  void
  writeInner(const MultiPoint    &g,
             boost::endian::order wkbOrder = boost::endian::order::native);
  ;
  void
  write(const MultiLineString &g,
        boost::endian::order   wkbOrder = boost::endian::order::native);
  ;
  void
  writeInner(const MultiLineString &g,
             boost::endian::order   wkbOrder = boost::endian::order::native);
  ;
  void
  write(const MultiPolygon  &g,
        boost::endian::order wkbOrder = boost::endian::order::native);
  ;
  void
  writeInner(const MultiPolygon  &g,
             boost::endian::order wkbOrder = boost::endian::order::native);
  ;
  void
  write(const MultiSolid    &g,
        boost::endian::order wkbOrder = boost::endian::order::native);
  ;
  void
  writeInner(const MultiSolid    &g,
             boost::endian::order wkbOrder = boost::endian::order::native);
  ;

  void
  write(const Triangle      &g,
        boost::endian::order wkbOrder = boost::endian::order::native);
  ;
  void
  writeInner(const Triangle      &g,
             boost::endian::order wkbOrder = boost::endian::order::native);
  ;

  void
  write(const TriangulatedSurface &g,
        boost::endian::order       wkbOrder = boost::endian::order::native);
  ;
  void
  writeInner(const TriangulatedSurface &g,
             boost::endian::order wkbOrder = boost::endian::order::native);
  ;

  void
  write(const PolyhedralSurface &g,
        boost::endian::order     wkbOrder = boost::endian::order::native);
  ;
  void
  writeInner(const PolyhedralSurface &g,
             boost::endian::order     wkbOrder = boost::endian::order::native);
  ;

  void
  write(const Solid         &g,
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

} // namespace io
} // namespace detail
} // namespace SFCGAL

#endif
