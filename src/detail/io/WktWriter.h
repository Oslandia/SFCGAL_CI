/**
 *   SFCGAL
 *
 *   Copyright (C) 2012-2013 Oslandia <infos@oslandia.com>
 *   Copyright (C) 2012-2013 IGN (http://www.ign.fr)
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Library General Public
 *   License as published by the Free Software Foundation; either
 *   version 2 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Library General Public License for more details.

 *   You should have received a copy of the GNU Library General Public
 *   License along with this library; if not, see
 <http://www.gnu.org/licenses/>.
 */

#ifndef _SFCGAL_IO_WKTWRITER_H_
#define _SFCGAL_IO_WKTWRITER_H_

#include <sstream>

#include <SFCGAL/Geometry.h>
#include <SFCGAL/config.h>

namespace SFCGAL {
namespace detail {
namespace io {

/**
 * Writer for WKT
 *
 * @warning Triangles are transformed into polygons
 */
class SFCGAL_API WktWriter {
public:
  WktWriter(std::ostream &s);

  /**
   * @todo replace with visitor dispatch
   */
  void
  write(const Geometry &g, bool exact = false);

protected:
  /**
   * write coordinate type (""|" Z"|" ZM")
   */
  void
  writeCoordinateType(const Geometry &g);

  void
  writeCoordinate(const Point &g);

  void
  write(const Point &g);
  void
  writeInner(const Point &g);

  void
  write(const LineString &g);
  void
  writeInner(const LineString &g);

  void
  write(const Polygon &g);
  void
  writeInner(const Polygon &g);

  void
  write(const GeometryCollection &g);

  void
  write(const MultiPoint &g);
  void
  write(const MultiLineString &g);
  void
  write(const MultiPolygon &g);
  void
  write(const MultiSolid &g);

  void
  write(const Triangle &g);
  void
  writeInner(const Triangle &g);

  void
  write(const TriangulatedSurface &g);

  void
  write(const PolyhedralSurface &g);
  void
  writeInner(const PolyhedralSurface &g);

  void
  write(const Solid &g);
  void
  writeInner(const Solid &g);

  // for recursive call use
  void
  writeRec(const Geometry &g);

private:
  std::ostream &_s;
  bool          _exactWrite;
};

} // namespace io
} // namespace detail
} // namespace SFCGAL

#endif
