// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2023, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/io/wkb.h>

#include <SFCGAL/Exception.h>
#include <SFCGAL/detail/io/WkbReader.h>
#include <SFCGAL/detail/io/WkbWriter.h>
#include <SFCGAL/detail/tools/CharArrayBuffer.h>

using namespace SFCGAL::detail::io;

namespace SFCGAL::io {

///
///
///
// auto
// readWkb(std::istream &s) -> std::unique_ptr<Geometry>
// {
//   WkbReader wkbReader(s);
//   return std::unique_ptr<Geometry>(wkbReader.readGeometry());
// }

///
///
///
auto
readWkb(const std::string &s) -> std::unique_ptr<Geometry>
{
  WkbReader                 wkbReader(s);
  std::unique_ptr<Geometry> geom(wkbReader.readWkb());

  return geom;
}

///
///
///
// auto
// readWkb(const char *str, size_t len) -> std::unique_ptr<Geometry>
// {
//   CharArrayBuffer           buf(str, str + len);
//   std::istream              istr(&buf);
//   WkbReader                 wkbReader(istr);
//   std::unique_ptr<Geometry> geom(wkbReader.readGeometry());
//   char                      extra = 0;
//   if (istr >> extra) {
//     std::string remaining(str + int(istr.tellg()) - 1, str + len);
//     throw WkbParseException("Extra characters in WKB: " + remaining);
//   }
//   return geom;
// }

} // namespace SFCGAL::io
