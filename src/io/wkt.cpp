// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/io/wkt.h>

#include <SFCGAL/Exception.h>
#include <SFCGAL/detail/io/WktReader.h>
#include <SFCGAL/detail/io/WktWriter.h>
#include <SFCGAL/detail/tools/CharArrayBuffer.h>

using namespace SFCGAL::detail::io;

namespace SFCGAL {
namespace io {

///
///
///
auto
readWkt(std::istream &s) -> std::unique_ptr<Geometry>
{
  WktReader wktReader(s);
  return std::unique_ptr<Geometry>(wktReader.readGeometry());
}

///
///
///
auto
readWkt(const std::string &s) -> std::unique_ptr<Geometry>
{
  std::istringstream        iss(s);
  WktReader                 wktReader(iss);
  std::unique_ptr<Geometry> geom(wktReader.readGeometry());

  char extra = 0;
  if (iss >> extra) {
    std::string remaining(s.substr(int(iss.tellg()) - 1));
    throw WktParseException("Extra characters in WKT: " + remaining);
  }
  return geom;
}

///
///
///
auto
readWkt(const char *str, size_t len) -> std::unique_ptr<Geometry>
{
  CharArrayBuffer           buf(str, str + len);
  std::istream              istr(&buf);
  WktReader                 wktReader(istr);
  std::unique_ptr<Geometry> geom(wktReader.readGeometry());
  char                      extra = 0;
  if (istr >> extra) {
    std::string remaining(str + int(istr.tellg()) - 1, str + len);
    throw WktParseException("Extra characters in WKT: " + remaining);
  }
  return geom;
}

} // namespace io
} // namespace SFCGAL
