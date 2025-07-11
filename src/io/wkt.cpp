// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/io/wkt.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/detail/io/WktReader.h"
#include "SFCGAL/detail/io/WktWriter.h"
#include "SFCGAL/detail/tools/CharArrayBuffer.h"

using namespace SFCGAL::detail::io;
using namespace SFCGAL::tools;

namespace SFCGAL::io {

auto
readWkt(std::istream &s) -> std::unique_ptr<Geometry>
{
  WktReader wktReader(s);
  return std::unique_ptr<Geometry>(wktReader.readGeometry());
}

auto
readWkt(const std::string &s) -> std::unique_ptr<Geometry>
{
  std::istringstream        iss(s);
  WktReader                 wktReader(iss);
  std::unique_ptr<Geometry> geom(wktReader.readGeometry());

  char extra = 0;
  if (iss >> extra) {
    std::string const remaining(s.substr(int(iss.tellg()) - 1));
    throw WktParseException("Extra characters in WKT: " + remaining);
  }
  return geom;
}

auto
readWkt(const char *str, size_t len) -> std::unique_ptr<Geometry>
{
  CharArrayBuffer           buf(str, str + len);
  std::istream              istr(&buf);
  WktReader                 wktReader(istr);
  std::unique_ptr<Geometry> geom(wktReader.readGeometry());
  char                      extra = 0;
  if (istr >> extra) {
    std::string const remaining(str + int(istr.tellg()) - 1, str + len);
    throw WktParseException("Extra characters in WKT: " + remaining);
  }
  return geom;
}

} // namespace SFCGAL::io
