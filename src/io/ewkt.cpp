// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/io/ewkt.h>

#include <SFCGAL/detail/io/WktReader.h>
#include <SFCGAL/detail/io/WktWriter.h>
#include <SFCGAL/detail/tools/CharArrayBuffer.h>

using namespace SFCGAL::detail::io;

namespace SFCGAL {
namespace io {

///
///
///
std::unique_ptr<PreparedGeometry>
readEwkt(std::istream &s)
{
  WktReader                         wktReader(s);
  srid_t                            srid = wktReader.readSRID();
  std::unique_ptr<Geometry>         g(wktReader.readGeometry());
  std::unique_ptr<PreparedGeometry> uptr(
      new PreparedGeometry(std::move(g), srid));
  return std::move(uptr);
}

///
///
///
std::unique_ptr<PreparedGeometry>
readEwkt(const std::string &s)
{
  std::istringstream        iss(s);
  WktReader                 wktReader(iss);
  srid_t                    srid = wktReader.readSRID();
  std::unique_ptr<Geometry> g(wktReader.readGeometry());
  return std::unique_ptr<PreparedGeometry>(
      new PreparedGeometry(std::move(g), srid));
}

///
///
///
std::unique_ptr<PreparedGeometry>
readEwkt(const char *str, size_t len)
{
  CharArrayBuffer           buf(str, str + len);
  std::istream              istr(&buf);
  WktReader                 wktReader(istr);
  srid_t                    srid = wktReader.readSRID();
  std::unique_ptr<Geometry> g(wktReader.readGeometry());
  return std::unique_ptr<PreparedGeometry>(
      new PreparedGeometry(std::move(g), srid));
}

} // namespace io
} // namespace SFCGAL
