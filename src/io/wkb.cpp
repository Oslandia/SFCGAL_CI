// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2023, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/io/wkb.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/detail/io/WkbReader.h"
#include "SFCGAL/detail/io/WkbWriter.h"
#include "SFCGAL/detail/tools/CharArrayBuffer.h"

using namespace SFCGAL::detail::io;

namespace SFCGAL::io {

///
///
///
auto
readWkb(std::istream &stream, bool asHexString) -> std::unique_ptr<Geometry>
{
  WkbReader wkbReader(stream, asHexString);
  wkbReader.readWkb();
  return wkbReader.geometry();
}

///
///
///
auto
readWkb(const std::string &s, bool asHexString) -> std::unique_ptr<Geometry>
{
  std::istringstream iss(s);
  WkbReader          wkbReader(iss, asHexString);
  wkbReader.readWkb();
  return wkbReader.geometry();
}

///
///
///
auto
readWkb(const char *str, size_t len, bool asHexString)
    -> std::unique_ptr<Geometry>
{
  CharArrayBuffer buf(str, str + len);
  std::istream    istr(&buf);

  return readWkb(istr, asHexString);
}

/**
 * Read a WKB geometry from an input stream
 */
auto
readEwkb(std::istream &stream, bool asHexString)
    -> std::unique_ptr<PreparedGeometry>
{

  WkbReader wkbReader(stream, asHexString);
  wkbReader.readWkb();
  return wkbReader.preparedGeometry();
}

auto
readEwkb(const std::string &s, bool asHexString)
    -> std::unique_ptr<PreparedGeometry>
{
  std::istringstream iss(s);
  WkbReader          wkbReader(iss, asHexString);
  wkbReader.readWkb();
  return wkbReader.preparedGeometry();
}

/**
 * Read a WKB geometry from a char*
 */
auto
readEwkb(const char *str, size_t len, bool asHexString)
    -> std::unique_ptr<PreparedGeometry>
{
  CharArrayBuffer buf(str, str + len);
  std::istream    istr(&buf);

  return readEwkb(istr, asHexString);
}

} // namespace SFCGAL::io
