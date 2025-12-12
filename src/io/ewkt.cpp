// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/io/ewkt.h"

#include "SFCGAL/detail/io/WktReader.h"
#include "SFCGAL/detail/io/WktWriter.h"
#include "SFCGAL/detail/tools/CharArrayBuffer.h"

#include <memory>

using namespace SFCGAL::detail::io;
using namespace SFCGAL::tools;

namespace SFCGAL::io {

auto
readEwkt(std::istream &s) -> std::unique_ptr<PreparedGeometry>
{
  WktReader    wktReader(s);
  srid_t const srid = wktReader.readSRID();
  auto         g    = wktReader.readGeometry();
  return std::make_unique<PreparedGeometry>(std::move(g), srid);
}

auto
readEwkt(const std::string &s) -> std::unique_ptr<PreparedGeometry>
{
  std::istringstream iss(s);
  WktReader          wktReader(iss);
  srid_t const       srid = wktReader.readSRID();
  auto               g    = wktReader.readGeometry();
  return std::make_unique<PreparedGeometry>(std::move(g), srid);
}

auto
readEwkt(const char *str, size_t len) -> std::unique_ptr<PreparedGeometry>
{
  CharArrayBuffer buf(str, str + len);
  std::istream    istr(&buf);
  WktReader       wktReader(istr);
  srid_t const    srid = wktReader.readSRID();
  auto            g    = wktReader.readGeometry();
  return std::make_unique<PreparedGeometry>(std::move(g), srid);
}

} // namespace SFCGAL::io
