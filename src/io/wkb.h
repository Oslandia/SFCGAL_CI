// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2023, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef _SFCGAL_IO_WKB_H_
#define _SFCGAL_IO_WKB_H_

#include <SFCGAL/config.h>

#include <memory>
#include <sstream>
#include <string>

namespace SFCGAL {
class Geometry;
class PreparedGeometry;
} // namespace SFCGAL

namespace SFCGAL {
namespace io {
/**
 * Read a WKB geometry from an input stream
 */
// SFCGAL_API std::unique_ptr<Geometry>
//            readWkb(std::istream &s);
/**
 * Read a WKB geometry from a string
 */
SFCGAL_API std::unique_ptr<Geometry>
           readWkb(const std::string &s);

SFCGAL_API auto
readEwkb(const std::string &s) -> std::unique_ptr<PreparedGeometry>;
/**
 * Read a WKB geometry from a char*
 */
// SFCGAL_API std::unique_ptr<Geometry>
//            readWkb(const char *, size_t);
} // namespace io
} // namespace SFCGAL

#endif
