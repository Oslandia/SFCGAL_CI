// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _SFCGAL_IO_WKT_H_
#define _SFCGAL_IO_WKT_H_

#include <SFCGAL/config.h>

#include <memory>
#include <sstream>
#include <string>

namespace SFCGAL {
class Geometry;
}

namespace SFCGAL {
namespace io {
/**
 * Read a WKT geometry from an input stream
 */
SFCGAL_API std::unique_ptr<Geometry>
           readWkt(std::istream &s);
/**
 * Read a WKT geometry from a string
 */
SFCGAL_API std::unique_ptr<Geometry>
           readWkt(const std::string &s);
/**
 * Read a WKT geometry from a char*
 */
SFCGAL_API std::unique_ptr<Geometry>
           readWkt(const char *, size_t);
} // namespace io
} // namespace SFCGAL

#endif
