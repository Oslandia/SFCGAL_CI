// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_IO_WKT_H_
#define SFCGAL_IO_WKT_H_

#include "SFCGAL/config.h"

#include <memory>
#include <sstream>
#include <string>

namespace SFCGAL {
class Geometry;
}

namespace SFCGAL::io {
/**
 * @brief Read a WKT geometry from an input stream
 * @param s The input stream to read from
 * @return A unique pointer to the geometry
 */
SFCGAL_API auto
readWkt(std::istream &s) -> std::unique_ptr<Geometry>;
/**
 * @brief Read a WKT geometry from a string
 * @param s The string containing WKT data
 * @return A unique pointer to the geometry
 */
SFCGAL_API auto
readWkt(const std::string &s) -> std::unique_ptr<Geometry>;
/**
 * @brief Read a WKT geometry from a char array
 * @param str The character array containing WKT data
 * @param len The length of the character array
 * @return A unique pointer to the geometry
 */
SFCGAL_API auto
readWkt(const char *str, size_t len) -> std::unique_ptr<Geometry>;
} // namespace SFCGAL::io

#endif
