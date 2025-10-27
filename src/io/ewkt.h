// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_IO_EWKT_H_
#define SFCGAL_IO_EWKT_H_

#include "SFCGAL/config.h"

#include <memory>
#include <sstream>
#include <string>

namespace SFCGAL {
class Geometry;
class PreparedGeometry;
} // namespace SFCGAL

namespace SFCGAL::io {
/**
 * @brief Read a EWKT prepared geometry from an input stream
 * @param s The input stream to read from
 * @return A unique pointer to the prepared geometry
 */
SFCGAL_API auto
readEwkt(std::istream &s) -> std::unique_ptr<PreparedGeometry>;
/**
 * @brief Read a EWKT geometry from a string
 * @param s The string containing EWKT data
 * @return A unique pointer to the prepared geometry
 */
SFCGAL_API auto
readEwkt(const std::string &s) -> std::unique_ptr<PreparedGeometry>;
/**
 * @brief Read a EWKT geometry from a char array
 * @param str The character array containing EWKT data
 * @param len The length of the character array
 * @return A unique pointer to the prepared geometry
 */
SFCGAL_API auto
readEwkt(const char *str, size_t len) -> std::unique_ptr<PreparedGeometry>;
} // namespace SFCGAL::io

#endif
