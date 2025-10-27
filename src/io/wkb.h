// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_IO_WKB_H_
#define SFCGAL_IO_WKB_H_

#include "SFCGAL/config.h"

#include <memory>
#include <sstream>
#include <string>

namespace SFCGAL {
class Geometry;
class PreparedGeometry;
} // namespace SFCGAL

namespace SFCGAL::io {

// WKB

/**
 * @brief Read a WKB geometry from an input stream
 * @param stream The input stream to read from
 * @param asHexString Whether the WKB data is encoded as hex string
 * @return A unique pointer to the geometry
 */
SFCGAL_API auto
readWkb(std::istream &stream, bool asHexString = false)
    -> std::unique_ptr<Geometry>;
/**
 * @brief Read a WKB geometry from a string
 * @param s The string containing WKB data
 * @param asHexString Whether the WKB data is encoded as hex string
 * @return A unique pointer to the geometry
 */
SFCGAL_API auto
readWkb(const std::string &s, bool asHexString = false)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Read a WKB geometry from a char array
 * @param str The character array containing WKB data
 * @param len The length of the character array
 * @param asHexString Whether the WKB data is encoded as hex string
 * @return A unique pointer to the geometry
 */
SFCGAL_API auto
readWkb(const char *str, size_t len, bool asHexString = false)
    -> std::unique_ptr<Geometry>;

// EWKB

/**
 * @brief Read a EWKB geometry from an input stream
 * @param stream The input stream to read from
 * @param asHexString Whether the EWKB data is encoded as hex string
 * @return A unique pointer to the prepared geometry
 */
SFCGAL_API auto
readEwkb(std::istream &stream, bool asHexString = false)
    -> std::unique_ptr<PreparedGeometry>;
/**
 * @brief Read a EWKB geometry from a string
 * @param s The string containing EWKB data
 * @param asHexString Whether the EWKB data is encoded as hex string
 * @return A unique pointer to the prepared geometry
 */
SFCGAL_API auto
readEwkb(const std::string &s, bool asHexString = false)
    -> std::unique_ptr<PreparedGeometry>;

/**
 * @brief Read a EWKB geometry from a char array
 * @param str The character array containing EWKB data
 * @param len The length of the character array
 * @param asHexString Whether the EWKB data is encoded as hex string
 * @return A unique pointer to the prepared geometry
 */
SFCGAL_API auto
readEwkb(const char *str, size_t len, bool asHexString = false)
    -> std::unique_ptr<PreparedGeometry>;

} // namespace SFCGAL::io

#endif
