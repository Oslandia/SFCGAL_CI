// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_IO_STL_H_
#define SFCGAL_IO_STL_H_

#include <SFCGAL/Geometry.h>
#include <ostream>
#include <string>

namespace SFCGAL::io::STL {

/**
 * @brief Saves a geometry to an STL format stream.
 * @ingroup io
 *
 * @param[in] geom The geometry to save
 * @param[out] out The output stream
 * @throws std::runtime_error If the geometry is invalid or unsupported
 */
SFCGAL_API auto
save(const Geometry &geom, std::ostream &out) -> void;

/**
 * @brief Saves a geometry to an STL file.
 * @ingroup io
 *
 * @param[in] geom The geometry to save
 * @param[in] filename The name of the file to save to
 * @throws std::runtime_error If the file cannot be opened or the geometry is
 * invalid
 */
SFCGAL_API auto
save(const Geometry &geom, const std::string &filename) -> void;

/**
 * @brief Saves a geometry to an STL format string.
 * @ingroup io
 *
 * @param[in] geom The geometry to save
 * @return The STL format string
 * @throws std::runtime_error If the geometry is invalid or unsupported
 */
SFCGAL_API auto
saveToString(const Geometry &geom) -> std::string;

/**
 * @brief Saves a geometry to an STL format buffer (C API).
 * @ingroup io
 *
 * @param[in] geom The geometry to save
 * @param[out] buffer The buffer to write to
 * @param[in,out] size On input, the size of the buffer. On output, the number
 * of bytes written (or required if buffer is null)
 * @throws std::runtime_error If the geometry is invalid or unsupported
 */
SFCGAL_API auto
saveToBuffer(const Geometry &geom, char *buffer, size_t *size) -> void;

} // namespace SFCGAL::io::STL

#endif // SFCGAL_IO_STL_H_
