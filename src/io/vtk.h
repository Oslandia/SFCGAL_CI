// Copyright (c) 2024-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later
#ifndef SFCGAL_IO_VTK_H_
#define SFCGAL_IO_VTK_H_

#include "SFCGAL/Geometry.h"
#include <ostream>
#include <string>

namespace SFCGAL {
namespace io {
namespace VTK {

/**
 * @brief Saves a geometry to a legacy VTK format stream.
 *
 * @param[in] geom The geometry to save
 * @param[out] out The output stream
 * @throws SFCGAL::Exception If the geometry is invalid or unsupported
 */
SFCGAL_API void
save(const Geometry &geom, std::ostream &out);

/**
 * @brief Saves a geometry to a legacy VTK file.
 *
 * @param[in] geom The geometry to save
 * @param[in] filename The name of the file to save to
 * @throws SFCGAL::Exception If the file cannot be opened or the geometry is
 * invalid
 */
SFCGAL_API void
save(const Geometry &geom, const std::string &filename);

/**
 * @brief Saves a geometry to a legacy VTK format string.
 *
 * @param[in] geom The geometry to save
 * @return The legacy VTK format string
 * @throws SFCGAL::Exception If the geometry is invalid or unsupported
 */
SFCGAL_API std::string
           saveToString(const Geometry &geom);

/**
 * @brief Saves a geometry to a legacy VTK format buffer (C API).
 *
 * @param[in] geom The geometry to save
 * @param[out] buffer The buffer to write to
 * @param[in,out] size On input, the size of the buffer. On output, the number
 * of bytes written (or required if buffer is null)
 * @throws SFCGAL::Exception If the geometry is invalid or unsupported
 */
SFCGAL_API void
saveToBuffer(const Geometry &geom, char *buffer, size_t *size);

} // namespace VTK
} // namespace io
} // namespace SFCGAL

#endif // SFCGAL_IO_VTK_H_
