// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_GEOMETRY_STREAMS_H_
#define SFCGAL_GEOMETRY_STREAMS_H_

#include <ostream>

#include "SFCGAL/config.h"

namespace SFCGAL {

class Envelope;
class Geometry;

/**
 * @brief Ostream operator for Envelope
 * @param os Output stream
 * @param envelope Envelope to output
 * @return Reference to the output stream
 */
SFCGAL_API auto
operator<<(std::ostream &os, const Envelope &envelope) -> std::ostream &;

/**
 * @brief Ostream operator for Geometry
 * @param os Output stream
 * @param geometry Geometry to output
 * @return Reference to the output stream
 */
SFCGAL_API auto
operator<<(std::ostream &os, const Geometry &geometry) -> std::ostream &;
} // namespace SFCGAL

#endif
