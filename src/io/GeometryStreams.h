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
 * Ostream operator for Envelope;
 */
SFCGAL_API std::ostream            &
operator<<(std::ostream &, const Envelope &);

/**
 * Ostream operator for Geometry;
 */
SFCGAL_API std::ostream            &
operator<<(std::ostream &, const Geometry &);
} // namespace SFCGAL

#endif
