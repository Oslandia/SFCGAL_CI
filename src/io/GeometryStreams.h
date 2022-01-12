// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _SFCGAL_GEOMETRY_STREAMS_H_
#define _SFCGAL_GEOMETRY_STREAMS_H_

#include <ostream>

#include <SFCGAL/config.h>

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
