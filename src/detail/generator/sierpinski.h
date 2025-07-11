// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_GENERATOR_SIERPINSKI_H_
#define SFCGAL_GENERATOR_SIERPINSKI_H_

#include "SFCGAL/config.h"

#include <memory>

namespace SFCGAL {
class MultiPolygon;
}

namespace SFCGAL {
namespace generator {

/**
 * generate hoch snowflake
 * @todo unittest
 */
SFCGAL_API std::unique_ptr<MultiPolygon>
           sierpinski(const unsigned int &order);

} // namespace generator
} // namespace SFCGAL

#endif
