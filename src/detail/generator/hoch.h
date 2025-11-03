// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_GENERATOR_HOCH_H_
#define SFCGAL_GENERATOR_HOCH_H_

#include "SFCGAL/config.h"

#include <memory>

namespace SFCGAL {
class Polygon;
}

namespace SFCGAL::generator {

/**
 * @brief generate hoch snowflake
 * @param order The order of the snowflake
 * @return A unique pointer to the generated polygon
 * @todo unittest
 */
SFCGAL_API auto
hoch(const unsigned int &order) -> std::unique_ptr<Polygon>;

} // namespace SFCGAL::generator

#endif
