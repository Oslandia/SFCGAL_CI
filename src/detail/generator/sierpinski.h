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
 * @brief generate sierpinski triangle
 * @param order The order of the sierpinski triangle
 * @return A unique pointer to the generated multipolygon
 * @todo unittest
 */
SFCGAL_API auto
sierpinski(const unsigned int &order) -> std::unique_ptr<MultiPolygon>;

} // namespace generator
} // namespace SFCGAL

#endif
