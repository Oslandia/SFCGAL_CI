// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _SFCGAL_GENERATOR_SIERPINSKI_H_
#define _SFCGAL_GENERATOR_SIERPINSKI_H_

#include <SFCGAL/config.h>

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
