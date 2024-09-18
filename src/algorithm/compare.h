// Copyright (c) 2024-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef _SFCGAL_ALGORITHM_COMPARE_H_
#define _SFCGAL_ALGORITHM_COMPARE_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/config.h"

namespace SFCGAL::algorithm::compare {

enum class FuzzyCompareMethod { Equal, DistanceEqual };

/**
 * Compares two geometries strictly
 * @ingroup public_api
 */
SFCGAL_API bool
strictCompare(const Geometry &ga, const Geometry &gb);

/**
 * Compares two geometries with sorting
 * @ingroup public_api
 */
SFCGAL_API bool
sortedCompare(const Geometry &ga, const Geometry &gb, bool useFuzzy = false,
              double             epsilon = 1e-8,
              FuzzyCompareMethod method  = FuzzyCompareMethod::Equal);

/**
 * Compares two geometries with a fuzzy tolerance
 * @ingroup public_api
 */
SFCGAL_API bool
fuzzyCompare(const Geometry &ga, const Geometry &gb, double epsilon = 1e-8,
             FuzzyCompareMethod method = FuzzyCompareMethod::Equal);

} // namespace SFCGAL::algorithm::compare

#endif
