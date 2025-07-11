// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_DETAIL_COMPLEXCOMPARATOR_H_
#define SFCGAL_DETAIL_COMPLEXCOMPARATOR_H_

#include "SFCGAL/config.h"

#include <complex>

namespace SFCGAL {
namespace detail {

/**
 * lexicographic order on complex
 */
struct SFCGAL_API ComplexComparator{
    template <typename T> inline bool operator()(const std::complex<T> &a,
                                                 const std::complex<T> &b)
        const {return (a.real() < b.real()) ||
                      (a.real() == b.real() && a.imag() < b.imag());
} // namespace detail
}; // namespace SFCGAL

} // namespace detail
} // namespace SFCGAL

#endif
