// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_UBLAS_H_
#define SFCGAL_UBLAS_H_

// Boost.uBLAS uses deprecated std::iterator (removed in C++17)
// Suppress warnings until Boost fixes this issue
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wdeprecated-declarations"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#elif defined(_MSC_VER)
  #pragma warning(push)
  #pragma warning(disable : 4996)
#endif

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#elif defined(_MSC_VER)
  #pragma warning(pop)
#endif

namespace SFCGAL::detail {
namespace ublas = boost::numeric::ublas;
} // namespace SFCGAL::detail

#endif
