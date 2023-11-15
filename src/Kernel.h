// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef _SFCGAL_KERNEL_H_
#define _SFCGAL_KERNEL_H_

#if defined(_SFCGAL_EXACT_)
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#else
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#endif

namespace SFCGAL {

/**
 * default Kernel
 */

#if defined(_SFCGAL_EXACT_)
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
#else
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
#endif // defined(_SFCGAL_EXACT_)

/**
 * Quotient type
 */
typedef CGAL::Gmpq QT;

} // namespace SFCGAL

#endif
