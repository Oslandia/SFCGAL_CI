// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _SFCGAL_KERNEL_H_
#define _SFCGAL_KERNEL_H_

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

namespace SFCGAL {

/**
 * default Kernel
 */
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

/**
 * Quotient type
 */
typedef CGAL::Gmpq QT;

} // namespace SFCGAL

#endif
