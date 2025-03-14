// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/Surface.h"
#include "SFCGAL/GeometryVisitor.h"

namespace SFCGAL {

Surface::~Surface() = default;

auto
Surface::dimension() const -> int
{
  return 2;
}

Surface::Surface() = default;

Surface::Surface(Surface const &other) = default;

/**
 * [OGC/SFS]"The area of this Surface, as measured in the spatial reference
 * system of this Surface"
 */
// virtual double area() const = 0 ;
/**
 * [OGC/SFS]"The mathematical centroid for this Surface as a Point. The result
 * in not guaranteed to be on this Surface"
 */
// virtual Point centroid() const = 0 ;
/**
 * [OGC/SFS]"A Point guaranteed to be on this Surface"
 * @warning empty point is isEmpty()
 */
// virtual Point pointOnSurface() const = 0 ;

} // namespace SFCGAL
