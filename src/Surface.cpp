/**
 *   SFCGAL
 *
 *   Copyright (C) 2012-2013 Oslandia <infos@oslandia.com>
 *   Copyright (C) 2012-2013 IGN (http://www.ign.fr)
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Library General Public
 *   License as published by the Free Software Foundation; either
 *   version 2 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Library General Public License for more details.

 *   You should have received a copy of the GNU Library General Public
 *   License along with this library; if not, see
 <http://www.gnu.org/licenses/>.
 */

#include <SFCGAL/GeometryVisitor.h>
#include <SFCGAL/Surface.h>

namespace SFCGAL {

///
///
///
Surface::~Surface() {}

///
///
///
int
Surface::dimension() const
{
  return 2;
}

///
///
///
Surface::Surface() : Geometry() {}

///
///
///
Surface::Surface(Surface const &other) : Geometry(other) {}

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
