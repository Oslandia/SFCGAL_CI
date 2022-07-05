// Copyright (c) 2023, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef _SFCGAL_TRIANGULATE_TRIANGULATEPOLYGONALDOMAIN_H_
#define _SFCGAL_TRIANGULATE_TRIANGULATEPOLYGONALDOMAIN_H_

#include <SFCGAL/config.h>

#include <SFCGAL/Geometry.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/TriangulatedSurface.h>
#include <SFCGAL/detail/triangulate/ConstraintDelaunayTriangulation.h>

#include <SFCGAL/Kernel.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/Surface.h>
namespace SFCGAL {
namespace triangulate {
/**
 * @brief Constraint Delaunay Triangulation for a polygonal domain
 * with constrained geometries.
 */
SFCGAL_API auto
constrainedPolygonalDomain(const Polygon         &poly,
                           const MultiLineString &multiline,
                           const MultiPoint &multipoint) -> TriangulatedSurface;

} // namespace triangulate
} // namespace SFCGAL

#endif
