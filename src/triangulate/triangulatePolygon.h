// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TRIANGULATE_TRIANGULATEPOLYGON_H_
#define SFCGAL_TRIANGULATE_TRIANGULATEPOLYGON_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Geometry.h"

namespace SFCGAL::triangulate {

/**
 * @brief Triangulate 3D polygons in a Geometry.
 *
 * @param g input geometry
 * @param triangulatedSurface resulting TriangulatedSurface
 */
SFCGAL_API auto
triangulatePolygon3D(const Geometry      &g,
                     TriangulatedSurface &triangulatedSurface) -> void;
/**
 * @brief Triangulate a 3D Polygon
 * @param g Input polygon to triangulate
 * @param triangulatedSurface Output triangulated surface
 * @todo unittest
 */
SFCGAL_API void
triangulatePolygon3D(const Polygon       &g,
                     TriangulatedSurface &triangulatedSurface);
/**
 * @brief Triangulate a 3D Triangle (copy triangle)
 * @param g Input triangle to triangulate
 * @param triangulatedSurface Output triangulated surface
 * @todo unittest
 */
SFCGAL_API void
triangulatePolygon3D(const Triangle      &g,
                     TriangulatedSurface &triangulatedSurface);
/**
 * @brief Triangulate a 3D TriangulatedSurface (copy triangles)
 * @param g Input triangulated surface to copy
 * @param triangulatedSurface Output triangulated surface
 * @todo unittest
 */
SFCGAL_API void
triangulatePolygon3D(const TriangulatedSurface &g,
                     TriangulatedSurface       &triangulatedSurface);
/**
 * @brief Triangulate a 3D MultiPolygon
 * @param g Input geometry collection to triangulate
 * @param triangulatedSurface Output triangulated surface
 * @todo unittest
 */
SFCGAL_API void
opentriangulatePolygon3D(const GeometryCollection &g,
                         TriangulatedSurface      &triangulatedSurface);
/**
 * @brief Triangulate 3D polygons in a PolyhedralSurface.
 * @param polyhedralSurface Input polyhedral surface to triangulate
 * @param triangulatedSurface Output triangulated surface
 */
SFCGAL_API void
triangulatePolygon3D(const PolyhedralSurface &polyhedralSurface,
                     TriangulatedSurface     &triangulatedSurface);
/**
 * @brief Triangulate a Solid
 * @param g Input solid to triangulate
 * @param triangulatedSurface Output triangulated surface
 * @todo unittest
 */
SFCGAL_API void
triangulatePolygon3D(const Solid &g, TriangulatedSurface &triangulatedSurface);

} // namespace SFCGAL::triangulate

#endif
