// Copyright (c) 2025-2025, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#pragma once

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>

#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <vector>

namespace SFCGAL {
namespace detail {

/**
 * @brief Simplifies a LineString
 */
auto
simplifyLineString(const LineString &lineString, double threshold,
                   bool preserveTopology, const SegmentStore &store)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Simplifies a MultiLineString
 */
auto
simplifyMultiLineString(const MultiLineString &multiLine, double threshold,
                        bool preserveTopology) -> std::unique_ptr<Geometry>;

/**
 * @brief Simplifies a Polygon
 */
auto
simplifyPolygon(const Polygon &polygon, double threshold, bool preserveTopology)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Simplifies a MultiPolygon
 */
auto
simplifyMultiPolygon(const MultiPolygon &multiPolygon, double threshold,
                     bool preserveTopology) -> std::unique_ptr<Geometry>;

/**
 * @brief Simplifies a GeometryCollection with topology preservation
 */
auto
simplifyGeometryCollectionTopology(const GeometryCollection &collection,
                                   double                    threshold)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Simplifies a PolyhedralSurface
 */
auto
simplifyPolyhedralSurface(const PolyhedralSurface &polySurface,
                          double threshold, bool preserveTopology)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Simplifies a GeometryCollection
 */
auto
simplifyGeometryCollection(const GeometryCollection &collection,
                           double threshold, bool preserveTopology)
    -> std::unique_ptr<Geometry>;

} // namespace detail
} // namespace SFCGAL
