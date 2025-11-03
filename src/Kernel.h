// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_KERNEL_H_
#define SFCGAL_KERNEL_H_

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

namespace SFCGAL {

/**
 * default Kernel
 */

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

/// @{
/// @privatesection
/**
 * Quotient type
 * @private
 */
using QT = CGAL::Gmpq;

using Point_2              = CGAL::Point_2<SFCGAL::Kernel>;
using Triangle_2           = CGAL::Triangle_2<SFCGAL::Kernel>;
using Polygon_2            = CGAL::Polygon_2<SFCGAL::Kernel>;
using Vector_2             = CGAL::Vector_2<SFCGAL::Kernel>;
using Segment_2            = CGAL::Segment_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<SFCGAL::Kernel>;
using Polygon_set_2        = CGAL::Polygon_set_2<SFCGAL::Kernel>;

using Point_3        = CGAL::Point_3<SFCGAL::Kernel>;
using Triangle_3     = CGAL::Triangle_3<SFCGAL::Kernel>;
using Plane_3        = CGAL::Plane_3<SFCGAL::Kernel>;
using Segment_3      = CGAL::Segment_3<Kernel>;
using Vector_3       = CGAL::Vector_3<SFCGAL::Kernel>;
using Polyhedron_3   = CGAL::Polyhedron_3<Kernel>;
using Surface_mesh_3 = CGAL::Surface_mesh<Point_3>;

/// @} end of private section

} // namespace SFCGAL

#endif
