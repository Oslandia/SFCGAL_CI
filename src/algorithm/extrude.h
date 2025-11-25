// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_EXTRUDE_H_
#define SFCGAL_ALGORITHM_EXTRUDE_H_

// SFCGAL
#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/config.h"

namespace SFCGAL::algorithm {

// Class forward declarations.
struct NoValidityCheck;

/**
 * @brief Returns a Geometry equal to the specified Geometry,
 *   extruded by the specified displacement.
 * @param geometry The specified Geometry.
 * @param deltaX The component of the specified displacement in
 *   the x-direction.
 * @param deltaY The component of the specified displacement in
 *   the y-direction.
 * @param deltaZ The component of the specified displacement in
 *   the z-direction.
 * @return A Geometry equal to geometry extruded by the displacement
 *   vector {deltaX, deltaY, deltaZ}.
 * @pre geometry must be a valid geometry.
 * @pre deltaX, deltaY and deltaZ must all be finite.
 * @note If geometry is such that geometry.isMeasured() is true, then,
 *   since there is no common expectation of the
 *   values of the measures on the returned Geometry,
 *   all measures from the result are removed.
 * @note When applied to NURBSCurve geometries, the extrusion
 *   is internally performed on a LineString obtained via
 *   toLineString() with its default parameters.
 */
SFCGAL_API auto
extrude(const Geometry &geometry, const Kernel::FT &deltaX,
        const Kernel::FT &deltaY, const Kernel::FT &deltaZ)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Extrude geometry without validity check.
 * @param inputGeom Input geometry to extrude.
 * @param displacementX Component of displacement in x-direction.
 * @param displacementY Component of displacement in y-direction.
 * @param displacementZ Component of displacement in z-direction.
 * @param nvc A NoValidityCheck object.
 * @return Extruded geometry.
 * @warning No actual validity check is conducted.
 */
SFCGAL_API auto
extrude(const Geometry &inputGeom, const Kernel::FT &displacementX,
        const Kernel::FT &displacementY, const Kernel::FT &displacementZ,
        NoValidityCheck &nvc) -> std::unique_ptr<Geometry>;

/**
 * @brief Extrude geometry with double precision parameters.
 * @param geom Input geometry to extrude.
 * @param displacementX Component of displacement in x-direction.
 * @param displacementY Component of displacement in y-direction.
 * @param displacementZ Component of displacement in z-direction.
 * @return Extruded geometry.
 * @warning No actual validity check is conducted.
 */
SFCGAL_API auto
extrude(const Geometry &geom, const double &displacementX,
        const double &displacementY, const double &displacementZ)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Extrude geometry using a displacement vector.
 * @param inputGeometry The specified Geometry.
 * @param vector The displacement vector.
 * @return A Geometry equal to inputGeometry extruded by the displacement
 * vector.
 * @pre inputGeometry must be a valid geometry.
 * @note When applied to NURBSCurve geometries, the extrusion
 *   is internally performed on a LineString obtained via
 *   toLineString() with its default parameters.
 * @todo Improve extrude for 3D surfaces - Extrude only
 *   faces whose scalar_product(vector,normal) > 0 and use
 *   Polyhedron union to get output geometries with a clean
 *   topology.
 */
SFCGAL_API auto
extrude(const Geometry &inputGeometry, const Kernel::Vector_3 &vector)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Extrude polygon with specified height.
 * @param polygon The polygon to extrude.
 * @param height The extrusion height.
 * @return Extruded geometry.
 */
SFCGAL_API auto
extrude(const Polygon &polygon, const double &height)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Extrude a footprint polygon vertically until it intersects a roof
 * surface.
 *
 * This algorithm performs vertical ray-casting from each vertex of the
 * footprint (at z=0) to find the intersection with the roof PolyhedralSurface.
 * The resulting Solid has:
 * - Bottom face: the original footprint at z=0
 * - Top face: the footprint vertices projected to roof intersection heights
 * - Lateral faces: vertical walls connecting bottom and top edges
 *
 * @param footprint The base polygon at z=0 (can include interior rings/holes).
 * @param roof The roof surface as a PolyhedralSurface with variable z values.
 * @return A Solid representing the extruded volume, or an empty Solid if
 *         any vertex fails to intersect the roof (with error message).
 *
 * @pre footprint must be a valid 2D polygon (z coordinates ignored/set to 0).
 * @pre roof must be a valid PolyhedralSurface above the footprint (z > 0).
 *
 * @note If all projected vertices form a nearly planar surface (within 0.01
 *       tolerance), the top face is a single polygon. Otherwise, the top
 *       surface is triangulated to handle non-planar cases.
 * @note Face orientations follow the right-hand rule: normals point outward
 *       from the Solid.
 * @note Uses SFCGAL::EPSILON (1e-8) for intersection tolerance.
 */
SFCGAL_API auto
extrudeUntil(const Polygon &footprint, const PolyhedralSurface &roof)
    -> std::unique_ptr<Solid>;

} // namespace SFCGAL::algorithm

#endif // ! SFCGAL_ALGORITHM_EXTRUDE_H_
