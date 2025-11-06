// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_ROOFGENERATION_H_
#define SFCGAL_ALGORITHM_ROOFGENERATION_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/config.h"

namespace SFCGAL {
class Point;
}

namespace SFCGAL::algorithm {

struct NoValidityCheck;

/**
 * @brief Enumeration of roof types supported by the roof generation algorithms.
 */
enum class RoofType {
  FLAT,     ///< Flat roof (simple extrusion)
  HIPPED,   ///< Hipped roof (straight skeleton extrusion)
  PITCHED,  ///< Pitched roof (single slope from ridge line)
  SKILLION, ///< Skillion roof (alias for pitched roof)
  GABLE     ///< Gable roof (dual symmetric slopes)
};

/**
 * @brief Position of the ridge line relative to the polygon footprint.
 */
enum class RidgePosition {
  INTERIOR, ///< Ridge line runs through the polygon interior (most common)
  EDGE,     ///< Ridge line coincides with polygon edge (shed roof)
  EXTERIOR  ///< Ridge line extends outside polygon boundaries
};

/**
 * @brief Parameters for roof generation algorithms.
 */
struct RoofParameters {
  RoofType      type          = RoofType::PITCHED;
  double        height        = 3.0;  ///< Maximum roof height
  double        slopeAngle    = 30.0; ///< Slope angle in degrees
  RidgePosition ridgePosition = RidgePosition::INTERIOR;
  bool          addHips       = false; ///< Add hip treatment for gable roofs
  bool          closeBase     = true;  ///< Include base polygon in output
};

/**
 * @brief Generate a pitched roof from a polygon footprint and ridge line.
 *
 * Creates a roof with a single sloped plane or two planes meeting at the ridge
 * line. The ridge line can be positioned inside, on the edge, or outside the
 * polygon.
 *
 * @param footprint The building footprint polygon
 * @param ridgeLine The ridge line defining the roof orientation
 * @param slopeAngle The roof slope angle in degrees (0-90)
 * @param ridgePosition Position of ridge relative to polygon
 * @return A PolyhedralSurface representing the pitched roof
 * @pre footprint must be a valid polygon
 * @pre ridgeLine must be a valid LineString
 * @pre slopeAngle must be between 0 and 90 degrees
 */
SFCGAL_API auto
generatePitchedRoof(const Polygon &footprint, const LineString &ridgeLine,
                    double        slopeAngle,
                    RidgePosition ridgePosition = RidgePosition::INTERIOR)
    -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Generate a pitched roof without validity check.
 * @param footprint The building footprint polygon
 * @param ridgeLine The ridge line defining the roof orientation
 * @param slopeAngle The roof slope angle in degrees
 * @param ridgePosition Position of ridge relative to polygon
 * @param nvc NoValidityCheck object
 * @return A PolyhedralSurface representing the pitched roof
 * @warning No actual validity check is conducted.
 */
SFCGAL_API auto
generatePitchedRoof(const Polygon &footprint, const LineString &ridgeLine,
                    double slopeAngle, RidgePosition ridgePosition,
                    NoValidityCheck &nvc) -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Generate a gable roof from a polygon footprint and ridge line.
 *
 * Creates a roof with two symmetric slopes meeting at the ridge line, with
 * triangular gable ends. Optionally adds hip treatment at the ends.
 *
 * @param footprint The building footprint polygon
 * @param ridgeLine The ridge line defining the roof orientation
 * @param slopeAngle The roof slope angle in degrees (0-90)
 * @param addHips Whether to add hip treatment at gable ends
 * @return A PolyhedralSurface representing the gable roof
 * @pre footprint must be a valid polygon
 * @pre ridgeLine must be a valid LineString
 * @pre slopeAngle must be between 0 and 90 degrees
 */
SFCGAL_API auto
generateGableRoof(const Polygon &footprint, const LineString &ridgeLine,
                  double slopeAngle, bool addHips = false)
    -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Generate a gable roof without validity check.
 * @param footprint The building footprint polygon
 * @param ridgeLine The ridge line defining the roof orientation
 * @param slopeAngle The roof slope angle in degrees
 * @param addHips Whether to add hip treatment at gable ends
 * @param nvc NoValidityCheck object
 * @return A PolyhedralSurface representing the gable roof
 * @warning No actual validity check is conducted.
 */
SFCGAL_API auto
generateGableRoof(const Polygon &footprint, const LineString &ridgeLine,
                  double slopeAngle, bool addHips, NoValidityCheck &nvc)
    -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Generate a skillion roof from a polygon footprint and ridge line.
 *
 * Creates a single-slope shed roof where the ridge line defines the high edge
 * and all other points slope down perpendicular to the ridge. This creates
 * the characteristic mono-pitch roof with consistent slope across the entire
 * surface.
 *
 * @param footprint The building footprint polygon
 * @param ridgeLine The ridge line defining the high edge direction
 * @param slopeAngle The roof slope angle in degrees (0-90)
 * @param addVerticalFaces Whether to add vertical triangular faces at roof ends
 * @return A PolyhedralSurface representing the skillion roof
 * @pre footprint must be a valid polygon
 * @pre ridgeLine must be a valid LineString
 * @pre slopeAngle must be between 0 and 90 degrees
 */
SFCGAL_API auto
generateSkillionRoof(const Polygon &footprint, const LineString &ridgeLine,
                     double slopeAngle, bool addVerticalFaces = false)
    -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Generate a skillion roof with optional building integration.
 *
 * Creates a single-slope shed roof, optionally combined with building walls.
 * When buildingHeight > 0, generates both the building structure and roof.
 *
 * @param footprint The building footprint polygon
 * @param ridgeLine The ridge line defining the high edge direction
 * @param slopeAngle The roof slope angle in degrees (0-90)
 * @param addVerticalFaces Whether to add vertical triangular faces at roof ends
 * @param buildingHeight Height of building walls (0 = roof only)
 * @return A PolyhedralSurface representing the skillion roof or complete
 * building
 * @pre footprint must be a valid polygon
 * @pre ridgeLine must be a valid LineString
 * @pre slopeAngle must be between 0 and 90 degrees
 * @pre buildingHeight must be non-negative
 */
SFCGAL_API auto
generateSkillionRoof(const Polygon &footprint, const LineString &ridgeLine,
                     double slopeAngle, double buildingHeight)
    -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Generate a skillion roof with optional vertical faces and building
 * integration.
 *
 * Creates a single-slope shed roof with full control over vertical faces and
 * building walls. When buildingHeight > 0, generates both the building
 * structure and roof.
 *
 * @param footprint The building footprint polygon
 * @param ridgeLine The ridge line defining the high edge direction
 * @param slopeAngle The roof slope angle in degrees (0-90)
 * @param addVerticalFaces Whether to add vertical triangular faces at roof ends
 * @param buildingHeight Height of building walls (0 = roof only)
 * @return A PolyhedralSurface representing the skillion roof or complete
 * building
 * @pre footprint must be a valid polygon
 * @pre ridgeLine must be a valid LineString
 * @pre slopeAngle must be between 0 and 90 degrees
 * @pre buildingHeight must be non-negative
 */
SFCGAL_API auto
generateSkillionRoof(const Polygon &footprint, const LineString &ridgeLine,
                     double slopeAngle, bool addVerticalFaces,
                     double buildingHeight)
    -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Generate a roof using unified parameters.
 *
 * Main entry point for roof generation supporting all roof types through
 * a unified parameter structure.
 *
 * @param footprint The building footprint polygon
 * @param ridgeLine The ridge line defining the roof orientation
 * @param params Roof generation parameters
 * @return A PolyhedralSurface representing the generated roof
 * @pre footprint must be a valid polygon
 * @pre ridgeLine must be a valid LineString (for non-flat roofs)
 */
SFCGAL_API auto
generateRoof(const Polygon &footprint, const LineString &ridgeLine,
             const RoofParameters &params)
    -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Generate a roof using unified parameters without validity check.
 * @param footprint The building footprint polygon
 * @param ridgeLine The ridge line defining the roof orientation
 * @param params Roof generation parameters
 * @param nvc NoValidityCheck object
 * @return A PolyhedralSurface representing the generated roof
 * @warning No actual validity check is conducted.
 */
SFCGAL_API auto
generateRoof(const Polygon &footprint, const LineString &ridgeLine,
             const RoofParameters &params, NoValidityCheck &nvc)
    -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Calculate the ridge height for a given slope angle and distance.
 *
 * Utility function to calculate the height of a ridge point given the
 * horizontal distance from the base and the slope angle.
 *
 * @param horizontalDistance Distance from base to ridge in horizontal plane
 * (must be non-negative)
 * @param slopeAngle Slope angle in degrees (0 < angle < 90)
 * @return Height of the ridge point using exact arithmetic
 *
 * @pre horizontalDistance >= 0
 * @pre 0 < slopeAngle < 90
 *
 * @example
 * ```cpp
 * double height = calculateRidgeHeight(3.0, 30.0); // height ≈ 1.73 for 30°
 * slope
 * ```
 */
SFCGAL_API auto
calculateRidgeHeight(double horizontalDistance, double slopeAngle) -> double;

/**
 * @brief Calculate the horizontal distance for a given height and slope angle.
 *
 * Utility function to calculate the horizontal distance required to achieve
 * a specific height at a given slope angle.
 *
 * @param height Target height (must be non-negative)
 * @param slopeAngle Slope angle in degrees (0 < angle < 90)
 * @return Required horizontal distance using exact arithmetic
 *
 * @pre height >= 0
 * @pre 0 < slopeAngle < 90
 *
 * @example
 * ```cpp
 * double distance = calculateHorizontalDistance(2.0, 45.0); // distance = 2.0
 * for 45° slope
 * ```
 */
SFCGAL_API auto
calculateHorizontalDistance(double height, double slopeAngle) -> double;

/**
 * @brief Generate a gable roof automatically using medial axis
 *
 * This function automatically generates a gable roof by computing the
 * medial axis of the polygon and using it as the ridge line. It can optionally
 * combine with building extrusion to create a complete building with gable
 * roof.
 *
 * @param footprint The building footprint polygon
 * @param slopeAngle The roof slope angle in degrees
 * @param addVerticalFaces Whether to add vertical triangular faces at ridge
 * endpoints
 * @param buildingHeight Height of the building walls (0 = roof only)
 * @return A PolyhedralSurface representing the gable roof or complete building
 * @pre footprint must be a valid polygon
 * @pre slopeAngle must be between 0 and 90 degrees
 * @pre buildingHeight must be non-negative
 */
SFCGAL_API auto
generateGableRoof(const Polygon &footprint, double slopeAngle,
                  bool addVerticalFaces = false, double buildingHeight = 0.0)
    -> std::unique_ptr<PolyhedralSurface>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_ROOFGENERATION_H_