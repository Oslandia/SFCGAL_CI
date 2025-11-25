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
 *
 * Input parameters:
 * - slopeAngle: Roof slope angle in degrees (0-90), used for gable and skillion
 * roofs
 * - roofHeight: Maximum roof height (>= 0), used for flat and hipped roofs
 * - buildingHeight: Height of building walls (>= 0), 0 = roof only mode
 * - closeBase: Whether to close the base (only if buildingHeight == 0)
 * - addVerticalFaces: Whether to add vertical gable/end faces (for gable and
 * skillion)
 *
 * Behavior:
 * - When buildingHeight == 0: generates roof only (PolyhedralSurface or Solid
 * if closed)
 * - When buildingHeight > 0: generates complete building with roof (always
 * Solid)
 * - closeBase is forced to false when buildingHeight > 0
 */
struct RoofParameters {
  RoofType type             = RoofType::GABLE;
  double   slopeAngle       = 30.0; ///< Slope angle in degrees (0-90) for
                                    ///< gable/skillion roofs
  double roofHeight = 3.0; ///< Maximum roof height for flat/hipped roofs (>= 0)
  double buildingHeight = 0.0;   ///< Building wall height (>= 0), 0 = roof only
  bool   closeBase      = false; ///< Close the base (only if buildingHeight ==
                                 ///< 0)
  bool addVerticalFaces =
      true; ///< Add vertical gable/end faces (for gable/skillion)
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
 * @brief Generate a skillion roof from a polygon footprint and ridge line.
 *
 * Creates a single-slope shed roof where the ridge line defines the high edge
 * and all other points slope down perpendicular to the ridge.
 *
 * @param footprint The building footprint polygon
 * @param ridgeLine The ridge line defining the high edge direction
 * @param slopeAngle The roof slope angle in degrees (0-90)
 * @param addVerticalFaces Whether to add vertical triangular faces at roof ends
 * @param buildingHeight Height of building walls (0 = roof only)
 * @param closeBase Whether to close the base (only used if buildingHeight == 0)
 * @return A PolyhedralSurface representing the skillion roof
 * @pre footprint must be a valid polygon
 * @pre ridgeLine must be a valid LineString
 * @pre slopeAngle must be between 0 and 90 degrees
 * @pre buildingHeight must be non-negative
 */
SFCGAL_API auto
generateSkillionRoof(const Polygon &footprint, const LineString &ridgeLine,
                     double slopeAngle, bool addVerticalFaces = true,
                     double buildingHeight = 0.0, bool closeBase = false)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Generate a roof using unified parameters.
 *
 * Main entry point for roof generation supporting all roof types.
 * - FLAT: Simple Z-translation or extrusion (slopeAngle, addVerticalFaces
 * ignored)
 * - HIPPED: Uses straight skeleton (slopeAngle, addVerticalFaces ignored)
 * - GABLE: Uses medial axis approach (slopeAngle, addVerticalFaces used)
 * - SKILLION: Uses ridge line (slopeAngle, addVerticalFaces used, requires
 * ridgeLine)
 *
 * @param footprint The building footprint polygon
 * @param ridgeLine The ridge line (required for SKILLION, ignored for others)
 * @param params Roof generation parameters
 * @return A Geometry representing the generated roof (PolyhedralSurface or
 * Solid)
 * @pre footprint must be a valid polygon
 * @pre ridgeLine must be a valid LineString for SKILLION type
 */
SFCGAL_API auto
generateRoof(const Polygon &footprint, const LineString &ridgeLine,
             const RoofParameters &params) -> std::unique_ptr<Geometry>;

/**
 * @brief Generate a roof using unified parameters without validity check.
 * @param footprint The building footprint polygon
 * @param ridgeLine The ridge line
 * @param params Roof generation parameters
 * @param nvc NoValidityCheck object
 * @return A Geometry representing the generated roof
 * @warning No actual validity check is conducted.
 */
SFCGAL_API auto
generateRoof(const Polygon &footprint, const LineString &ridgeLine,
             const RoofParameters &params, NoValidityCheck &nvc)
    -> std::unique_ptr<Geometry>;

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
 * medial axis of the polygon and using it as the ridge line.
 *
 * @param footprint The building footprint polygon
 * @param slopeAngle The roof slope angle in degrees (0-90)
 * @param addVerticalFaces Whether to add vertical gable faces at ridge endpoints
 * @param buildingHeight Height of the building walls (>= 0, 0 = roof only)
 * @param closeBase Whether to close the base (only used if buildingHeight == 0)
 * @return A Geometry representing the gable roof (PolyhedralSurface or Solid)
 * @pre footprint must be a valid polygon
 * @pre slopeAngle must be between 0 and 90 degrees
 * @pre buildingHeight must be non-negative
 */
SFCGAL_API auto
generateGableRoof(const Polygon &footprint, double slopeAngle,
                  bool addVerticalFaces = true, double buildingHeight = 0.0,
                  bool closeBase = false) -> std::unique_ptr<Geometry>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_ROOFGENERATION_H_