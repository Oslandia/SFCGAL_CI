// Copyright (c) 2025-2025, SFCGAL team.
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
enum class RoofType : std::uint8_t {
  FLAT,     ///< Flat roof (simple extrusion)
  HIPPED,   ///< Hipped roof (straight skeleton extrusion)
  SKILLION, ///< Skillion roof (alias for pitched roof)
  GABLE     ///< Gable roof (dual symmetric slopes)
};

/**
 * @brief Parameters for roof generation algorithms.
 *
 * Input parameters:
 * - type: Roof type (cf. RoofType)
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
  RoofType type       = RoofType::GABLE; ///< Roof Type (cf. RoofType)
  double   slopeAngle = 30.0;            ///< Slope angle in degrees (0-90) for
                                         ///< gable/skillion roofs
  double roofHeight = 3.0; ///< Maximum roof height for flat/hipped roofs (>= 0)
  std::vector<std::vector<Kernel::FT>>
         angles; //< angles for each edge or each polygon ring for hipped roofs
  double buildingHeight = 0.0;   ///< Building wall height (>= 0), 0 = roof only
  bool   closeBase      = false; ///< Close the base (only if buildingHeight ==
                                 ///< 0)
  bool addVerticalFaces =
      true; ///< Add vertical gable/end faces (for gable/skillion)
};

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
 * @brief Generate a gable roof automatically using medial axis
 *
 * This function automatically generates a gable roof by computing the
 * medial axis of the polygon and using it as the ridge line.
 *
 * @param footprint The building footprint polygon
 * @param slopeAngle The roof slope angle in degrees (0-90)
 * @param addVerticalFaces Whether to add vertical gable faces at ridge
 * endpoints
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

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_ROOFGENERATION_H_
