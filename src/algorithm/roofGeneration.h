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
  RoofType      type           = RoofType::PITCHED;
  double        buildingHeight = 0.0;    ///< Building extrusion height (base level)
  double        roofHeight     = 3.0;    ///< Maximum roof height above building
  double        slopeAngle     = 30.0;   ///< Slope angle in degrees
  RidgePosition ridgePosition  = RidgePosition::INTERIOR;
  bool          addHips        = false;  ///< Add hip treatment for gable roofs
  bool          closeBase      = true;   ///< Include base polygon in output
  bool          generateSolid  = true;   ///< Generate closed Solid instead of PolyhedralSurface
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
                    double slopeAngle,
                    RidgePosition ridgePosition = RidgePosition::INTERIOR)
    -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Generate a pitched roof with building height and roof height.
 *
 * Creates a complete building with base extrusion and roof structure.
 * The building is extruded to buildingHeight, then the roof is added
 * with a maximum height of roofHeight above the building.
 *
 * @param footprint The building footprint polygon
 * @param ridgeLine The ridge line defining the roof orientation
 * @param buildingHeight Height to extrude the building base
 * @param roofHeight Maximum height of the roof above building
 * @param slopeAngle The roof slope angle in degrees (0-90)
 * @param ridgePosition Position of ridge relative to polygon
 * @return A Solid representing the complete building with roof
 * @pre footprint must be a valid polygon
 * @pre ridgeLine must be a valid LineString
 * @pre buildingHeight >= 0
 * @pre roofHeight > 0
 * @pre slopeAngle must be between 0 and 90 degrees
 */
SFCGAL_API auto
generatePitchedRoof(const Polygon &footprint, const LineString &ridgeLine,
                    double buildingHeight, double roofHeight, double slopeAngle,
                    RidgePosition ridgePosition = RidgePosition::INTERIOR)
    -> std::unique_ptr<Geometry>;

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
 * @brief Generate a gable roof with building height and roof height.
 *
 * Creates a complete building with base extrusion and gable roof structure.
 *
 * @param footprint The building footprint polygon
 * @param ridgeLine The ridge line defining the roof orientation
 * @param buildingHeight Height to extrude the building base
 * @param roofHeight Maximum height of the roof above building
 * @param slopeAngle The roof slope angle in degrees (0-90)
 * @param addHips Whether to add hip treatment at gable ends
 * @return A Solid representing the complete building with gable roof
 */
SFCGAL_API auto
generateGableRoof(const Polygon &footprint, const LineString &ridgeLine,
                  double buildingHeight, double roofHeight, double slopeAngle,
                  bool addHips = false)
    -> std::unique_ptr<Geometry>;

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
 * @brief Generate a gable roof using automatic ridge line from medial axis.
 *
 * Creates a roof with two symmetric slopes meeting at the automatically computed
 * ridge line from the polygon's medial axis.
 *
 * @param footprint The building footprint polygon
 * @param slopeAngle The roof slope angle in degrees (0-90)
 * @return A PolyhedralSurface representing the gable roof
 * @pre footprint must be a valid polygon
 * @pre slopeAngle must be between 0 and 90 degrees
 */
SFCGAL_API auto
generateGableRoof(const Polygon &footprint, double slopeAngle)
    -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Generate a gable roof using automatic ridge line with vertical faces option.
 *
 * Creates a roof with two symmetric slopes meeting at the automatically computed
 * ridge line from the polygon's medial axis, with optional vertical faces at ends.
 *
 * @param footprint The building footprint polygon
 * @param slopeAngle The roof slope angle in degrees (0-90)
 * @param addVerticalFaces Whether to add vertical faces at ridge line ends
 * @return A PolyhedralSurface representing the gable roof
 * @pre footprint must be a valid polygon
 * @pre slopeAngle must be between 0 and 90 degrees
 */
SFCGAL_API auto
generateGableRoof(const Polygon &footprint, double slopeAngle, bool addVerticalFaces = false)
    -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Generate a gable roof using automatic ridge line with specified roof height.
 *
 * Creates a roof with two symmetric slopes meeting at the automatically computed
 * ridge line from the polygon's medial axis, using a direct roof height.
 *
 * @param footprint The building footprint polygon
 * @param roofHeight The maximum height of the roof
 * @param addVerticalFaces Whether to add vertical faces at ridge line ends
 * @return A PolyhedralSurface representing the gable roof
 * @pre footprint must be a valid polygon
 * @pre roofHeight must be positive
 */
SFCGAL_API auto
generateGableRoofWithHeight(const Polygon &footprint, double roofHeight, bool addVerticalFaces = false)
    -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Generate a gable roof with building height and roof height using automatic ridge line.
 *
 * Creates a complete building with base extrusion and gable roof structure using
 * the polygon's medial axis as the ridge line.
 *
 * @param footprint The building footprint polygon
 * @param buildingHeight Height to extrude the building base
 * @param roofHeight Maximum height of the roof above building
 * @param slopeAngle The roof slope angle in degrees (0-90)
 * @param addHips Whether to add hip treatment at gable ends
 * @return A Solid representing the complete building with gable roof
 */
SFCGAL_API auto
generateGableRoof(const Polygon &footprint, double buildingHeight, double roofHeight,
                  double slopeAngle, bool addHips = false)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Generate a gable roof with building height and roof height using automatic ridge line with vertical faces option.
 *
 * Creates a complete building with base extrusion and gable roof structure using
 * the polygon's medial axis as the ridge line, with optional vertical faces at ends.
 *
 * @param footprint The building footprint polygon
 * @param buildingHeight Height to extrude the building base
 * @param roofHeight Maximum height of the roof above building
 * @param slopeAngle The roof slope angle in degrees (0-90)
 * @param addHips Whether to add hip treatment at gable ends
 * @param addVerticalFaces Whether to add vertical faces at ridge line ends
 * @return A Solid representing the complete building with gable roof
 */
SFCGAL_API auto
generateGableRoof(const Polygon &footprint, double buildingHeight, double roofHeight,
                  double slopeAngle, bool addHips, bool addVerticalFaces)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Generate a skillion roof from a polygon footprint and ridge line.
 *
 * Creates a simple single-slope roof, essentially an alias for pitched roof
 * with edge ridge position.
 *
 * @param footprint The building footprint polygon
 * @param ridgeLine The ridge line (typically on polygon edge)
 * @param slopeAngle The roof slope angle in degrees (0-90)
 * @return A PolyhedralSurface representing the skillion roof
 * @pre footprint must be a valid polygon
 * @pre ridgeLine must be a valid LineString
 * @pre slopeAngle must be between 0 and 90 degrees
 */
SFCGAL_API auto
generateSkillionRoof(const Polygon &footprint, const LineString &ridgeLine,
                     double slopeAngle)
    -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Generate a skillion roof with building height and roof height.
 *
 * Creates a complete building with base extrusion and skillion roof structure.
 *
 * @param footprint The building footprint polygon
 * @param ridgeLine The ridge line (typically on polygon edge)
 * @param buildingHeight Height to extrude the building base
 * @param roofHeight Maximum height of the roof above building
 * @param slopeAngle The roof slope angle in degrees (0-90)
 * @return A Solid representing the complete building with skillion roof
 */
SFCGAL_API auto
generateSkillionRoof(const Polygon &footprint, const LineString &ridgeLine,
                     double buildingHeight, double roofHeight, double slopeAngle)
    -> std::unique_ptr<Geometry>;

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
             const RoofParameters &params) -> std::unique_ptr<PolyhedralSurface>;

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
 * @param slopeAngle Slope angle in degrees
 * @return Height of the ridge point
 */
SFCGAL_API auto
calculateRidgeHeight(double horizontalDistance, double slopeAngle) -> double;

/**
 * @brief Calculate the horizontal distance for a given height and slope angle.
 *
 * Utility function to calculate the horizontal distance required to achieve
 * a specific height at a given slope angle.
 *
 * @param height Target height
 * @param slopeAngle Slope angle in degrees
 * @return Required horizontal distance
 */
SFCGAL_API auto
calculateHorizontalDistance(double height, double slopeAngle) -> double;

/**
 * @brief Generate a complete building with roof using unified parameters.
 *
 * Enhanced version of generateRoof that creates complete buildings with
 * separate building and roof heights, and can generate closed Solid geometry.
 *
 * @param footprint The building footprint polygon
 * @param ridgeLine The ridge line defining the roof orientation
 * @param params Roof generation parameters including building and roof heights
 * @return A Geometry (PolyhedralSurface or Solid) representing the complete building
 * @pre footprint must be a valid polygon
 * @pre ridgeLine must be a valid LineString (for non-flat roofs)
 */
SFCGAL_API auto
generateBuildingWithRoof(const Polygon &footprint, const LineString &ridgeLine,
                        const RoofParameters &params) -> std::unique_ptr<Geometry>;

/**
 * @brief Generate a complete building with roof using unified parameters without validity check.
 * @param footprint The building footprint polygon
 * @param ridgeLine The ridge line defining the roof orientation
 * @param params Roof generation parameters including building and roof heights
 * @param nvc NoValidityCheck object
 * @return A Geometry representing the complete building
 * @warning No actual validity check is conducted.
 */
SFCGAL_API auto
generateBuildingWithRoof(const Polygon &footprint, const LineString &ridgeLine,
                        const RoofParameters &params, NoValidityCheck &nvc)
    -> std::unique_ptr<Geometry>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_ROOFGENERATION_H_