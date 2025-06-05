// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/simplification.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/detail/ConstraintInfo.h"
#include "SFCGAL/detail/SegmentStore.h"

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

namespace PS = CGAL::Polyline_simplification_2;

using K   = CGAL::Projection_traits_xy_3<SFCGAL::Kernel>;
using Vb  = PS::Vertex_base_2<K>;
using Fb  = CGAL::Constrained_triangulation_face_base_2<K>;
using TDS = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using CDT =
    CGAL::Constrained_Delaunay_triangulation_2<K, TDS,
                                               CGAL::Exact_predicates_tag>;
using CT                  = CGAL::Constrained_triangulation_plus_2<CDT>;
using CGALPoint           = CT::Point;
using Constraint_id       = CT::Constraint_id;
using Constraint_iterator = CT::Constraint_iterator;
using Points_in_constraint_iterator = CT::Points_in_constraint_iterator;
using Stop                          = PS::Stop_above_cost_threshold;
using Cost                          = PS::Squared_distance_cost;

// Type alias for our constraint order info with the specific Constraint_id type
using ConstraintInfoType = ConstraintOrderInfo<Constraint_id>;

/**
 * @brief Extracts points from a LineString to create a CGAL point vector
 * @param lineString The source LineString
 * @return Vector of CGALPoints
 */
static std::vector<CGALPoint>
extractPointsFromLineString(const LineString &lineString)
{
  std::vector<CGALPoint> points;
  points.reserve(lineString.numPoints());

  for (std::size_t i = 0; i < lineString.numPoints(); ++i) {
    const auto &pt = lineString.pointN(i);
    points.emplace_back(CGAL::to_double(pt.x()), CGAL::to_double(pt.y()),
                        pt.is3D() ? CGAL::to_double(pt.z()) : 0.0);
  }

  return points;
}

/**
 * @brief Extracts points from a ring to create a CGAL point vector
 * @param ring The source ring (LineString)
 * @return Vector of CGALPoints
 */
static std::vector<CGALPoint>
extractPointsFromRing(const LineString &ring)
{
  return extractPointsFromLineString(ring);
}

/**
 * @brief Builds a LineString from constraint points
 * @param ct The constrained triangulation
 * @param cid The constraint ID
 * @param store The segment store for interpolation
 * @param dimension The coordinate dimension type
 * @return A LineString built from the simplified constraint
 */
static LineString
buildLineStringFromConstraint(const CT &ct, const Constraint_id &cid,
                              const SegmentStore &store,
                              CoordinateType      dimension)
{
  LineString result;

  for (auto vit = ct.points_in_constraint_begin(cid);
       vit != ct.points_in_constraint_end(cid); ++vit) {
    double x = CGAL::to_double(vit->x());
    double y = CGAL::to_double(vit->y());

    result.addPoint(store.createPoint(x, y, dimension));
  }

  return result;
}

/**
 * @brief Extracts LineString constraints and adds them to the triangulation
 * @param lineString The source LineString
 * @param ct The constrained triangulation to add to
 * @param constraintInfos Collection of constraint infos
 * @param geomIdx The geometry index
 * @return true if constraint was added, false otherwise
 */
static bool
extractLineStringConstraints(const LineString &lineString, CT &ct,
                             std::vector<ConstraintInfoType> &constraintInfos,
                             std::size_t                      geomIdx)
{
  if (lineString.numPoints() < 2)
    return false;

  std::vector<CGALPoint> points = extractPointsFromLineString(lineString);

  // Insert constraint with explicit LINESTRING type
  Constraint_id cid = ct.insert_constraint(points.begin(), points.end());
  constraintInfos.emplace_back(cid, geomIdx, 0, 0, 0,
                               ConstraintInfo::Type::LINESTRING);
  return true;
}

/**
 * @brief Extracts polygon exterior ring constraints and adds them to the
 * triangulation
 * @param ring The exterior ring
 * @param ct The constrained triangulation to add to
 * @param constraintInfos Collection of constraint infos
 * @param geomIdx The geometry index
 * @param type The constraint type
 * @return true if constraint was added, false otherwise
 */
static bool
extractPolygonExteriorConstraint(
    const LineString &ring, CT &ct,
    std::vector<ConstraintInfoType> &constraintInfos, std::size_t geomIdx,
    std::size_t polyIdx, ConstraintInfo::Type type)
{
  if (ring.numPoints() < 4)
    return false;

  std::vector<CGALPoint> points = extractPointsFromRing(ring);

  // Insert constraint
  Constraint_id cid = ct.insert_constraint(points.begin(), points.end());
  constraintInfos.emplace_back(cid, geomIdx, polyIdx, 0, 0, type);
  return true;
}

/**
 * @brief Extracts polygon interior ring constraints and adds them to the
 * triangulation
 * @param polygon The source polygon
 * @param ct The constrained triangulation to add to
 * @param constraintInfos Collection of constraint infos
 * @param geomIdx The geometry index
 * @param polyIdx The polygon index in a collection
 * @param type The constraint type for interior rings
 */
static void
extractPolygonInteriorConstraints(
    const Polygon &polygon, CT &ct,
    std::vector<ConstraintInfoType> &constraintInfos, std::size_t geomIdx,
    std::size_t polyIdx, ConstraintInfo::Type type)
{
  // Add interior rings
  for (std::size_t ringIdx = 0; ringIdx < polygon.numInteriorRings();
       ++ringIdx) {
    const auto &ring = polygon.interiorRingN(ringIdx);
    if (ring.numPoints() < 4)
      continue;

    std::vector<CGALPoint> points = extractPointsFromRing(ring);

    // Insert constraint
    Constraint_id cid = ct.insert_constraint(points.begin(), points.end());
    constraintInfos.emplace_back(cid, geomIdx, polyIdx, ringIdx + 1, 0, type);
  }
}

/**
 * @brief Extracts polygon constraints and adds them to the triangulation
 * @param polygon The source polygon
 * @param ct The constrained triangulation to add to
 * @param constraintInfos Collection of constraint infos
 * @param geomIdx The geometry index
 */
static void
extractPolygonConstraints(const Polygon &polygon, CT &ct,
                          std::vector<ConstraintInfoType> &constraintInfos,
                          std::size_t                      geomIdx)
{
  // Add exterior ring
  if (extractPolygonExteriorConstraint(
          polygon.exteriorRing(), ct, constraintInfos, geomIdx, 0,
          ConstraintInfo::Type::POLYGON_EXTERIOR)) {

    // Add interior rings if exterior was added
    extractPolygonInteriorConstraints(polygon, ct, constraintInfos, geomIdx, 0,
                                      ConstraintInfo::Type::POLYGON_INTERIOR);
  }
}

/**
 * @brief Extracts MultiPolygon constraints and adds them to the triangulation
 * @param multiPolygon The source MultiPolygon
 * @param ct The constrained triangulation to add to
 * @param constraintInfos Collection of constraint infos
 * @param geomIdx The geometry index
 */
static void
extractMultiPolygonConstraints(const MultiPolygon &multiPolygon, CT &ct,
                               std::vector<ConstraintInfoType> &constraintInfos,
                               std::size_t                      geomIdx)
{
  // Process each polygon in the MultiPolygon
  for (std::size_t polyIdx = 0; polyIdx < multiPolygon.numGeometries();
       ++polyIdx) {
    const auto &polygon = multiPolygon.polygonN(polyIdx);

    // Add exterior ring
    if (extractPolygonExteriorConstraint(
            polygon.exteriorRing(), ct, constraintInfos, geomIdx, polyIdx,
            ConstraintInfo::Type::MULTIPOLYGON_EXTERIOR)) {

      // Add interior rings if exterior was added
      extractPolygonInteriorConstraints(
          polygon, ct, constraintInfos, geomIdx, polyIdx,
          ConstraintInfo::Type::MULTIPOLYGON_INTERIOR);
    }
  }
}

/**
 * @brief Extracts PolyhedralSurface constraints and adds them to the
 * triangulation
 * @param surface The source PolyhedralSurface
 * @param ct The constrained triangulation to add to
 * @param constraintInfos Collection of constraint infos
 * @param geomIdx The geometry index
 */
static void
extractPolyhedralSurfaceConstraints(
    const PolyhedralSurface &surface, CT &ct,
    std::vector<ConstraintInfoType> &constraintInfos, std::size_t geomIdx)
{
  // Process each polygon in the PolyhedralSurface
  for (std::size_t polyIdx = 0; polyIdx < surface.numPatches(); ++polyIdx) {
    const auto &polygon = surface.patchN(polyIdx);

    // Add exterior ring
    if (extractPolygonExteriorConstraint(
            polygon.exteriorRing(), ct, constraintInfos, geomIdx, polyIdx,
            ConstraintInfo::Type::POLYHEDRALSURFACE_EXTERIOR)) {

      // Add interior rings if exterior was added
      extractPolygonInteriorConstraints(
          polygon, ct, constraintInfos, geomIdx, polyIdx,
          ConstraintInfo::Type::POLYHEDRALSURFACE_INTERIOR);
    }
  }
}

/**
 * @brief Extracts constraints from all geometry types in a collection
 * @param collection The source GeometryCollection
 * @param ct The constrained triangulation to add to
 * @param constraintInfos Collection of constraint infos
 */
static void
extractAllConstraints(const GeometryCollection &collection, CT &ct,
                      std::vector<ConstraintInfoType> &constraintInfos)
{
  for (std::size_t geomIdx = 0; geomIdx < collection.numGeometries();
       ++geomIdx) {
    const auto &geom = collection.geometryN(geomIdx);

    switch (geom.geometryTypeId()) {
    case TYPE_LINESTRING:
      extractLineStringConstraints(static_cast<const LineString &>(geom), ct,
                                   constraintInfos, geomIdx);
      break;

    case TYPE_POLYGON:
      extractPolygonConstraints(static_cast<const Polygon &>(geom), ct,
                                constraintInfos, geomIdx);
      break;

    case TYPE_MULTIPOLYGON:
      extractMultiPolygonConstraints(static_cast<const MultiPolygon &>(geom),
                                     ct, constraintInfos, geomIdx);
      break;

    case TYPE_POLYHEDRALSURFACE:
      extractPolyhedralSurfaceConstraints(
          static_cast<const PolyhedralSurface &>(geom), ct, constraintInfos,
          geomIdx);
      break;

    default:
      // Other types are not processed for simplification with topology
      // preservation
      break;
    }
  }
}

/**
 * @brief Organizes simplified constraints into corresponding data structures
 * @param ct The constrained triangulation with simplified constraints
 * @param constraintInfos Collection of constraint infos
 * @param store The segment store for interpolation
 * @param dimension The coordinate dimension type
 * @param linestrings Output map for linestrings
 * @param polygonRings Output map for polygon rings
 * @param multiPolygonRings Output map for multipolygon rings
 * @param polyhedralSurfaceRings Output map for polyhedral surface rings
 */
static void
organizeSimplifiedConstraints(
    const CT &ct, const std::vector<ConstraintInfoType> &constraintInfos,
    const SegmentStore &store, CoordinateType dimension,
    std::map<size_t, std::unique_ptr<LineString>>  &linestrings,
    std::map<size_t, std::map<size_t, LineString>> &polygonRings,
    std::map<size_t, std::map<size_t, std::map<size_t, LineString>>>
        &multiPolygonRings,
    std::map<size_t, std::map<size_t, std::map<size_t, LineString>>>
        &polyhedralSurfaceRings)
{

  for (const auto &info : constraintInfos) {
    // Build the simplified linestring from the constraint
    LineString simplifiedRing =
        buildLineStringFromConstraint(ct, info.cid, store, dimension);

    // Organize based on the constraint type
    switch (info.type) {
    case ConstraintInfo::Type::LINESTRING:
      linestrings[info.geomIndex] =
          std::make_unique<LineString>(simplifiedRing);
      break;

    case ConstraintInfo::Type::POLYGON_EXTERIOR:
    case ConstraintInfo::Type::POLYGON_INTERIOR:
      polygonRings[info.geomIndex][info.ringIndex] = simplifiedRing;
      break;

    case ConstraintInfo::Type::POLYHEDRALSURFACE_EXTERIOR:
    case ConstraintInfo::Type::POLYHEDRALSURFACE_INTERIOR:
      polyhedralSurfaceRings[info.geomIndex][info.polyIndex][info.ringIndex] =
          simplifiedRing;
      break;

    case ConstraintInfo::Type::MULTIPOLYGON_EXTERIOR:
    case ConstraintInfo::Type::MULTIPOLYGON_INTERIOR:
      multiPolygonRings[info.geomIndex][info.polyIndex][info.ringIndex] =
          simplifiedRing;
      break;

    default:
      std::cerr << "Unknown constraint type encountered during simplification."
                << std::endl;
      break;
    }
  }
}

/**
 * @brief Reconstructs a LineString in a GeometryCollection from simplified
 * constraints
 * @param linestrings Map of simplified linestrings
 * @param collection The original GeometryCollection
 * @param geomIdx The geometry index
 * @param result The output GeometryCollection
 */
static void
reconstructLineString(
    const std::map<size_t, std::unique_ptr<LineString>> &linestrings,
    const GeometryCollection &collection, std::size_t geomIdx,
    GeometryCollection &result)
{

  auto it = linestrings.find(geomIdx);
  if (it != linestrings.end()) {
    // Use the simplified linestring
    result.addGeometry(*(it->second));
  } else {
    // No matching simplified linestring, use the original
    result.addGeometry(collection.geometryN(geomIdx));
  }
}

/**
 * @brief Reconstructs a Polygon in a GeometryCollection from simplified rings
 * @param polygonRings Map of simplified polygon rings
 * @param collection The original GeometryCollection
 * @param geomIdx The geometry index
 * @param result The output GeometryCollection
 */
static void
reconstructPolygon(
    const std::map<size_t, std::map<size_t, LineString>> &polygonRings,
    const GeometryCollection &collection, std::size_t geomIdx,
    GeometryCollection &result)
{

  auto it = polygonRings.find(geomIdx);
  if (it != polygonRings.end() && !it->second.empty()) {
    // We have at least one simplified ring
    auto polygonPtr = std::make_unique<Polygon>();

    // Set the exterior ring if it exists
    auto exteriorIt = it->second.find(0);
    if (exteriorIt != it->second.end()) {
      polygonPtr->setExteriorRing(exteriorIt->second);
    }

    // Add interior rings
    for (const auto &[ringIdx, ring] : it->second) {
      if (ringIdx > 0) {
        polygonPtr->addInteriorRing(ring);
      }
    }

    result.addGeometry(*polygonPtr);
  } else {
    // No matching simplified polygon, use the original
    result.addGeometry(collection.geometryN(geomIdx));
  }
}

/**
 * @brief Reconstructs a MultiPolygon in a GeometryCollection from simplified
 * rings, preserving dimensions
 * @param multiPolygonRings Map of simplified multipolygon rings
 * @param collection The original GeometryCollection
 * @param geomIdx The geometry index
 * @param result The output GeometryCollection
 */
static void
reconstructMultiPolygon(
    const std::map<size_t, std::map<size_t, std::map<size_t, LineString>>>
                             &multiPolygonRings,
    const GeometryCollection &collection, std::size_t geomIdx,
    GeometryCollection &result)
{

  auto        it = multiPolygonRings.find(geomIdx);
  const auto &originalMultiPolygon =
      static_cast<const MultiPolygon &>(collection.geometryN(geomIdx));

  if (it != multiPolygonRings.end() && !it->second.empty()) {
    // We have at least one simplified polygon in the MultiPolygon
    auto multiPolygonPtr = std::make_unique<MultiPolygon>();

    // Process each polygon in the original MultiPolygon
    for (std::size_t polyIdx = 0;
         polyIdx < originalMultiPolygon.numGeometries(); ++polyIdx) {
      auto polyIt = it->second.find(polyIdx);

      if (polyIt != it->second.end() && !polyIt->second.empty()) {
        // We have at least one simplified ring for this polygon
        Polygon simplifiedPolygon;

        // Set exterior ring if it exists
        auto exteriorIt = polyIt->second.find(0);
        if (exteriorIt != polyIt->second.end()) {
          // Make sure the ring has the right dimension type
          simplifiedPolygon.setExteriorRing(exteriorIt->second);
        }

        // Add interior rings
        for (const auto &[ringIdx, ring] : polyIt->second) {
          if (ringIdx > 0) {
            LineString ringCopy(ring);
            if (polyIdx < originalMultiPolygon.numGeometries() &&
                ringIdx - 1 <
                    originalMultiPolygon.polygonN(polyIdx).numInteriorRings()) {
            }
            simplifiedPolygon.addInteriorRing(ringCopy);
          }
        }

        multiPolygonPtr->addGeometry(simplifiedPolygon);
      } else {
        // No simplified version of this polygon, use the original
        multiPolygonPtr->addGeometry(originalMultiPolygon.polygonN(polyIdx));
      }
    }

    result.addGeometry(*multiPolygonPtr);
  } else {
    // No matching simplified MultiPolygon, use the original
    result.addGeometry(originalMultiPolygon);
  }
}

/**
 * @brief Reconstructs a PolyhedralSurface in a GeometryCollection from
 * simplified rings, preserving dimensions
 * @param polyhedralSurfaceRings Map of simplified polyhedral surface rings
 * @param collection The original GeometryCollection
 * @param geomIdx The geometry index
 * @param result The output GeometryCollection
 */
static void
reconstructPolyhedralSurface(
    const std::map<size_t, std::map<size_t, std::map<size_t, LineString>>>
                             &polyhedralSurfaceRings,
    const GeometryCollection &collection, std::size_t geomIdx,
    GeometryCollection &result)
{

  auto        it = polyhedralSurfaceRings.find(geomIdx);
  const auto &originalSurface =
      static_cast<const PolyhedralSurface &>(collection.geometryN(geomIdx));

  if (it != polyhedralSurfaceRings.end() && !it->second.empty()) {
    // We have at least one simplified polygon in the PolyhedralSurface
    auto surfacePtr = std::make_unique<PolyhedralSurface>();

    // Process each polygon in the original PolyhedralSurface
    for (std::size_t polyIdx = 0; polyIdx < originalSurface.numPatches();
         ++polyIdx) {
      auto polyIt = it->second.find(polyIdx);

      if (polyIt != it->second.end() && !polyIt->second.empty()) {
        // We have at least one simplified ring for this polygon
        Polygon simplifiedPolygon;

        // Set exterior ring if it exists
        auto exteriorIt = polyIt->second.find(0);
        if (exteriorIt != polyIt->second.end()) {
          simplifiedPolygon.setExteriorRing(exteriorIt->second);
        }

        // Add interior rings
        for (const auto &[ringIdx, ring] : polyIt->second) {
          if (ringIdx > 0) {
            LineString ringCopy(ring);
            if (polyIdx < originalSurface.numPatches() &&
                ringIdx - 1 <
                    originalSurface.patchN(polyIdx).numInteriorRings()) {
            }
            simplifiedPolygon.addInteriorRing(ringCopy);
          }
        }

        surfacePtr->addPatch(simplifiedPolygon);
      } else {
        // No simplified version of this polygon, use the original
        surfacePtr->addPatch(originalSurface.patchN(polyIdx));
      }
    }

    result.addGeometry(*surfacePtr);
  } else {
    // No matching simplified PolyhedralSurface, use the original
    result.addGeometry(originalSurface);
  }
}

/**
 * @brief Reconstructs all geometries from simplified constraints, preserving
 * dimensions
 * @param collection The original GeometryCollection
 * @param linestrings Map of simplified linestrings
 * @param polygonRings Map of simplified polygon rings
 * @param multiPolygonRings Map of simplified multipolygon rings
 * @param polyhedralSurfaceRings Map of simplified polyhedral surface rings
 * @return A unique_ptr to a GeometryCollection with reconstructed geometries
 */
static std::unique_ptr<GeometryCollection>
reconstructAllGeometries(
    const GeometryCollection                             &collection,
    const std::map<size_t, std::unique_ptr<LineString>>  &linestrings,
    const std::map<size_t, std::map<size_t, LineString>> &polygonRings,
    const std::map<size_t, std::map<size_t, std::map<size_t, LineString>>>
        &multiPolygonRings,
    const std::map<size_t, std::map<size_t, std::map<size_t, LineString>>>
        &polyhedralSurfaceRings)
{

  auto result = std::make_unique<GeometryCollection>();

  for (std::size_t geomIdx = 0; geomIdx < collection.numGeometries();
       ++geomIdx) {
    const auto &geom = collection.geometryN(geomIdx);

    switch (geom.geometryTypeId()) {
    case TYPE_LINESTRING:
      reconstructLineString(linestrings, collection, geomIdx, *result);
      break;

    case TYPE_POLYGON:
      reconstructPolygon(polygonRings, collection, geomIdx, *result);
      break;

    case TYPE_MULTIPOLYGON:
      reconstructMultiPolygon(multiPolygonRings, collection, geomIdx, *result);
      break;

    case TYPE_POLYHEDRALSURFACE:
      reconstructPolyhedralSurface(polyhedralSurfaceRings, collection, geomIdx,
                                   *result);
      break;

    default:
      // For other geometry types, just clone the original
      result->addGeometry(geom);
      break;
    }
  }

  return result;
}

/**
 * @brief Simplifies a LineString
 * @param lineString The LineString to simplify
 * @param threshold The simplification threshold
 * @param preserveTopology Whether to preserve topology
 * @param store The segment store for interpolation
 * @return A simplified geometry
 */
auto
simplifyLineString(const LineString &lineString, double           threshold,
                   bool /*preserveTopology*/, const SegmentStore &store)
    -> std::unique_ptr<Geometry>
{
  CoordinateType dimension = lineString.getCoordinateType();

  // Create triangulation for the LineString
  CT                     ct;
  std::vector<CGALPoint> points = extractPointsFromLineString(lineString);

  // Insert constraint and simplify
  Constraint_id cid = ct.insert_constraint(points.begin(), points.end());
  PS::simplify(ct, Cost(), Stop(threshold));

  // Reconstruct simplified LineString
  auto result = std::make_unique<LineString>();

  for (auto vit = ct.points_in_constraint_begin(cid);
       vit != ct.points_in_constraint_end(cid); ++vit) {

    double x = CGAL::to_double(vit->x());
    double y = CGAL::to_double(vit->y());

    // Create point with interpolated ZM values
    result->addPoint(store.createPoint(x, y, dimension));
  }

  return result;
}

/**
 * @brief Simplifies a MultiLineString
 * @param multiLine The MultiLineString to simplify
 * @param threshold The simplification threshold
 * @param preserveTopology Whether to preserve topology
 * @return A simplified geometry
 */
auto
simplifyMultiLineString(const MultiLineString &multiLine, double threshold,
                        bool preserveTopology) -> std::unique_ptr<Geometry>
{
  // Extract dimension info
  CoordinateType dimension = multiLine.getCoordinateType();

  // Extract segments for interpolation
  SegmentStore store;
  store.extractSegments(multiLine);

  if (preserveTopology) {
    // Topology-preserving mode: simplify all LineStrings together
    CT                              ct;
    std::vector<ConstraintInfoType> constraintInfos;

    // Add all LineStrings as constraints
    for (size_t i = 0; i < multiLine.numGeometries(); ++i) {
      extractLineStringConstraints(multiLine.lineStringN(i), ct,
                                   constraintInfos, i);
    }

    // Simplify all constraints
    PS::simplify(ct, Cost(), Stop(threshold));

    // Sort constraints by original order
    std::sort(constraintInfos.begin(), constraintInfos.end(),
              ConstraintInfoCompare<Constraint_id>());

    // Reconstruct simplified MultiLineString
    auto result = std::make_unique<MultiLineString>();

    for (const auto &info : constraintInfos) {
      auto ls = std::make_unique<LineString>(
          buildLineStringFromConstraint(ct, info.cid, store, dimension));
      result->addGeometry(*ls);
    }

    return result;
  } else {
    // Non-topology preserving mode: simplify each LineString independently
    auto result = std::make_unique<MultiLineString>();

    for (size_t i = 0; i < multiLine.numGeometries(); ++i) {
      const auto &ls    = multiLine.lineStringN(i);
      auto simplifiedLs = simplifyLineString(ls, threshold, false, store);
      if (auto *linestring =
              dynamic_cast<const LineString *>(simplifiedLs.get())) {
        result->addGeometry(*linestring);
      }
    }

    return result;
  }
}

/**
 * @brief Simplifies a Polygon
 * @param polygon The Polygon to simplify
 * @param threshold The simplification threshold
 * @param preserveTopology Whether to preserve topology
 * @return A simplified geometry
 */
auto
simplifyPolygon(const Polygon &polygon, double threshold,
                bool /*preserveTopology*/) -> std::unique_ptr<Geometry>
{
  // Extract dimension info
  CoordinateType dimension = polygon.getCoordinateType();

  // Extract segments for interpolation
  SegmentStore store;
  store.extractSegments(polygon);

  // Create triangulation
  CT                              ct;
  std::vector<ConstraintInfoType> constraintInfos;

  // Add polygon constraints
  extractPolygonConstraints(polygon, ct, constraintInfos, 0);

  // Simplify with squared threshold
  PS::simplify(ct, Cost(), Stop(threshold));

  // Reconstruct the polygon with correct ring order
  auto result = std::make_unique<Polygon>();

  bool exteriorRingSet = false;
  // Process constraints in order
  for (const auto &info : constraintInfos) {
    LineString simplifiedRing =
        buildLineStringFromConstraint(ct, info.cid, store, dimension);

    if (info.type == ConstraintInfo::Type::POLYGON_EXTERIOR &&
        !exteriorRingSet) {
      // Exterior ring
      result->setExteriorRing(simplifiedRing);
      exteriorRingSet = true;
    } else {
      // Interior ring
      result->addInteriorRing(simplifiedRing);
    }
  }

  return result;
}

/**
 * @brief Simplifies a MultiPolygon
 * @param multiPolygon The MultiPolygon to simplify
 * @param threshold The simplification threshold
 * @param preserveTopology Whether to preserve topology
 * @return A simplified geometry
 */
auto
simplifyMultiPolygon(const MultiPolygon &multiPolygon, double threshold,
                     bool preserveTopology) -> std::unique_ptr<Geometry>
{
  // Extract dimension info
  CoordinateType dimension = multiPolygon.getCoordinateType();

  // Extract segments for interpolation
  SegmentStore store;
  store.extractSegments(multiPolygon);

  if (preserveTopology) {
    // Topology-preserving mode: simplify all polygons together
    CT                              ct;
    std::vector<ConstraintInfoType> constraintInfos;

    // Extract all polygon constraints
    extractMultiPolygonConstraints(multiPolygon, ct, constraintInfos, 0);

    // Simplify
    PS::simplify(ct, Cost(), Stop(threshold));

    // Sort constraints by original order
    std::sort(constraintInfos.begin(), constraintInfos.end(),
              ConstraintInfoCompare<Constraint_id>());

    // Create a GeometryCollection to use our reconstructAllGeometries function
    GeometryCollection tempCollection;
    tempCollection.addGeometry(multiPolygon);

    // Create empty data structures for other geometry types
    std::map<size_t, std::unique_ptr<LineString>>  emptyLinestrings;
    std::map<size_t, std::map<size_t, LineString>> emptyPolygonRings;
    std::map<size_t, std::map<size_t, std::map<size_t, LineString>>>
        multiPolygonRings;
    std::map<size_t, std::map<size_t, std::map<size_t, LineString>>>
        emptyPolyhedralRings;

    // Organize the constraints for this single MultiPolygon
    // Process constraints in order
    for (const auto &info : constraintInfos) {
      // Build simplified ring
      LineString simplifiedRing =
          buildLineStringFromConstraint(ct, info.cid, store, dimension);

      // Add to appropriate data structure
      multiPolygonRings[0][info.polyIndex][info.ringIndex] = simplifiedRing;
    }

    // Reconstruct collection with just this MultiPolygon
    auto result = reconstructAllGeometries(tempCollection, emptyLinestrings,
                                           emptyPolygonRings, multiPolygonRings,
                                           emptyPolyhedralRings);

    // Extract the MultiPolygon from the collection
    if (result->numGeometries() > 0) {
      auto *mp = dynamic_cast<const MultiPolygon *>(&(result->geometryN(0)));
      if (mp) {
        return std::unique_ptr<Geometry>(mp->clone());
      }
    }

    // Fallback if reconstruction failed
    return std::unique_ptr<Geometry>(multiPolygon.clone());
  } else {
    // Non-topology preserving mode: simplify each polygon independently
    auto result = std::make_unique<MultiPolygon>();

    for (size_t i = 0; i < multiPolygon.numGeometries(); ++i) {
      const auto &polygon = multiPolygon.polygonN(i);
      auto simplified = simplifyPolygon(polygon, threshold, preserveTopology);
      if (auto *poly = dynamic_cast<const Polygon *>(simplified.get())) {
        result->addGeometry(*poly);
      }
    }

    return result;
  }
}

/**
 * @brief Simplifies a PolyhedralSurface
 * @param polySurface The PolyhedralSurface to simplify
 * @param threshold The simplification threshold
 * @param preserveTopology Whether to preserve topology
 * @return A simplified geometry
 */
auto
simplifyPolyhedralSurface(const PolyhedralSurface &polySurface,
                          double threshold, bool preserveTopology)
    -> std::unique_ptr<Geometry>
{
  // Extract dimension info
  CoordinateType dimension = polySurface.getCoordinateType();

  // Extract segments for interpolation
  SegmentStore store;
  store.extractSegments(polySurface);

  if (preserveTopology) {
    // Topology-preserving mode: use direct approach
    CT                              ct;
    std::vector<ConstraintInfoType> constraintInfos;

    // Extract all constraints
    extractPolyhedralSurfaceConstraints(polySurface, ct, constraintInfos, 0);

    // Simplify
    PS::simplify(ct, Cost(), Stop(threshold));

    // Sort constraints by original order
    std::sort(constraintInfos.begin(), constraintInfos.end(),
              ConstraintInfoCompare<Constraint_id>());

    // Create a GeometryCollection to use our reconstructAllGeometries function
    GeometryCollection tempCollection;
    tempCollection.addGeometry(polySurface);

    // Create empty data structures for other geometry types
    std::map<size_t, std::unique_ptr<LineString>>  emptyLinestrings;
    std::map<size_t, std::map<size_t, LineString>> emptyPolygonRings;
    std::map<size_t, std::map<size_t, std::map<size_t, LineString>>>
        emptyMultiPolygonRings;
    std::map<size_t, std::map<size_t, std::map<size_t, LineString>>>
        polyhedralSurfaceRings;

    // Organize the constraints for this single PolyhedralSurface
    // Process constraints in order
    for (const auto &info : constraintInfos) {
      // Build simplified ring
      LineString simplifiedRing =
          buildLineStringFromConstraint(ct, info.cid, store, dimension);

      // Add to appropriate data structure
      polyhedralSurfaceRings[0][info.polyIndex][info.ringIndex] =
          simplifiedRing;
    }

    // Reconstruct collection with just this PolyhedralSurface
    auto result = reconstructAllGeometries(
        tempCollection, emptyLinestrings, emptyPolygonRings,
        emptyMultiPolygonRings, polyhedralSurfaceRings);

    // Extract the PolyhedralSurface from the collection
    if (result->numGeometries() > 0) {
      auto *ps =
          dynamic_cast<const PolyhedralSurface *>(&(result->geometryN(0)));
      if (ps) {
        return std::unique_ptr<Geometry>(ps->clone());
      }
    }

    // Fallback if reconstruction failed
    return std::unique_ptr<Geometry>(polySurface.clone());
  } else {
    // Create a MultiPolygon containing the same polygons and preserving
    // coordinate type
    MultiPolygon multiPolygon;

    for (size_t i = 0; i < polySurface.numPatches(); ++i) {
      Polygon polygon = polySurface.patchN(i);
      multiPolygon.addGeometry(polygon);
    }

    // Simplify the MultiPolygon
    auto simplifiedMP =
        simplifyMultiPolygon(multiPolygon, threshold, preserveTopology);
    auto *simplifiedMultiPolygon =
        dynamic_cast<const MultiPolygon *>(simplifiedMP.get());

    if (!simplifiedMultiPolygon) {
      // Fallback if casting fails
      return std::unique_ptr<Geometry>(polySurface.clone());
    }

    // Create a new PolyhedralSurface from the simplified MultiPolygon
    auto result = std::make_unique<PolyhedralSurface>();

    for (size_t i = 0; i < simplifiedMultiPolygon->numGeometries(); ++i) {
      const auto &polygon = simplifiedMultiPolygon->polygonN(i);
      Polygon     polygonCopy(polygon);
      result->addPatch(polygonCopy);
    }

    return result;
  }
}

/**
 * @brief Simplifies a GeometryCollection with topology preservation
 * @param collection The GeometryCollection to simplify
 * @param threshold The simplification threshold
 * @return A simplified geometry
 */
auto
simplifyGeometryCollectionTopology(const GeometryCollection &collection,
                                   double                    threshold)
    -> std::unique_ptr<Geometry>
{
  // Extract dimension info
  CoordinateType dimension = collection.getCoordinateType();

  // Extract segments for interpolation
  SegmentStore store;
  store.extractSegments(collection);

  // Create triangulation to hold all constraints
  CT                              ct;
  std::vector<ConstraintInfoType> constraintInfos;

  // Step 1: Extract all constraints from geometries
  extractAllConstraints(collection, ct, constraintInfos);

  // Step 2: Simplify all constraints together
  PS::simplify(ct, Cost(), Stop(threshold));

  // Step 3: Sort constraints by original order
  std::sort(constraintInfos.begin(), constraintInfos.end(),
            ConstraintInfoCompare<Constraint_id>());

  // Step 4: Create structures to hold matched geometry components
  std::map<size_t, std::unique_ptr<LineString>>  linestrings;
  std::map<size_t, std::map<size_t, LineString>> polygonRings;
  std::map<size_t, std::map<size_t, std::map<size_t, LineString>>>
      multiPolygonRings;
  std::map<size_t, std::map<size_t, std::map<size_t, LineString>>>
      polyhedralSurfaceRings;

  // Step 5: Organize simplified constraints by geometry type
  organizeSimplifiedConstraints(ct, constraintInfos, store, dimension,
                                linestrings, polygonRings, multiPolygonRings,
                                polyhedralSurfaceRings);

  // Step 6: Reconstruct geometries in original order
  return reconstructAllGeometries(collection, linestrings, polygonRings,
                                  multiPolygonRings, polyhedralSurfaceRings);
}

/**
 * @brief Wrapper for the algorithm::simplify function
 * @param geometry The geometry to simplify
 * @param threshold The simplification threshold
 * @param preserveTopology Whether to preserve topology
 * @return A simplified geometry
 */
inline auto
simplify(const Geometry &geometry, double threshold, bool preserveTopology)
    -> std::unique_ptr<Geometry>
{
  return algorithm::simplify(geometry, threshold, preserveTopology);
}

/**
 * @brief Simplifies a GeometryCollection
 * @param collection The GeometryCollection to simplify
 * @param threshold The simplification threshold
 * @param preserveTopology Whether to preserve topology
 * @return A simplified geometry
 */
auto
simplifyGeometryCollection(const GeometryCollection &collection,
                           double threshold, bool preserveTopology)
    -> std::unique_ptr<Geometry>
{
  if (preserveTopology) {
    return simplifyGeometryCollectionTopology(collection, threshold);
  } else {
    // Non-topology preserving: simplify each geometry independently
    auto result = std::make_unique<GeometryCollection>();

    for (std::size_t i = 0; i < collection.numGeometries(); ++i) {
      const auto &geom       = collection.geometryN(i);
      auto        simplified = simplify(geom, threshold, preserveTopology);
      result->addGeometry(*simplified);
    }

    return result;
  }
}

} // namespace detail
} // namespace SFCGAL
