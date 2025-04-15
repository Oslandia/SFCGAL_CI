#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/algorithm/simplification.h"
#include "SFCGAL/detail/SegmentStore.h"
#include "SFCGAL/detail/ConstraintInfo.h"

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
 * @brief Simplifies a LineString
 */
auto
simplifyLineString(const LineString &lineString, double           threshold,
                   bool /*preserveTopology*/, const SegmentStore &store)
    -> std::unique_ptr<Geometry>
{
  CoordinateType dimension = lineString.getCoordinateType();

  // Create triangulation for the LineString
  CT                     ct;
  std::vector<CGALPoint> points;
  points.reserve(lineString.numPoints());

  // Extract all points
  for (size_t i = 0; i < lineString.numPoints(); ++i) {
    const auto &pt = lineString.pointN(i);
    points.emplace_back(CGAL::to_double(pt.x()), CGAL::to_double(pt.y()),
                        pt.is3D() ? CGAL::to_double(pt.z()) : 0.0);
  }

  // Insert constraint and simplify
  Constraint_id cid = ct.insert_constraint(points.begin(), points.end());
  PS::simplify(ct, Cost(), Stop(threshold));

  // Reconstruct simplified LineString
  auto result = std::make_unique<LineString>();

  for (auto vit = ct.points_in_constraint_begin(cid);
       vit != ct.points_in_constraint_end(cid); ++vit) {

    double x = CGAL::to_double(vit->x());
    double y = CGAL::to_double(vit->y());
    double z = CGAL::to_double(vit->z());

    // Create point with interpolated ZM values
    result->addPoint(createPoint(x, y, z, store, dimension));
  }

  return result;
}

/**
 * @brief Simplifies a MultiLineString
 */
auto
simplifyMultiLineString(const MultiLineString &multiLine, double threshold,
                        bool preserveTopology) -> std::unique_ptr<Geometry>
{
  // Extract dimension info
  CoordinateType dimension = multiLine.getCoordinateType();

  // Extract segments for interpolation
  SegmentStore store;
  extractSegments(multiLine, store);

  if (preserveTopology) {
    // Topology-preserving mode: simplify all LineStrings together
    CT                               ct;
    std::vector<ConstraintInfoType> constraintInfos;

    // Add all LineStrings as constraints
    for (size_t i = 0; i < multiLine.numGeometries(); ++i) {
      const auto            &ls = multiLine.lineStringN(i);
      std::vector<CGALPoint> points;
      points.reserve(ls.numPoints());

      for (size_t j = 0; j < ls.numPoints(); ++j) {
        const auto &pt = ls.pointN(j);
        points.emplace_back(CGAL::to_double(pt.x()), CGAL::to_double(pt.y()),
                            pt.is3D() ? CGAL::to_double(pt.z()) : 0.0);
      }

      // Insert constraint and record order
      Constraint_id cid = ct.insert_constraint(points.begin(), points.end());
      constraintInfos.emplace_back(cid, i, 0, 0, 0,
                                   ConstraintInfo::Type::LINESTRING);
    }

    // Simplify all constraints
    PS::simplify(ct, Cost(), Stop(threshold));

    // Sort constraints by original order
    std::sort(constraintInfos.begin(), constraintInfos.end(),
          [](const ConstraintInfoType& a, const ConstraintInfoType& b) {
            return compareConstraintInfo<Constraint_id>(a, b);
          });

    // Reconstruct simplified MultiLineString
    auto result = std::make_unique<MultiLineString>();

    for (const auto &info : constraintInfos) {
      auto ls = std::make_unique<LineString>();

      for (auto vit = ct.points_in_constraint_begin(info.cid);
           vit != ct.points_in_constraint_end(info.cid); ++vit) {

        double x = CGAL::to_double(vit->x());
        double y = CGAL::to_double(vit->y());
        double z = CGAL::to_double(vit->z());

        // Create point with interpolated ZM values
        ls->addPoint(createPoint(x, y, z, store, dimension));
      }

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
 */
auto
simplifyPolygon(const Polygon &polygon, double threshold,
                bool /*preserveTopology*/) -> std::unique_ptr<Geometry>
{
  // Extract dimension info
  CoordinateType dimension = polygon.getCoordinateType();

  // Extract segments for interpolation
  SegmentStore store;
  extractSegments(polygon, store);

  // Create triangulation
  CT                               ct;
  std::vector<ConstraintInfoType> constraintInfos;

  // Add exterior ring
  {
    const auto            &ring = polygon.exteriorRing();
    std::vector<CGALPoint> points;
    points.reserve(ring.numPoints());

    for (size_t i = 0; i < ring.numPoints(); ++i) {
      const auto &pt = ring.pointN(i);
      points.emplace_back(CGAL::to_double(pt.x()), CGAL::to_double(pt.y()),
                          pt.is3D() ? CGAL::to_double(pt.z()) : 0.0);
    }

    Constraint_id cid = ct.insert_constraint(points.begin(), points.end());
    constraintInfos.emplace_back(cid, 0, 0, 0, 0,
                                 ConstraintInfo::Type::POLYGON_EXTERIOR);
  }

  // Add interior rings
  for (size_t ringIdx = 0; ringIdx < polygon.numInteriorRings(); ++ringIdx) {
    const auto            &ring = polygon.interiorRingN(ringIdx);
    std::vector<CGALPoint> points;
    points.reserve(ring.numPoints());

    for (size_t i = 0; i < ring.numPoints(); ++i) {
      const auto &pt = ring.pointN(i);
      points.emplace_back(CGAL::to_double(pt.x()), CGAL::to_double(pt.y()),
                          pt.is3D() ? CGAL::to_double(pt.z()) : 0.0);
    }

    Constraint_id cid = ct.insert_constraint(points.begin(), points.end());
    constraintInfos.emplace_back(cid, 0, 0, ringIdx + 1, 0,
                                 ConstraintInfo::Type::POLYGON_INTERIOR);
  }

  // Simplify with squared threshold
  PS::simplify(ct, Cost(), Stop(threshold));

  // Reconstruct the polygon with correct ring order
  auto result = std::make_unique<Polygon>();

  bool exteriorRingSet = false;

  // Process constraints in order
  for (const auto &info : constraintInfos) {
    LineString simplifiedRing;

    for (auto vit = ct.points_in_constraint_begin(info.cid);
         vit != ct.points_in_constraint_end(info.cid); ++vit) {
      double x = CGAL::to_double(vit->x());
      double y = CGAL::to_double(vit->y());
      double z = CGAL::to_double(vit->z());
      simplifiedRing.addPoint(
          createPoint(x, y, z, store, dimension));
    }

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
 */
auto
simplifyMultiPolygon(const MultiPolygon &multiPolygon, double threshold,
                     bool preserveTopology) -> std::unique_ptr<Geometry>
{
  // Extract dimension info
  CoordinateType dimension = multiPolygon.getCoordinateType();

  // Extract segments for interpolation
  SegmentStore store;
  extractSegments(multiPolygon, store);

  if (preserveTopology) {
    // Topology-preserving mode: simplify all polygons together
    CT                               ct;
    std::vector<ConstraintInfoType> constraintInfos;

    // For each polygon in the MultiPolygon
    for (size_t polyIdx = 0; polyIdx < multiPolygon.numGeometries();
         ++polyIdx) {
      const auto &polygon = multiPolygon.polygonN(polyIdx);

      // Add exterior ring
      {
        const auto            &ring = polygon.exteriorRing();
        std::vector<CGALPoint> points;
        points.reserve(ring.numPoints());

        for (size_t i = 0; i < ring.numPoints(); ++i) {
          const auto &pt = ring.pointN(i);
          points.emplace_back(CGAL::to_double(pt.x()), CGAL::to_double(pt.y()),
                              pt.is3D() ? CGAL::to_double(pt.z()) : 0.0);
        }

        // Insert and record order
        Constraint_id cid = ct.insert_constraint(points.begin(), points.end());
        constraintInfos.emplace_back(cid, polyIdx, 0, 0, 0,
                                     ConstraintInfo::Type::POLYGON_EXTERIOR);
      }

      // Add interior rings
      for (size_t ringIdx = 0; ringIdx < polygon.numInteriorRings();
           ++ringIdx) {
        const auto            &ring = polygon.interiorRingN(ringIdx);
        std::vector<CGALPoint> points;
        points.reserve(ring.numPoints());

        for (size_t i = 0; i < ring.numPoints(); ++i) {
          const auto &pt = ring.pointN(i);
          points.emplace_back(CGAL::to_double(pt.x()), CGAL::to_double(pt.y()),
                              pt.is3D() ? CGAL::to_double(pt.z()) : 0.0);
        }

        // Insert and record order
        Constraint_id cid = ct.insert_constraint(points.begin(), points.end());
        constraintInfos.emplace_back(cid, polyIdx, 0, ringIdx + 1, 0,
                                     ConstraintInfo::Type::POLYGON_INTERIOR);
      }
    }

    // Simplify
    PS::simplify(ct, Cost(), Stop(threshold));

    // Sort constraints by original order
    std::sort(constraintInfos.begin(), constraintInfos.end(),
          [](const ConstraintInfoType& a, const ConstraintInfoType& b) {
            return compareConstraintInfo<Constraint_id>(a, b);
          });

    // Reconstruct simplified MultiPolygon
    auto result = std::make_unique<MultiPolygon>();

    // Temporary map to collect polygons in order
    std::map<size_t, std::unique_ptr<Polygon>> polygons;

    // Process constraints in order
    for (const auto &info : constraintInfos) {
      // Check if we need a new polygon
      if (polygons.find(info.geomIndex) == polygons.end()) {
        polygons[info.geomIndex] = std::make_unique<Polygon>();
      }

      // Build simplified ring
      LineString simplifiedRing;

      for (auto vit = ct.points_in_constraint_begin(info.cid);
           vit != ct.points_in_constraint_end(info.cid); ++vit) {
        double x = CGAL::to_double(vit->x());
        double y = CGAL::to_double(vit->y());
        double z = CGAL::to_double(vit->z());
        simplifiedRing.addPoint(
            createPoint(x, y, z, store, dimension));
      }

      // Add ring to appropriate polygon
      if (info.type == ConstraintInfo::Type::POLYGON_EXTERIOR) {
        // Exterior ring
        polygons[info.geomIndex]->setExteriorRing(simplifiedRing);
      } else {
        // Interior ring
        polygons[info.geomIndex]->addInteriorRing(simplifiedRing);
      }
    }

    // Add polygons to MultiPolygon in order
    for (size_t i = 0; i < multiPolygon.numGeometries(); ++i) {
      auto it = polygons.find(i);
      if (it != polygons.end()) {
        result->addGeometry(*(it->second));
      }
    }

    return result;
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
 * @brief Simplifies a GeometryCollection with topology preservation
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
  extractSegments(collection, store);

  // Create triangulation to hold all constraints
  CT                               ct;
  std::vector<ConstraintInfoType> constraintInfos;

  // Extract all LineStrings and polygon boundaries
  for (std::size_t geomIdx = 0; geomIdx < collection.numGeometries();
       ++geomIdx) {
    const auto &geom = collection.geometryN(geomIdx);

    switch (geom.geometryTypeId()) {
    case TYPE_LINESTRING: {
      const auto &lineString = static_cast<const LineString &>(geom);
      if (lineString.numPoints() < 2)
        continue;

      std::vector<CGALPoint> points;
      points.reserve(lineString.numPoints());

      for (std::size_t i = 0; i < lineString.numPoints(); ++i) {
        const auto &pt = lineString.pointN(i);
        points.emplace_back(CGAL::to_double(pt.x()), CGAL::to_double(pt.y()),
                            pt.is3D() ? CGAL::to_double(pt.z()) : 0.0);
      }

      // Insert constraint with explicit LINESTRING type
      Constraint_id cid = ct.insert_constraint(points.begin(), points.end());
      constraintInfos.emplace_back(cid, geomIdx, 0, 0, 0,
                                   ConstraintInfo::Type::LINESTRING);
      break;
    }

    case TYPE_POLYGON: {
      const auto &polygon = static_cast<const Polygon &>(geom);

      // Add exterior ring
      {
        const auto &ring = polygon.exteriorRing();
        if (ring.numPoints() < 4)
          continue;

        std::vector<CGALPoint> points;
        points.reserve(ring.numPoints());

        for (std::size_t i = 0; i < ring.numPoints(); ++i) {
          const auto &pt = ring.pointN(i);
          points.emplace_back(CGAL::to_double(pt.x()), CGAL::to_double(pt.y()),
                              pt.is3D() ? CGAL::to_double(pt.z()) : 0.0);
        }

        // Insert constraint with explicit POLYGON_EXTERIOR type
        Constraint_id cid = ct.insert_constraint(points.begin(), points.end());
        constraintInfos.emplace_back(cid, geomIdx, 0, 0, 0,
                                     ConstraintInfo::Type::POLYGON_EXTERIOR);
      }

      // Add interior rings
      for (std::size_t ringIdx = 0; ringIdx < polygon.numInteriorRings();
           ++ringIdx) {
        const auto &ring = polygon.interiorRingN(ringIdx);
        if (ring.numPoints() < 4)
          continue;

        std::vector<CGALPoint> points;
        points.reserve(ring.numPoints());

        for (std::size_t i = 0; i < ring.numPoints(); ++i) {
          const auto &pt = ring.pointN(i);
          points.emplace_back(CGAL::to_double(pt.x()), CGAL::to_double(pt.y()),
                              pt.is3D() ? CGAL::to_double(pt.z()) : 0.0);
        }

        // Insert constraint with explicit POLYGON_INTERIOR type
        Constraint_id cid = ct.insert_constraint(points.begin(), points.end());
        constraintInfos.emplace_back(cid, geomIdx, 0, ringIdx + 1, 0,
                                     ConstraintInfo::Type::POLYGON_INTERIOR);
      }
      break;
    }

    case TYPE_MULTIPOLYGON: {
      const auto &multiPolygon = static_cast<const MultiPolygon &>(geom);

      // Process each polygon in the MultiPolygon
      for (std::size_t polyIdx = 0; polyIdx < multiPolygon.numGeometries();
           ++polyIdx) {
        const auto &polygon = multiPolygon.polygonN(polyIdx);

        // Add exterior ring
        {
          const auto &ring = polygon.exteriorRing();
          if (ring.numPoints() < 4)
            continue;

          std::vector<CGALPoint> points;
          points.reserve(ring.numPoints());

          for (std::size_t i = 0; i < ring.numPoints(); ++i) {
            const auto &pt = ring.pointN(i);
            points.emplace_back(CGAL::to_double(pt.x()),
                                CGAL::to_double(pt.y()),
                                pt.is3D() ? CGAL::to_double(pt.z()) : 0.0);
          }

          // Insert constraint with explicit MULTIPOLYGON_EXTERIOR type
          Constraint_id cid =
              ct.insert_constraint(points.begin(), points.end());
          constraintInfos.emplace_back(
              cid, geomIdx, polyIdx, 0, 0,
              ConstraintInfo::Type::MULTIPOLYGON_EXTERIOR);
        }

        // Add interior rings
        for (std::size_t ringIdx = 0; ringIdx < polygon.numInteriorRings();
             ++ringIdx) {
          const auto &ring = polygon.interiorRingN(ringIdx);
          if (ring.numPoints() < 4)
            continue;

          std::vector<CGALPoint> points;
          points.reserve(ring.numPoints());

          for (std::size_t i = 0; i < ring.numPoints(); ++i) {
            const auto &pt = ring.pointN(i);
            points.emplace_back(CGAL::to_double(pt.x()),
                                CGAL::to_double(pt.y()),
                                pt.is3D() ? CGAL::to_double(pt.z()) : 0.0);
          }

          // Insert constraint with explicit MULTIPOLYGON_INTERIOR type
          Constraint_id cid =
              ct.insert_constraint(points.begin(), points.end());
          constraintInfos.emplace_back(
              cid, geomIdx, polyIdx, ringIdx + 1, 0,
              ConstraintInfo::Type::MULTIPOLYGON_INTERIOR);
        }
      }
      break;
    }

    case TYPE_POLYHEDRALSURFACE: {
      const auto &surface = static_cast<const PolyhedralSurface &>(geom);

      // Process each polygon in the PolyhedralSurface
      for (std::size_t polyIdx = 0; polyIdx < surface.numPolygons();
           ++polyIdx) {
        const auto &polygon = surface.polygonN(polyIdx);

        // Add exterior ring
        {
          const auto &ring = polygon.exteriorRing();
          if (ring.numPoints() < 4)
            continue;

          std::vector<CGALPoint> points;
          points.reserve(ring.numPoints());

          for (std::size_t i = 0; i < ring.numPoints(); ++i) {
            const auto &pt = ring.pointN(i);
            points.emplace_back(CGAL::to_double(pt.x()),
                                CGAL::to_double(pt.y()),
                                pt.is3D() ? CGAL::to_double(pt.z()) : 0.0);
          }

          // Insert constraint with POLYHEDRALSURFACE_EXTERIOR type
          Constraint_id cid =
              ct.insert_constraint(points.begin(), points.end());
          constraintInfos.emplace_back(
              cid, geomIdx, polyIdx, 0, 0,
              ConstraintInfo::Type::POLYHEDRALSURFACE_EXTERIOR);
        }

        // Add interior rings
        for (std::size_t ringIdx = 0; ringIdx < polygon.numInteriorRings();
             ++ringIdx) {
          const auto &ring = polygon.interiorRingN(ringIdx);
          if (ring.numPoints() < 4)
            continue;

          std::vector<CGALPoint> points;
          points.reserve(ring.numPoints());

          for (std::size_t i = 0; i < ring.numPoints(); ++i) {
            const auto &pt = ring.pointN(i);
            points.emplace_back(CGAL::to_double(pt.x()),
                                CGAL::to_double(pt.y()),
                                pt.is3D() ? CGAL::to_double(pt.z()) : 0.0);
          }

          // Insert constraint with POLYHEDRALSURFACE_INTERIOR type
          Constraint_id cid =
              ct.insert_constraint(points.begin(), points.end());
          constraintInfos.emplace_back(
              cid, geomIdx, polyIdx, ringIdx + 1, 0,
              ConstraintInfo::Type::POLYHEDRALSURFACE_INTERIOR);
        }
      }
      break;
    }

    default:
      // Other types are not processed for simplification with topology
      // preservation
      break;
    }
  }

  // Simplify all constraints together
  PS::simplify(ct, Cost(), Stop(threshold));

  // Sort constraints by original order
  std::sort(constraintInfos.begin(), constraintInfos.end(),
          [](const ConstraintInfoType& a, const ConstraintInfoType& b) {
            return compareConstraintInfo<Constraint_id>(a, b);
          });

  // Create structures to hold matched geometry components
  std::map<size_t, std::unique_ptr<LineString>> linestrings;
  std::map<size_t, std::map<size_t, LineString>>
      polygonRings; // geomIdx -> (ringIdx -> ring)
  std::map<size_t, std::map<size_t, std::map<size_t, LineString>>>
      multiPolygonRings; // geomIdx -> polyIdx -> ringIdx -> ring
  std::map<size_t, std::map<size_t, std::map<size_t, LineString>>>
      polyhedralSurfaceRings; // geomIdx -> polyIdx -> ringIdx -> ring

  // Reconstruct geometries from simplified constraints
  for (const auto &info : constraintInfos) {
    LineString simplifiedRing;

    for (auto vit = ct.points_in_constraint_begin(info.cid);
         vit != ct.points_in_constraint_end(info.cid); ++vit) {

      double x = CGAL::to_double(vit->x());
      double y = CGAL::to_double(vit->y());
      double z = CGAL::to_double(vit->z());

      simplifiedRing.addPoint(
          createPoint(x, y, z, store, dimension));
    }

    // Check geometry type based on constraint info
    if (info.type == ConstraintInfo::Type::LINESTRING) {
      // LineString
      linestrings[info.geomIndex] =
          std::make_unique<LineString>(simplifiedRing);
    } else if (info.type == ConstraintInfo::Type::POLYGON_EXTERIOR ||
               info.type == ConstraintInfo::Type::POLYGON_INTERIOR) {
      // Polygon
      polygonRings[info.geomIndex][info.ringIndex] = simplifiedRing;
    } else if (info.type == ConstraintInfo::Type::POLYHEDRALSURFACE_EXTERIOR ||
               info.type == ConstraintInfo::Type::POLYHEDRALSURFACE_INTERIOR) {
      // PolyhedralSurface
      polyhedralSurfaceRings[info.geomIndex][info.polyIndex][info.ringIndex] =
          simplifiedRing;
    } else if (info.type == ConstraintInfo::Type::MULTIPOLYGON_EXTERIOR ||
               info.type == ConstraintInfo::Type::MULTIPOLYGON_INTERIOR) {
      // MultiPolygon
      multiPolygonRings[info.geomIndex][info.polyIndex][info.ringIndex] =
          simplifiedRing;
    } else {
      // Unknown type - shouldn't happen in normal operation
      std::cerr << "Unknown constraint type encountered during simplification."
                << std::endl;
    }
  }

  // Reconstruct geometries in original order
  auto result = std::make_unique<GeometryCollection>();

  for (std::size_t geomIdx = 0; geomIdx < collection.numGeometries();
       ++geomIdx) {
    const auto &geom = collection.geometryN(geomIdx);

    switch (geom.geometryTypeId()) {
    case TYPE_LINESTRING: {
      auto it = linestrings.find(geomIdx);
      if (it != linestrings.end()) {
        // Use the simplified linestring
        result->addGeometry(*(it->second));
      } else {
        // No matching simplified linestring, use the original
        result->addGeometry(geom);
      }
      break;
    }

    case TYPE_POLYGON: {
      auto it = polygonRings.find(geomIdx);
      if (it != polygonRings.end() && !it->second.empty()) {
        // We have at least one simplified ring
        auto polygonPtr = std::make_unique<Polygon>();

        // Set the exterior ring if it exists
        if (it->second.find(0) != it->second.end()) {
          polygonPtr->setExteriorRing(it->second[0]);
        }

        // Add interior rings
        for (const auto &ringEntry : it->second) {
          if (ringEntry.first > 0) {
            polygonPtr->addInteriorRing(ringEntry.second);
          }
        }

        result->addGeometry(*polygonPtr);
      } else {
        // No matching simplified polygon, use the original
        result->addGeometry(geom);
      }
      break;
    }

    case TYPE_MULTIPOLYGON: {
      auto it = multiPolygonRings.find(geomIdx);
      if (it != multiPolygonRings.end() && !it->second.empty()) {
        // We have at least one simplified polygon in the MultiPolygon
        auto        multiPolygonPtr = std::make_unique<MultiPolygon>();
        const auto &originalMultiPolygon =
            static_cast<const MultiPolygon &>(geom);

        // Process each polygon in the original MultiPolygon
        for (std::size_t polyIdx = 0;
             polyIdx < originalMultiPolygon.numGeometries(); ++polyIdx) {
          auto polyIt = it->second.find(polyIdx);

          if (polyIt != it->second.end() && !polyIt->second.empty()) {
            // We have at least one simplified ring for this polygon
            Polygon simplifiedPolygon;

            // Set exterior ring if it exists
            if (polyIt->second.find(0) != polyIt->second.end()) {
              simplifiedPolygon.setExteriorRing(polyIt->second[0]);
            }

            // Add interior rings
            for (const auto &ringEntry : polyIt->second) {
              if (ringEntry.first > 0) {
                simplifiedPolygon.addInteriorRing(ringEntry.second);
              }
            }

            multiPolygonPtr->addGeometry(simplifiedPolygon);
          } else {
            // No simplified version of this polygon, use the original
            multiPolygonPtr->addGeometry(
                originalMultiPolygon.polygonN(polyIdx));
          }
        }

        result->addGeometry(*multiPolygonPtr);
      } else {
        // No matching simplified MultiPolygon, use the original
        result->addGeometry(geom);
      }
      break;
    }

    case TYPE_POLYHEDRALSURFACE: {
      auto it = polyhedralSurfaceRings.find(geomIdx);
      if (it != polyhedralSurfaceRings.end() && !it->second.empty()) {
        // We have at least one simplified polygon in the PolyhedralSurface
        auto        surfacePtr = std::make_unique<PolyhedralSurface>();
        const auto &originalSurface =
            static_cast<const PolyhedralSurface &>(geom);

        // Process each polygon in the original PolyhedralSurface
        for (std::size_t polyIdx = 0; polyIdx < originalSurface.numPolygons();
             ++polyIdx) {
          auto polyIt = it->second.find(polyIdx);

          if (polyIt != it->second.end() && !polyIt->second.empty()) {
            // We have at least one simplified ring for this polygon
            Polygon simplifiedPolygon;

            // Set exterior ring if it exists
            if (polyIt->second.find(0) != polyIt->second.end()) {
              simplifiedPolygon.setExteriorRing(polyIt->second[0]);
            }

            // Add interior rings
            for (const auto &ringEntry : polyIt->second) {
              if (ringEntry.first > 0) {
                simplifiedPolygon.addInteriorRing(ringEntry.second);
              }
            }

            surfacePtr->addPolygon(simplifiedPolygon);
          } else {
            // No simplified version of this polygon, use the original
            surfacePtr->addPolygon(originalSurface.polygonN(polyIdx));
          }
        }

        result->addGeometry(*surfacePtr);
      } else {
        // No matching simplified PolyhedralSurface, use the original
        result->addGeometry(geom);
      }
      break;
    }

    default:
      // For other geometry types, just clone the original
      result->addGeometry(geom);
      break;
    }
  }

  return result;
}

/**
 * @brief Simplifies a PolyhedralSurface
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
  extractSegments(polySurface, store);

  // Create a MultiPolygon containing the same polygons and preserving
  // coordinate type
  MultiPolygon multiPolygon;

  for (size_t i = 0; i < polySurface.numPolygons(); ++i) {
    Polygon polygon = polySurface.polygonN(i);
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
    result->addPolygon(polygonCopy);
  }

  return result;
}

/**
 * @brief Wrapper for the algorithm::simplify function
 */
inline auto
simplify(const Geometry &geometry, double threshold, bool preserveTopology)
    -> std::unique_ptr<Geometry>
{
  return algorithm::simplify(geometry, threshold, preserveTopology);
}

/**
 * @brief Simplifies a GeometryCollection
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
