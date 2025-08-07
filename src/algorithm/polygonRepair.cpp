// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <algorithm>
#include <vector>

#include "SFCGAL/Exception.h"
#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/algorithm/polygonRepair.h"
#include "SFCGAL/detail/tools/Log.h"

#include <CGAL/Multipolygon_with_holes_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_repair/repair.h>
#include <CGAL/Polygon_with_holes_2.h>

using Kernel                    = SFCGAL::Kernel;
using Point_2                   = Kernel::Point_2;
using Polygon_2                 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2      = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;

// Anonymous namespace for internal utilities
namespace {

/**
 * @brief Converts CGAL Multipolygon_with_holes_2 back to SFCGAL Geometry
 */
auto
multipolygonWithHolesToSFCGAL(const Multipolygon_with_holes_2 &mp)
    -> std::unique_ptr<SFCGAL::Geometry>
{
  if (mp.number_of_polygons_with_holes() == 0) {
    return std::make_unique<SFCGAL::GeometryCollection>();
  }

  if (mp.number_of_polygons_with_holes() == 1) {
    // Single polygon - return as Polygon
    const auto &pwh    = mp.polygons_with_holes().front();
    auto        result = std::make_unique<SFCGAL::Polygon>();

    // Convert exterior ring
    if (!pwh.outer_boundary().is_empty()) {
      SFCGAL::LineString exterior;
      for (auto vertex_it = pwh.outer_boundary().vertices_begin();
           vertex_it != pwh.outer_boundary().vertices_end(); ++vertex_it) {
        exterior.addPoint(SFCGAL::Point(CGAL::to_double(vertex_it->x()),
                                        CGAL::to_double(vertex_it->y())));
      }
      // Close the ring
      if (!exterior.isEmpty()) {
        exterior.addPoint(exterior.startPoint());
        result->setExteriorRing(exterior);
      }
    }

    // Convert holes
    for (auto hole_it = pwh.holes_begin(); hole_it != pwh.holes_end();
         ++hole_it) {
      SFCGAL::LineString hole_ring;
      for (auto vertex_it = hole_it->vertices_begin();
           vertex_it != hole_it->vertices_end(); ++vertex_it) {
        hole_ring.addPoint(SFCGAL::Point(CGAL::to_double(vertex_it->x()),
                                         CGAL::to_double(vertex_it->y())));
      }
      // Close the ring
      if (!hole_ring.isEmpty()) {
        hole_ring.addPoint(hole_ring.startPoint());
        result->addInteriorRing(hole_ring);
      }
    }

    return result;
  } else {
    // Multiple polygons - return as MultiPolygon
    auto result = std::make_unique<SFCGAL::MultiPolygon>();

    for (const auto &pwh : mp.polygons_with_holes()) {
      auto polygon = std::make_unique<SFCGAL::Polygon>();

      // Convert exterior ring
      if (!pwh.outer_boundary().is_empty()) {
        SFCGAL::LineString exterior;
        for (auto vertex_it = pwh.outer_boundary().vertices_begin();
             vertex_it != pwh.outer_boundary().vertices_end(); ++vertex_it) {
          exterior.addPoint(SFCGAL::Point(CGAL::to_double(vertex_it->x()),
                                          CGAL::to_double(vertex_it->y())));
        }
        // Close the ring
        if (!exterior.isEmpty()) {
          exterior.addPoint(exterior.startPoint());
          polygon->setExteriorRing(exterior);
        }
      }

      // Convert holes
      for (auto hole_it = pwh.holes_begin(); hole_it != pwh.holes_end();
           ++hole_it) {
        SFCGAL::LineString hole_ring;
        for (auto vertex_it = hole_it->vertices_begin();
             vertex_it != hole_it->vertices_end(); ++vertex_it) {
          hole_ring.addPoint(SFCGAL::Point(CGAL::to_double(vertex_it->x()),
                                           CGAL::to_double(vertex_it->y())));
        }
        // Close the ring
        if (!hole_ring.isEmpty()) {
          hole_ring.addPoint(hole_ring.startPoint());
          polygon->addInteriorRing(hole_ring);
        }
      }

      result->addGeometry(polygon.release());
    }

    return result;
  }
}

/**
 * @brief Applies polygon repair with specified rule
 */
template <typename RepairRuleType>
auto
applyRepairWithRule(const Multipolygon_with_holes_2 &mp, RepairRuleType rule)
    -> Multipolygon_with_holes_2
{
  return CGAL::Polygon_repair::repair(mp, rule);
}

/**
 * @brief Dispatcher for repair rules
 */
auto
repairWithRule(const Multipolygon_with_holes_2     &mp,
               SFCGAL::algorithm::PolygonRepairRule rule)
    -> Multipolygon_with_holes_2
{
  switch (rule) {
  case SFCGAL::algorithm::PolygonRepairRule::EVEN_ODD_RULE:
    return applyRepairWithRule(mp, CGAL::Polygon_repair::Even_odd_rule());

#if CGAL_VERSION_MAJOR == 6 && CGAL_VERSION_MINOR >= 1
  case SFCGAL::algorithm::PolygonRepairRule::NON_ZERO_RULE:
    return applyRepairWithRule(mp, CGAL::Polygon_repair::Non_zero_rule());

  case SFCGAL::algorithm::PolygonRepairRule::UNION_RULE:
    return applyRepairWithRule(mp, CGAL::Polygon_repair::Union_rule());

  case SFCGAL::algorithm::PolygonRepairRule::INTERSECTION_RULE:
    return applyRepairWithRule(mp, CGAL::Polygon_repair::Intersection_rule());
#endif
  default:
    BOOST_THROW_EXCEPTION(SFCGAL::Exception("Unknown polygon repair rule"));
  }
}

/**
 * @brief Convert geometry to Multipolygon_with_holes_2
 */
auto
geometryToMultipolygonWithHoles(const SFCGAL::Geometry &geometry)
    -> Multipolygon_with_holes_2
{
  switch (geometry.geometryTypeId()) {
  case SFCGAL::TYPE_POLYGON: {
    const auto               &polygon = geometry.as<SFCGAL::Polygon>();
    Multipolygon_with_holes_2 mp;
    mp.add_polygon_with_holes(polygon.toPolygon_with_holes_2(true));
    return mp;
  }

  case SFCGAL::TYPE_MULTIPOLYGON: {
    const auto &multipolygon = geometry.as<SFCGAL::MultiPolygon>();
    return multipolygon.toMultipolygon_with_holes_2(true);
  }

  default:
    BOOST_THROW_EXCEPTION(SFCGAL::Exception(
        (boost::format("polygonRepair is not implemented for %s") %
         geometry.geometryType())
            .str()));
  }
}

} // anonymous namespace

namespace SFCGAL {
namespace algorithm {

auto
polygonRepair(const Geometry &geometry, PolygonRepairRule repairRule)
    -> std::unique_ptr<Geometry>
{
  // Input validation
  if (geometry.isEmpty()) {
    return std::unique_ptr<Geometry>(geometry.clone());
  }

  try {
    // Convert to CGAL Multipolygon_with_holes_2
    auto cgal_mp = geometryToMultipolygonWithHoles(geometry);

    // Apply polygon repair
    auto repaired_mp = repairWithRule(cgal_mp, repairRule);

    // Convert back to SFCGAL
    return multipolygonWithHolesToSFCGAL(repaired_mp);

  } catch (const CGAL::Precondition_exception &e) {
    BOOST_THROW_EXCEPTION(Exception(
        (boost::format("CGAL precondition error in polygon repair: %s") %
         e.what())
            .str()));
  } catch (const CGAL::Assertion_exception &e) {
    BOOST_THROW_EXCEPTION(Exception(
        (boost::format("CGAL assertion error in polygon repair: %s") % e.what())
            .str()));
  } catch (const std::exception &e) {
    BOOST_THROW_EXCEPTION(Exception(
        (boost::format("Error in polygonRepair: %s") % e.what()).str()));
  }
}

} // namespace algorithm
} // namespace SFCGAL
