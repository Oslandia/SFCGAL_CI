// Copyright (c) 2025-2025, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/ShapeRegularization.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/Triangle.h"

#include "SFCGAL/Kernel.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/detail/SegmentStore.h"

#include <CGAL/Shape_regularization.h>
#ifdef CGAL_USE_OSQP
#include <CGAL/OSQP_quadratic_program_traits.h>
#endif

#include <algorithm>
#include <vector>

namespace SFCGAL::algorithm {

namespace {

using Seg2 = Kernel::Segment_2;
using Pt2  = Kernel::Point_2;

static void
collect_segments_2d(const Geometry &geom, std::vector<Seg2> &out)
{
  switch (geom.geometryTypeId()) {
  case TYPE_LINESTRING: {
    const auto &ls = geom.as<LineString>();
    if (ls.numPoints() < 2) {
      return;
    }
    for (size_t i = 1; i < ls.numPoints(); ++i) {
      const Pt2 p0 = ls.pointN(i - 1).toPoint_2();
      const Pt2 p1 = ls.pointN(i).toPoint_2();
      if (p0 != p1) {
        out.emplace_back(p0, p1);
      }
    }
    break;
  }

  case TYPE_POLYGON: {
    const auto &pg = geom.as<Polygon>();
    const auto &er = pg.exteriorRing();
    if (er.numPoints() >= 2) {
      for (size_t i = 1; i < er.numPoints(); ++i) {
        const Pt2 p0 = er.pointN(i - 1).toPoint_2();
        const Pt2 p1 = er.pointN(i).toPoint_2();
        if (p0 != p1) {
          out.emplace_back(p0, p1);
        }
      }
    }
    for (size_t r = 0; r < pg.numInteriorRings(); ++r) {
      const auto &ir = pg.interiorRingN(r);
      if (ir.numPoints() >= 2) {
        for (size_t i = 1; i < ir.numPoints(); ++i) {
          const Pt2 p0 = ir.pointN(i - 1).toPoint_2();
          const Pt2 p1 = ir.pointN(i).toPoint_2();
          if (p0 != p1) {
            out.emplace_back(p0, p1);
          }
        }
      }
    }
    break;
  }

  case TYPE_MULTILINESTRING: {
    const auto &mls = geom.as<MultiLineString>();
    for (size_t i = 0; i < mls.numGeometries(); ++i) {
      collect_segments_2d(mls.lineStringN(i), out);
    }
    break;
  }

  case TYPE_MULTIPOLYGON: {
    const auto &mp = geom.as<MultiPolygon>();
    for (size_t i = 0; i < mp.numGeometries(); ++i) {
      collect_segments_2d(mp.polygonN(i), out);
    }
    break;
  }

  case TYPE_GEOMETRYCOLLECTION: {
    const auto &gc = geom.as<GeometryCollection>();
    for (size_t i = 0; i < gc.numGeometries(); ++i) {
      collect_segments_2d(gc.geometryN(i), out);
    }
    break;
  }

  case TYPE_TRIANGLE: {
    const auto &tri = geom.as<Triangle>();
    LineString    ring = tri.toPolygon().exteriorRing();
    if (ring.numPoints() >= 2) {
      for (size_t i = 1; i < ring.numPoints(); ++i) {
        const Pt2 p0 = ring.pointN(i - 1).toPoint_2();
        const Pt2 p1 = ring.pointN(i).toPoint_2();
        if (p0 != p1) {
          out.emplace_back(p0, p1);
        }
      }
    }
    break;
  }

  default:
    break;
  }
}

// Build a MultiLineString of 2-point segments from regularized CGAL segments.
static std::unique_ptr<Geometry>
segments_to_multilinestring(const std::vector<Seg2> &segments,
                            const Geometry         &source)
{
  auto result = std::make_unique<MultiLineString>();

  // Build a store for Z/M interpolation from the source geometry
  detail::SegmentStore store;
  store.extractSegments(source);

  const CoordinateType dim = source.getCoordinateType();

  for (const auto &s : segments) {
    const double x0 = CGAL::to_double(s.source().x());
    const double y0 = CGAL::to_double(s.source().y());
    const double x1 = CGAL::to_double(s.target().x());
    const double y1 = CGAL::to_double(s.target().y());

    LineString seg;
    seg.addPoint(store.createPoint(x0, y0, dim));
    seg.addPoint(store.createPoint(x1, y1, dim));
    result->addGeometry(new LineString(seg));
  }

  return result;
}

// Convert a LineString to a vector of 2D points.
static std::vector<Pt2>
linestring_to_points_2d(const LineString &ls)
{
  std::vector<Pt2> pts;
  pts.reserve(ls.numPoints());
  for (size_t i = 0; i < ls.numPoints(); ++i) {
    pts.emplace_back(ls.pointN(i).toPoint_2());
  }
  return pts;
}

// Rebuild a LineString from 2D points, interpolating Z/M from a SegmentStore.
static LineString
points2d_to_linestring(const std::vector<Pt2> &pts,
                       const detail::SegmentStore &store, CoordinateType dim)
{
  LineString ret;
  ret.reserve(pts.size());
  for (const auto &p : pts) {
    const double x = CGAL::to_double(p.x());
    const double y = CGAL::to_double(p.y());
    ret.addPoint(store.createPoint(x, y, dim));
  }
  return ret;
}

static void
regularize_contour_via_segments(std::vector<Pt2> &contour, bool apply_angles = true)
{
  if (contour.size() < 4) {
    return;
  }
  
  std::vector<Seg2> segments;
  segments.reserve(contour.size() - 1);
  
  for (size_t i = 0; i < contour.size() - 1; ++i) {
    const Seg2 seg(contour[i], contour[i + 1]);
    if (seg.source() != seg.target()) {
      segments.push_back(seg);
    }
  }
  
  if (segments.empty()) {
    return;
  }
  
#if defined(CGAL_USE_OSQP) && CGAL_VERSION_NR >= 1060000000
  using QPSolver = CGAL::OSQP_quadratic_program_traits;
  
  if (apply_angles) {
    CGAL::Shape_regularization::Segments::
        regularize_angles<Kernel, std::vector<Seg2>, QPSolver>(
            segments,
            CGAL::Shape_regularization::Segments::delaunay_neighbor_query_2(segments),
            CGAL::Shape_regularization::Segments::make_angle_regularization_2(segments),
            QPSolver());
  }
#endif
  
  if (segments.size() >= contour.size() - 1) {
    for (size_t i = 0; i < segments.size() && i < contour.size() - 1; ++i) {
      contour[i] = segments[i].source();
      if (i == segments.size() - 1) {
        contour[i + 1] = segments[i].target();
      }
    }
  }
}

static bool
regularize_closed_contour_points(std::vector<Pt2> &contour,
                                 ShapeRegularization::DirectionEstimator estimator)
{
  using namespace CGAL::Shape_regularization::Contours;

  if (contour.size() < 4) {
    return true;
  }

  // Remove duplicate consecutive points
  std::vector<Pt2> clean_contour;
  clean_contour.reserve(contour.size());
  clean_contour.push_back(contour[0]);
  
  for (size_t i = 1; i < contour.size(); ++i) {
    if (contour[i] != contour[i-1]) {
      clean_contour.push_back(contour[i]);
    }
  }
  
  if (clean_contour.size() < 4) {
    return true;
  }
  
  // Check for sufficient area to avoid degenerate contours
  double total_area = 0.0;
  for (size_t i = 2; i < clean_contour.size(); ++i) {
    const auto area = CGAL::area(clean_contour[0], clean_contour[i-1], clean_contour[i]);
    total_area += CGAL::to_double(CGAL::abs(area));
  }
  
  if (total_area < 0.0001) {
    return true;
  }
  
  // Check for numerical stability
  for (const auto &p : clean_contour) {
    const double x = CGAL::to_double(p.x());
    const double y = CGAL::to_double(p.y());
    if (!std::isfinite(x) || !std::isfinite(y) ||
        std::abs(x) > 1e10 || std::abs(y) > 1e10) {
      return true;
    }
  }

  // Ensure proper closure
  if (clean_contour.front() != clean_contour.back()) {
    if (clean_contour.back() == clean_contour[0]) {
      clean_contour.pop_back();
    }
    clean_contour.push_back(clean_contour[0]);
  }

  // Use segment-based regularization to avoid CGAL stability issues
  // CGAL's regularize_closed_contour has known SIGFPE bugs
  regularize_contour_via_segments(clean_contour, true);
  contour = std::move(clean_contour);
  
  (void)estimator; // Suppress unused parameter warning
  
  return true;
}

static bool
regularize_open_contour_points(std::vector<Pt2> &contour,
                               ShapeRegularization::DirectionEstimator estimator)
{
  using namespace CGAL::Shape_regularization::Contours;

  if (contour.size() < 2) {
    return true; // nothing to do
  }

  try {
    if (estimator == ShapeRegularization::DirectionEstimator::MULTIPLE) {
      std::vector<Pt2> out;
      out.reserve(contour.size());
      Multiple_directions_2<Kernel, std::vector<Pt2>> dirs(contour, false);
      regularize_open_contour(contour, dirs, std::back_inserter(out));
      if (out.size() >= 2) {
        contour = std::move(out);
        return true;
      }
    } else {
      std::vector<Pt2> out;
      out.reserve(contour.size());
      Longest_direction_2<Kernel, std::vector<Pt2>> dirs(contour, false);
      regularize_open_contour(contour, dirs, std::back_inserter(out));
      if (out.size() >= 2) {
        contour = std::move(out);
        return true;
      }
    }
  } catch (const std::exception &e) {
    // CGAL threw an exception - return false to signal failure
    return false;
  } catch (...) {
    // Catch any other exceptions
    return false;
  }
  
  return false;
}

} // namespace

// ----------------------------------------------------------------------------------
// Public API
// ----------------------------------------------------------------------------------

auto
ShapeRegularization::regularizeSegments(const Geometry &geometry,
                                        bool           applyAngles,
                                        bool           applyOffsets,
                                        double         angleBoundDeg,
                                        double         offsetBound,
                                        double         angleTolerance,
                                        NeighborQuery  query)
    -> std::unique_ptr<Geometry>
{
  // Validity check in 2D unless 3D-specific is needed.
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(geometry);

  auto result = regularizeSegments(geometry, applyAngles, applyOffsets,
                                   angleBoundDeg, offsetBound, angleTolerance,
                                   query, NoValidityCheck());
  propagateValidityFlag(*result, true);
  return result;
}

auto
ShapeRegularization::regularizeSegments(const Geometry &geometry,
                                        bool           applyAngles,
                                        bool           applyOffsets,
                                        double         angleBoundDeg,
                                        double         offsetBound,
                                        double         angleTolerance,
                                        NeighborQuery  query,
                                        NoValidityCheck /*unused*/)
    -> std::unique_ptr<Geometry>
{
  if (geometry.isEmpty()) {
    return std::unique_ptr<Geometry>(geometry.clone());
  }

  std::vector<Seg2> segments;
  collect_segments_2d(geometry, segments);

  if (segments.empty()) {
    return std::unique_ptr<Geometry>(geometry.clone());
  }

  using namespace CGAL::Shape_regularization::Segments;
  using Segment_map = CGAL::Identity_property_map<Seg2>;

  // Neighbor query over segments
  Delaunay_neighbor_query_2<Kernel, std::vector<Seg2>, Segment_map> nq(
      segments);

#ifdef CGAL_USE_OSQP
  using FT = Kernel::FT;
  CGAL::OSQP_quadratic_program_traits<FT> qp;

  // Apply selected regularizations (in-place on segments)
  if (applyAngles) {
    Angle_regularization_2<Kernel, std::vector<Seg2>, Segment_map> angle(
        segments, CGAL::parameters::maximum_angle(angleBoundDeg));
    regularize_segments(segments, nq, angle, qp);
  }

  if (applyOffsets) {
    // Compute parallel groups to guide the offset regularization
    std::vector<std::vector<std::size_t>> par_groups;
    parallel_groups(segments, std::back_inserter(par_groups));

    Offset_regularization_2<Kernel, std::vector<Seg2>, Segment_map> offset(
        segments, CGAL::parameters::maximum_offset(offsetBound));

    nq.clear();
    for (const auto &g : par_groups) {
      nq.add_group(g);
      offset.add_group(g);
    }
    regularize_segments(segments, nq, offset, qp);
  }
#else
  // Without OSQP support in CGAL, leave segments unchanged.
  (void)applyAngles;
  (void)applyOffsets;
  (void)angleBoundDeg;
  (void)offsetBound;
#endif

  return segments_to_multilinestring(segments, geometry);
}

auto
ShapeRegularization::regularizeClosedContour(const Geometry &geometry,
                                             DirectionEstimator estimator,
                                             double angleToleranceDeg,
                                             double minSegmentLength)
    -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(geometry);
  auto ret = regularizeClosedContour(geometry, estimator, angleToleranceDeg,
                                     minSegmentLength, NoValidityCheck());
  propagateValidityFlag(*ret, true);
  return ret;
}

auto
ShapeRegularization::regularizeClosedContour(const Geometry &geometry,
                                             DirectionEstimator estimator,
                                             double angleToleranceDeg,
                                             double minSegmentLength,
                                             NoValidityCheck /*unused*/)
    -> std::unique_ptr<Geometry>
{
  if (geometry.isEmpty()) {
    return std::unique_ptr<Geometry>(geometry.clone());
  }

  detail::SegmentStore store;
  store.extractSegments(geometry);
  const CoordinateType dim = geometry.getCoordinateType();

  switch (geometry.geometryTypeId()) {
  case TYPE_LINESTRING: {
    const auto &ls = geometry.as<LineString>();
    // Expect closed ring
    if (ls.numPoints() < 4 || ls.startPoint() != ls.endPoint()) {
      return std::unique_ptr<Geometry>(geometry.clone());
    }
    auto pts = linestring_to_points_2d(ls);
    if (!regularize_closed_contour_points(pts, estimator)) {
      // Regularization failed, return original geometry
      return std::unique_ptr<Geometry>(geometry.clone());
    }
    LineString out = points2d_to_linestring(pts, store, dim);
    return std::unique_ptr<Geometry>(new LineString(out));
  }

  case TYPE_POLYGON: {
    const auto &poly = geometry.as<Polygon>();
    auto        out  = std::make_unique<Polygon>();

    // Exterior ring
    {
      auto er_pts = linestring_to_points_2d(poly.exteriorRing());
      if (!regularize_closed_contour_points(er_pts, estimator)) {
        // Regularization failed on exterior ring, return original
        return std::unique_ptr<Geometry>(geometry.clone());
      }
      out->setExteriorRing(points2d_to_linestring(er_pts, store, dim).clone());
    }

    // Interior rings
    for (size_t i = 0; i < poly.numInteriorRings(); ++i) {
      auto ir_pts = linestring_to_points_2d(poly.interiorRingN(i));
      if (!regularize_closed_contour_points(ir_pts, estimator)) {
        // Regularization failed on interior ring, return original
        return std::unique_ptr<Geometry>(geometry.clone());
      }
      out->addInteriorRing(points2d_to_linestring(ir_pts, store, dim).clone());
    }

    return out;
  }

  case TYPE_MULTIPOLYGON: {
    const auto &mp = geometry.as<MultiPolygon>();
    auto        out = std::make_unique<MultiPolygon>();
    for (size_t i = 0; i < mp.numGeometries(); ++i) {
      auto p = regularizeClosedContour(mp.polygonN(i), estimator, 
                                       angleToleranceDeg, minSegmentLength,
                                       NoValidityCheck());
      out->addGeometry(p.release());
    }
    return out;
  }

  default:
    return std::unique_ptr<Geometry>(geometry.clone());
  }
}

auto
ShapeRegularization::regularizeOpenContour(const Geometry &geometry,
                                           DirectionEstimator estimator,
                                           double angleToleranceDeg,
                                           double minSegmentLength)
    -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(geometry);
  auto ret = regularizeOpenContour(geometry, estimator, angleToleranceDeg,
                                   minSegmentLength, NoValidityCheck());
  propagateValidityFlag(*ret, true);
  return ret;
}

auto
ShapeRegularization::regularizeOpenContour(const Geometry &geometry,
                                           DirectionEstimator estimator,
                                           double angleToleranceDeg,
                                           double minSegmentLength,
                                           NoValidityCheck /*unused*/)
    -> std::unique_ptr<Geometry>
{
  if (geometry.isEmpty()) {
    return std::unique_ptr<Geometry>(geometry.clone());
  }

  detail::SegmentStore store;
  store.extractSegments(geometry);
  const CoordinateType dim = geometry.getCoordinateType();

  switch (geometry.geometryTypeId()) {
  case TYPE_LINESTRING: {
    const auto &ls  = geometry.as<LineString>();
    const bool  is_closed = (ls.numPoints() >= 2 && ls.startPoint() == ls.endPoint());
    if (is_closed) {
      // Defer to closed contour variant (which handles errors internally)
      return regularizeClosedContour(geometry, estimator, angleToleranceDeg,
                                     minSegmentLength, NoValidityCheck());
    }

    auto pts = linestring_to_points_2d(ls);
    if (!regularize_open_contour_points(pts, estimator)) {
      // Regularization failed, return original geometry
      return std::unique_ptr<Geometry>(geometry.clone());
    }
    LineString out = points2d_to_linestring(pts, store, dim);
    return std::unique_ptr<Geometry>(new LineString(out));
  }

  case TYPE_MULTILINESTRING: {
    const auto &mls = geometry.as<MultiLineString>();
    auto        out = std::make_unique<MultiLineString>();
    for (size_t i = 0; i < mls.numGeometries(); ++i) {
      auto geom = regularizeOpenContour(mls.lineStringN(i), estimator,
                                        angleToleranceDeg, minSegmentLength,
                                        NoValidityCheck());
      out->addGeometry(geom.release());
    }
    return out;
  }

  default:
    return std::unique_ptr<Geometry>(geometry.clone());
  }
}

} // namespace SFCGAL::algorithm
