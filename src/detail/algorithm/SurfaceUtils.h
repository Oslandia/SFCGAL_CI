// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_DETAIL_ALGORITHM_SURFACEUTILS_H_
#define SFCGAL_DETAIL_ALGORITHM_SURFACEUTILS_H_

#include "SFCGAL/Coordinate.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/algorithm/area.h"

#include <map>
#include <utility>
#include <vector>

namespace SFCGAL::detail::algorithm {

/**
 * @brief Normalized edge key for consistent edge comparison.
 *
 * The edge is stored with the lexicographically smaller coordinate first,
 * ensuring that edges (A,B) and (B,A) map to the same key.
 */
struct NormalizedEdge {
  Point point1;
  Point point2;

  NormalizedEdge() = default;
  NormalizedEdge(const Point &pt1, const Point &pt2)
  {
    if (comparePoints(pt1, pt2)) {
      point1 = pt1;
      point2 = pt2;
    } else {
      point1 = pt2;
      point2 = pt1;
    }
  }

  auto
  operator<(const NormalizedEdge &other) const -> bool
  {
    if (point1.x() != other.point1.x()) {
      return point1.x() < other.point1.x();
    }
    if (point1.y() != other.point1.y()) {
      return point1.y() < other.point1.y();
    }
    if (point1.z() != other.point1.z()) {
      return point1.z() < other.point1.z();
    }
    if (point2.x() != other.point2.x()) {
      return point2.x() < other.point2.x();
    }
    if (point2.y() != other.point2.y()) {
      return point2.y() < other.point2.y();
    }
    return point2.z() < other.point2.z();
  }

private:
  static auto
  comparePoints(const Point &p1, const Point &p2) -> bool
  {
    if (p1.x() != p2.x()) {
      return p1.x() < p2.x();
    }
    if (p1.y() != p2.y()) {
      return p1.y() < p2.y();
    }
    return p1.z() < p2.z();
  }
};

/**
 * @brief Oriented edge with original direction preserved.
 */
using OrientedEdge = std::pair<Point, Point>;

/**
 * @brief Result of boundary edge extraction.
 */
struct BoundaryEdgesResult {
  std::vector<OrientedEdge> edges;
  bool                      success{true};
};

/**
 * @brief Extract boundary edges from a collection of polygons.
 *
 * Boundary edges are edges that appear exactly once across all polygons.
 * The orientation of each edge is preserved from the original polygon.
 *
 * @param polygons Vector of polygons to extract boundary from.
 * @return BoundaryEdgesResult containing oriented boundary edges.
 */
inline auto
extractBoundaryEdges(const std::vector<Polygon> &polygons)
    -> BoundaryEdgesResult
{
  BoundaryEdgesResult result;

  // Count edge occurrences
  std::map<NormalizedEdge, int>          edgeCount;
  std::map<NormalizedEdge, OrientedEdge> edgeOrientation;

  for (const auto &poly : polygons) {
    for (size_t ringIdx = 0; ringIdx < poly.numRings(); ++ringIdx) {
      const LineString &ring = poly.ringN(ringIdx);
      for (size_t i = 0; i < ring.numPoints() - 1; ++i) {
        const Point   &p1 = ring.pointN(i);
        const Point   &p2 = ring.pointN(i + 1);
        NormalizedEdge key(p1, p2);
        edgeCount[key]++;
        if (edgeCount[key] == 1) {
          edgeOrientation[key] = {p1, p2};
        }
      }
    }
  }

  // Collect boundary edges (count == 1) with correct orientation
  for (const auto &[edge, count] : edgeCount) {
    if (count == 1) {
      // Find original orientation and use it directly
      const auto &origEdge = edgeOrientation[edge];
      result.edges.emplace_back(origEdge.first, origEdge.second);
    }
  }

  return result;
}

/**
 * @brief Extract boundary edges from a triangulated surface.
 *
 * @param surface The triangulated surface.
 * @return BoundaryEdgesResult containing oriented boundary edges.
 */
inline auto
extractBoundaryEdges(const TriangulatedSurface &surface) -> BoundaryEdgesResult
{
  BoundaryEdgesResult result;

  std::map<NormalizedEdge, int>          edgeCount;
  std::map<NormalizedEdge, OrientedEdge> edgeOrientation;

  for (size_t i = 0; i < surface.numPatches(); ++i) {
    const Triangle &tri = surface.patchN(i);
    for (int e = 0; e < 3; ++e) {
      const Point   &p1 = tri.vertex(e);
      const Point   &p2 = tri.vertex((e + 1) % 3);
      NormalizedEdge key(p1, p2);
      edgeCount[key]++;
      if (edgeCount[key] == 1) {
        edgeOrientation[key] = {p1, p2};
      }
    }
  }

  for (const auto &[edge, count] : edgeCount) {
    if (count == 1) {
      const auto &origEdge = edgeOrientation[edge];
      result.edges.emplace_back(origEdge.first, origEdge.second);
    }
  }

  return result;
}

/**
 * @brief Extract boundary edges from a PolyhedralSurface (triangles only).
 *
 * @param surface The polyhedral surface containing triangles.
 * @return BoundaryEdgesResult containing oriented boundary edges.
 */
inline auto
extractBoundaryEdgesFromTriangles(const PolyhedralSurface &surface)
    -> BoundaryEdgesResult
{
  BoundaryEdgesResult result;

  std::map<NormalizedEdge, int>          edgeCount;
  std::map<NormalizedEdge, OrientedEdge> edgeOrientation;

  for (size_t i = 0; i < surface.numPatches(); ++i) {
    const auto &patch = surface.patchN(i);
    if (const auto *tri = dynamic_cast<const Triangle *>(&patch)) {
      for (int e = 0; e < 3; ++e) {
        const Point   &p1 = tri->vertex(e);
        const Point   &p2 = tri->vertex((e + 1) % 3);
        NormalizedEdge key(p1, p2);
        edgeCount[key]++;
        if (edgeCount[key] == 1) {
          edgeOrientation[key] = {p1, p2};
        }
      }
    }
  }

  for (const auto &[edge, count] : edgeCount) {
    if (count == 1) {
      const auto &origEdge = edgeOrientation[edge];
      result.edges.emplace_back(origEdge.first, origEdge.second);
    }
  }

  return result;
}

/**
 * @brief Build connected loops from boundary edges.
 *
 * Groups boundary edges into connected loops (useful for handling holes).
 * Uses 2D projection (x,y) for connectivity.
 *
 * @param edges Vector of oriented boundary edges.
 * @return Vector of loops, where each loop is a vector of oriented edges.
 */
inline auto
buildConnectedLoops(std::vector<OrientedEdge> edges)
    -> std::vector<std::vector<OrientedEdge>>
{
  std::vector<std::vector<OrientedEdge>> allLoops;

  while (!edges.empty()) {
    std::vector<OrientedEdge> currentLoop;
    currentLoop.push_back(edges[0]);
    edges.erase(edges.begin());

    bool loopClosed = false;
    while (!loopClosed && !edges.empty()) {
      const Point &currentEnd = currentLoop.back().second;
      const Point &loopStart  = currentLoop.front().first;

      // Check if loop is closed (using 2D comparison)
      if (currentEnd.x() == loopStart.x() && currentEnd.y() == loopStart.y()) {
        loopClosed = true;
        break;
      }

      // Find edge starting at currentEnd
      bool foundNext = false;
      for (auto it = edges.begin(); it != edges.end(); ++it) {
        if (it->first.x() == currentEnd.x() &&
            it->first.y() == currentEnd.y()) {
          currentLoop.push_back(*it);
          edges.erase(it);
          foundNext = true;
          break;
        }
      }
      if (!foundNext) {
        break;
      }
    }

    if (!currentLoop.empty()) {
      allLoops.push_back(std::move(currentLoop));
    }
  }

  return allLoops;
}

/**
 * @brief Find exterior loop index (largest absolute area).
 *
 * @param loops Vector of boundary loops.
 * @return Index of the exterior loop.
 */
inline auto
findExteriorLoopIndex(const std::vector<std::vector<OrientedEdge>> &loops)
    -> size_t
{
  if (loops.size() <= 1) {
    return 0;
  }

  size_t     exteriorIdx = 0;
  Kernel::FT maxArea(0);

  for (size_t i = 0; i < loops.size(); ++i) {
    // Convert loop edges to LineString for area calculation
    LineString ring;
    for (const auto &[start, end] : loops[i]) {
      ring.addPoint(Point(start.x(), start.y(), 0));
    }
    if (!ring.isEmpty()) {
      ring.addPoint(ring.pointN(0)); // Close the ring
    }

    Kernel::FT loopArea = CGAL::abs(SFCGAL::algorithm::signedArea(ring));
    if (loopArea > maxArea) {
      maxArea     = loopArea;
      exteriorIdx = i;
    }
  }

  return exteriorIdx;
}

/**
 * @brief Create vertical wall triangles from a boundary edge.
 *
 * Creates wall triangles connecting a roof edge to the base (z=0).
 * The wall orientation is consistent with the roof edge direction.
 *
 * @param roofPt1 First point of the roof edge.
 * @param roofPt2 Second point of the roof edge.
 * @param baseZ The z-coordinate of the base (default 0).
 * @return Vector of triangles forming the wall (0-2 triangles).
 */
inline auto
createWallTriangles(const Point &roofPt1, const Point &roofPt2,
                    Kernel::FT baseZ = Kernel::FT(0)) -> std::vector<Triangle>
{
  std::vector<Triangle> walls;

  Point basePt1(roofPt1.x(), roofPt1.y(), baseZ);
  Point basePt2(roofPt2.x(), roofPt2.y(), baseZ);

  bool hasHeight1 = roofPt1.z() != baseZ;
  bool hasHeight2 = roofPt2.z() != baseZ;

  if (!hasHeight1 && !hasHeight2) {
    return walls; // No wall needed
  }

  if (hasHeight1 && hasHeight2) {
    // Two triangles for quadrilateral wall
    walls.emplace_back(basePt2, basePt1, roofPt1);
    walls.emplace_back(basePt2, roofPt1, roofPt2);
  } else if (hasHeight1) {
    walls.emplace_back(basePt2, basePt1, roofPt1);
  } else {
    walls.emplace_back(basePt1, roofPt2, basePt2);
  }

  return walls;
}

} // namespace SFCGAL::detail::algorithm

#endif // SFCGAL_DETAIL_ALGORITHM_SURFACEUTILS_H_
