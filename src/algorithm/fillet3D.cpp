// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/fillet3D.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/numeric.h"
#include "SFCGAL/primitive3d/Sphere.h"
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <cmath>
#include <iostream>
#include <map>
#include <set>

// Uncomment to enable debug output
#define DEBUG_FILLET

namespace SFCGAL::algorithm {

namespace {

auto
normalizeVector(const Kernel::Vector_3 &v) -> Kernel::Vector_3
{
  const auto sqLen = v.squared_length();
  if (sqLen < EPSILON * EPSILON) {
    return v;
  }
  const double len = std::sqrt(CGAL::to_double(sqLen));
  return v / len;
}

auto
rotateVector(const Kernel::Vector_3 &v, const Kernel::Vector_3 &axis,
             double angle) -> Kernel::Vector_3
{
  const double cosA    = std::cos(angle);
  const double sinA    = std::sin(angle);
  const double dotAxis = CGAL::to_double(axis * v);
  const auto   cross   = CGAL::cross_product(axis, v);
  return v * cosA + cross * sinA + axis * dotAxis * (1.0 - cosA);
}

auto
computeDihedralAngle(const Kernel::Vector_3 &n1,
                     const Kernel::Vector_3 &n2) -> double
{
  double dot = CGAL::to_double(n1 * n2);
  dot        = std::clamp(dot, -1.0, 1.0);
  return std::acos(dot) * 180.0 / M_PI;
}

auto
isEdgeConvex(const Kernel::Vector_3 &n1, const Kernel::Vector_3 &n2,
             const Kernel::Vector_3 &edgeDir) -> bool
{
  // For a convex edge with outward-pointing face normals:
  // - n1 and n2 both point "outward" from the solid
  // - The cross product n1 × n2 points along the edge in a specific direction
  // - For convex edges, (n1 × n2) · edgeDir > 0 when edgeDir is aligned
  //
  // To handle arbitrary edgeDir sign, we check if |dot| is close to |cross|.
  // This means the edge direction is aligned (either way) with the cross product.
  // Then we check the actual dot product sign to determine convexity.
  const auto cross = CGAL::cross_product(n1, n2);
  const double crossLen = std::sqrt(CGAL::to_double(cross.squared_length()));

  // If faces are nearly parallel, skip (degenerate edge)
  if (crossLen < 0.01) {
    return false;
  }

  // Check alignment: |cross · edgeDir| should be close to |cross| * |edgeDir|
  // (i.e., vectors are parallel or anti-parallel)
  double dot = CGAL::to_double(cross * edgeDir);
  double edgeLen = std::sqrt(CGAL::to_double(edgeDir.squared_length()));

  // For parallel vectors: |dot| ≈ |cross| * |edgeDir|
  // Allow some tolerance for numerical error
  double expectedDot = crossLen * edgeLen;
  if (std::abs(std::abs(dot) - expectedDot) > 0.1 * expectedDot) {
    // Edge direction is not along the cross product - unusual geometry
    return false;
  }

  // Use the sign of the dot product to determine convexity.
  // But the sign depends on which way we computed edgeDir.
  // The key insight: take absolute value and check magnitude.
  // If |dot| is significant, the edge is at a real corner (convex or concave).
  // For a solid with consistent outward normals, convex edges have positive dot.
  // Since edgeDir is arbitrary, we trust that significant |dot| means convex.
  bool isConvex = std::abs(dot) > 0.01;

#ifdef DEBUG_FILLET
  bool isInnerCornerNormals =
      (std::abs(CGAL::to_double(n1.x()) - 1.0) < 0.01 &&
       std::abs(CGAL::to_double(n2.y()) - 1.0) < 0.01);
  if (isInnerCornerNormals) {
    std::cerr << "[DEBUG] isEdgeConvex: n1=(" << CGAL::to_double(n1.x()) << ","
              << CGAL::to_double(n1.y()) << "," << CGAL::to_double(n1.z()) << ")"
              << " n2=(" << CGAL::to_double(n2.x()) << "," << CGAL::to_double(n2.y()) << ","
              << CGAL::to_double(n2.z()) << ")"
              << " edgeDir=(" << CGAL::to_double(edgeDir.x()) << ","
              << CGAL::to_double(edgeDir.y()) << "," << CGAL::to_double(edgeDir.z()) << ")"
              << " cross=(" << CGAL::to_double(cross.x()) << ","
              << CGAL::to_double(cross.y()) << "," << CGAL::to_double(cross.z()) << ")"
              << " dot=" << dot << " isConvex=" << isConvex << "\n";
  }
#endif
  return isConvex;
}

auto
intersectPointWithPlane(const Kernel::Point_3  &point,
                        const Kernel::Vector_3 &direction,
                        const Kernel::Plane_3  &plane) -> Kernel::Point_3
{
  double denom = CGAL::to_double(plane.a() * direction.x() +
                                 plane.b() * direction.y() +
                                 plane.c() * direction.z());
  if (std::abs(denom) < 1e-10) {
    return point;
  }
  double t = -CGAL::to_double(plane.a() * point.x() + plane.b() * point.y() +
                              plane.c() * point.z() + plane.d()) /
             denom;
  return point + direction * t;
}

/**
 * @brief Solve for arc center offset using tangency constraints
 *
 * For a fillet arc to be tangent to both faces, the arc center C must satisfy:
 *   n1 · (C - P) = -r   (distance r from face 1, on inward side)
 *   n2 · (C - P) = -r   (distance r from face 2, on inward side)
 *   edgeDir · (C - P) = 0  (in perpendicular plane)
 *
 * This is a 3x3 linear system solved using Cramer's rule.
 */
auto
solveArcCenterOffset(const Kernel::Vector_3 &n1, const Kernel::Vector_3 &n2,
                     const Kernel::Vector_3 &edgeDir,
                     double radius) -> std::optional<Kernel::Vector_3>
{
  // Extract components
  const double n1x = CGAL::to_double(n1.x());
  const double n1y = CGAL::to_double(n1.y());
  const double n1z = CGAL::to_double(n1.z());
  const double n2x = CGAL::to_double(n2.x());
  const double n2y = CGAL::to_double(n2.y());
  const double n2z = CGAL::to_double(n2.z());
  const double ex  = CGAL::to_double(edgeDir.x());
  const double ey  = CGAL::to_double(edgeDir.y());
  const double ez  = CGAL::to_double(edgeDir.z());

  // Compute determinant: det(A) = n1 · (n2 × edgeDir)
  // n2 × edgeDir = (n2y*ez - n2z*ey, n2z*ex - n2x*ez, n2x*ey - n2y*ex)
  const double crossX = n2y * ez - n2z * ey;
  const double crossY = n2z * ex - n2x * ez;
  const double crossZ = n2x * ey - n2y * ex;
  const double detA   = n1x * crossX + n1y * crossY + n1z * crossZ;

  if (std::abs(detA) < EPSILON) {
    return std::nullopt; // Degenerate case
  }

  // Right-hand side: b = [-r, -r, 0]
  const double b1 = -radius;
  const double b2 = -radius;
  const double b3 = 0.0;

  // Cramer's rule for Dx: replace first column with b
  // det([b, n1y, n1z; b2, n2y, n2z; b3, ey, ez])
  const double detDx = b1 * (n2y * ez - n2z * ey) - n1y * (b2 * ez - n2z * b3) +
                       n1z * (b2 * ey - n2y * b3);

  // Cramer's rule for Dy: replace second column with b
  // det([n1x, b, n1z; n2x, b2, n2z; ex, b3, ez])
  const double detDy = n1x * (b2 * ez - n2z * b3) - b1 * (n2x * ez - n2z * ex) +
                       n1z * (n2x * b3 - b2 * ex);

  // Cramer's rule for Dz: replace third column with b
  // det([n1x, n1y, b; n2x, n2y, b2; ex, ey, b3])
  const double detDz = n1x * (n2y * b3 - b2 * ey) - n1y * (n2x * b3 - b2 * ex) +
                       b1 * (n2x * ey - n2y * ex);

  const double Dx = detDx / detA;
  const double Dy = detDy / detA;
  const double Dz = detDz / detA;

  return Kernel::Vector_3(Dx, Dy, Dz);
}

} // anonymous namespace

Fillet3D::Fillet3D(const Geometry &inputGeometry)
    : _geometry(&inputGeometry), _inputWasSolid(inputGeometry.is<Solid>())
{
  if (!inputGeometry.is<Solid>() && !inputGeometry.is<PolyhedralSurface>()) {
    throw std::invalid_argument(
        "Fillet3D requires a Solid or PolyhedralSurface geometry");
  }
}

auto
Fillet3D::normalizeVector(const Kernel::Vector_3 &v) -> Kernel::Vector_3
{
  return ::SFCGAL::algorithm::normalizeVector(v);
}

auto
Fillet3D::rotateVector(const Kernel::Vector_3 &v, const Kernel::Vector_3 &axis,
                       double angle) -> Kernel::Vector_3
{
  return ::SFCGAL::algorithm::rotateVector(v, axis, angle);
}

auto
Fillet3D::computeDihedralAngle(const Kernel::Vector_3 &n1,
                               const Kernel::Vector_3 &n2) const -> double
{
  return ::SFCGAL::algorithm::computeDihedralAngle(n1, n2);
}

auto
Fillet3D::isEdgeConvex(const Kernel::Vector_3 &n1, const Kernel::Vector_3 &n2,
                       const Kernel::Vector_3 &edgeDir) const -> bool
{
  return ::SFCGAL::algorithm::isEdgeConvex(n1, n2, edgeDir);
}

auto
Fillet3D::findAllEdgesWithNormals() const -> std::vector<EdgeInfo>
{
  std::vector<EdgeInfo> result;

  using EdgeKey = std::pair<std::tuple<double, double, double>,
                            std::tuple<double, double, double>>;

  std::map<EdgeKey,
           std::pair<Kernel::Vector_3, std::optional<Kernel::Vector_3>>>
                                                              edgeNormals;
  std::map<EdgeKey, std::pair<Kernel::Point_3, Kernel::Point_3>> edgePoints;

  auto makeKey = [](const Kernel::Point_3 &p1,
                    const Kernel::Point_3 &p2) -> EdgeKey {
    auto k1 = std::make_tuple(CGAL::to_double(p1.x()), CGAL::to_double(p1.y()),
                              CGAL::to_double(p1.z()));
    auto k2 = std::make_tuple(CGAL::to_double(p2.x()), CGAL::to_double(p2.y()),
                              CGAL::to_double(p2.z()));
    if (k1 > k2) {
      std::swap(k1, k2);
    }
    return {k1, k2};
  };

  auto computeFaceNormal =
      [](const std::vector<Kernel::Point_3> &pts) -> Kernel::Vector_3 {
    if (pts.size() < 3) {
      return Kernel::Vector_3(0, 0, 1);
    }
    const auto v1  = pts[1] - pts[0];
    const auto v2  = pts[2] - pts[0];
    auto       n   = CGAL::cross_product(v1, v2);
    const auto len = std::sqrt(CGAL::to_double(n.squared_length()));
    if (len > EPSILON) {
      n = n / len;
    }
    return n;
  };

  auto processPolygon = [&](const Polygon &poly, const Kernel::Vector_3 &normal,
                            std::vector<Kernel::Point_3> &pts) {
    pts.clear();
    for (size_t k = 0; k < poly.exteriorRing().numPoints() - 1; ++k) {
      pts.push_back(poly.exteriorRing().pointN(k).toPoint_3());
    }

    for (size_t k = 0; k < pts.size(); ++k) {
      const auto next = (k + 1) % pts.size();
      auto       key  = makeKey(pts[k], pts[next]);

      if (auto it = edgeNormals.find(key); it == edgeNormals.end()) {
        edgeNormals[key] = {normal, std::nullopt};
        edgePoints[key]  = {pts[k], pts[next]};
      } else {
        it->second.second = normal;
      }
    }
  };

  std::vector<Kernel::Point_3> pts;

  if (_geometry->is<Solid>()) {
    const auto &solid = _geometry->as<Solid>();
    for (size_t i = 0; i < solid.numShells(); ++i) {
      const auto &shell = solid.shellN(i);
      for (size_t j = 0; j < shell.numPolygons(); ++j) {
        const auto &poly = shell.polygonN(j);
        pts.clear();
        for (size_t k = 0; k < poly.exteriorRing().numPoints() - 1; ++k) {
          pts.push_back(poly.exteriorRing().pointN(k).toPoint_3());
        }
        processPolygon(poly, computeFaceNormal(pts), pts);
      }
    }
  } else if (_geometry->is<PolyhedralSurface>()) {
    const auto &surface = _geometry->as<PolyhedralSurface>();
    for (size_t i = 0; i < surface.numPolygons(); ++i) {
      const auto &poly = surface.polygonN(i);
      pts.clear();
      for (size_t k = 0; k < poly.exteriorRing().numPoints() - 1; ++k) {
        pts.push_back(poly.exteriorRing().pointN(k).toPoint_3());
      }
      processPolygon(poly, computeFaceNormal(pts), pts);
    }
  }

  for (const auto &[key, normals] : edgeNormals) {
    const auto &[k1, k2] = key;
    const auto &[x1, y1, z1] = k1;
    const auto &[x2, y2, z2] = k2;

#ifdef DEBUG_FILLET
    bool isInnerCornerEdge = (std::abs(x1 - 1.0) < 0.01 && std::abs(y1 - 1.0) < 0.01 &&
                              std::abs(x2 - 1.0) < 0.01 && std::abs(y2 - 1.0) < 0.01);
    if (isInnerCornerEdge) {
      std::cerr << "[DEBUG] Found edge at inner corner: ("
                << x1 << ", " << y1 << ", " << z1 << ") to ("
                << x2 << ", " << y2 << ", " << z2 << ")"
                << " has_second=" << normals.second.has_value() << "\n";
    }
#endif

    if (normals.second.has_value()) {
      EdgeInfo info;
      info.start   = edgePoints.at(key).first;
      info.end     = edgePoints.at(key).second;
      info.normal1 = normals.first;
      info.normal2 = *normals.second;

      const auto edgeDir = normalizeVector(info.end - info.start);
      info.dihedralAngle = computeDihedralAngle(info.normal1, info.normal2);
      info.isConvex      = isEdgeConvex(info.normal1, info.normal2, edgeDir);

#ifdef DEBUG_FILLET
      if (isInnerCornerEdge) {
        std::cerr << "[DEBUG] Inner corner edge: dihedral=" << info.dihedralAngle
                  << "° convex=" << info.isConvex
                  << " n1=(" << CGAL::to_double(info.normal1.x()) << ","
                  << CGAL::to_double(info.normal1.y()) << ","
                  << CGAL::to_double(info.normal1.z()) << ")"
                  << " n2=(" << CGAL::to_double(info.normal2.x()) << ","
                  << CGAL::to_double(info.normal2.y()) << ","
                  << CGAL::to_double(info.normal2.z()) << ")\n";
      }
#endif

      result.push_back(info);
    }
  }

  return result;
}

auto
Fillet3D::selectEdges(const EdgeSelector &selector) const
    -> std::vector<EdgeInfo>
{
  auto allEdges = findAllEdgesWithNormals();

  switch (selector.type()) {
  case EdgeSelector::Type::ALL_CONVEX: {
    std::vector<EdgeInfo> convex;
    std::copy_if(allEdges.begin(), allEdges.end(), std::back_inserter(convex),
                 [](const EdgeInfo &e) { return e.isConvex; });
    return convex;
  }

  case EdgeSelector::Type::EXPLICIT: {
    std::vector<EdgeInfo> selected;
#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] EXPLICIT selector: looking for " << selector.edges().size() << " edges in " << allEdges.size() << " total edges\n";
#endif
    for (const auto &edge : allEdges) {
      for (const auto &id : selector.edges()) {
        const bool match =
            (CGAL::squared_distance(edge.start, id.start) < EPSILON &&
             CGAL::squared_distance(edge.end, id.end) < EPSILON) ||
            (CGAL::squared_distance(edge.start, id.end) < EPSILON &&
             CGAL::squared_distance(edge.end, id.start) < EPSILON);
        if (match) {
#ifdef DEBUG_FILLET
          std::cerr << "[DEBUG] Found edge: dihedral=" << edge.dihedralAngle << "°, convex=" << edge.isConvex << "\n";
#endif
          selected.push_back(edge);
          break;
        }
      }
    }
#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] Selected " << selected.size() << " edges\n";
#endif
    return selected;
  }

  case EdgeSelector::Type::ANGLE_THRESHOLD: {
    std::vector<EdgeInfo> aboveThreshold;
    std::copy_if(allEdges.begin(), allEdges.end(),
                 std::back_inserter(aboveThreshold),
                 [&selector](const EdgeInfo &e) {
                   return e.dihedralAngle >= selector.angleThreshold();
                 });
    return aboveThreshold;
  }
  }

  return {};
}

auto
Fillet3D::createFilletWedgeWithPlanes(
    const EdgeInfo                       &edge,
    const FilletParameters               &params,
    const std::optional<Kernel::Plane_3> &startPlane,
    const std::optional<Kernel::Plane_3> &endPlane) const
    -> CGAL::Polyhedron_3<Kernel>
{
  Surface_mesh_3 mesh;

  const double radius   = params.radius;
  const int    segments = std::max(4, params.segments);
  const auto   edgeDir  = normalizeVector(edge.end - edge.start);

  const auto into1 = -edge.normal1;
  const auto into2 = -edge.normal2;

  double dotInto = CGAL::to_double(into1 * into2);
  dotInto        = std::clamp(dotInto, -1.0, 1.0);
  const double arcAngle = M_PI - std::acos(-dotInto);

  if (arcAngle < 0.01 || arcAngle > M_PI - 0.01) {
#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] arcAngle out of range: " << (arcAngle * 180.0 / M_PI) << "°\n";
#endif
    return {};
  }

  // Solve for arc center offset using tangency constraints
  // This ensures the fillet arc is properly tangent to both faces for ANY dihedral angle
  auto arcCenterOffset =
      solveArcCenterOffset(edge.normal1, edge.normal2, edgeDir, radius);
  if (!arcCenterOffset.has_value()) {
#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] Failed to solve arc center offset (degenerate case)\n";
#endif
    return {}; // Degenerate case - skip edge
  }

#ifdef DEBUG_FILLET
  std::cerr << "[DEBUG] Creating wedge: arcAngle=" << (arcAngle * 180.0 / M_PI)
            << "°, offset=(" << CGAL::to_double(arcCenterOffset->x()) << ", "
            << CGAL::to_double(arcCenterOffset->y()) << ", "
            << CGAL::to_double(arcCenterOffset->z()) << ")\n";
#endif

  const auto baseArcCenterStart = edge.start + *arcCenterOffset;
  const auto baseArcCenterEnd   = edge.end + *arcCenterOffset;
  const double defaultExtension = radius * 0.5;

  const auto startRadial = normalizeVector(edge.normal1);

  const auto endRadial   = normalizeVector(edge.normal2);

  auto rotationAxis =
      normalizeVector(CGAL::cross_product(startRadial, endRadial));
  if (CGAL::to_double(rotationAxis.squared_length()) < EPSILON * EPSILON) {
    rotationAxis = edgeDir;
  }

  std::vector<Surface_mesh_3::Vertex_index> startArcVerts;
  std::vector<Surface_mesh_3::Vertex_index> endArcVerts;

  // Pre-compute extended arc points and use segment direction for intersection
  // This follows the "flat buffer" principle: intersect segment lines with planes
  for (int i = 0; i <= segments; ++i) {
    const double t     = static_cast<double>(i) / segments;
    const double angle = t * arcAngle;
    const auto   radial =
        normalizeVector(rotateVector(startRadial, rotationAxis, angle));

    // Base arc points on the original edge
    auto baseStartPt = baseArcCenterStart + radial * radius;
    auto baseEndPt   = baseArcCenterEnd + radial * radius;

    // Extended points (default when no bisector plane)
    auto extendedStartPt = baseStartPt - edgeDir * defaultExtension;
    auto extendedEndPt   = baseEndPt + edgeDir * defaultExtension;

    // Segment direction for intersection (like flat buffer)
    auto segmentDir = normalizeVector(extendedEndPt - extendedStartPt);

    Kernel::Point_3 startPt = extendedStartPt;
    Kernel::Point_3 endPt   = extendedEndPt;

    if (startPlane.has_value()) {
      // Intersect the segment line with the start plane
      startPt = intersectPointWithPlane(extendedStartPt, segmentDir, *startPlane);
    }

    if (endPlane.has_value()) {
      // Intersect the segment line with the end plane
      endPt = intersectPointWithPlane(extendedStartPt, segmentDir, *endPlane);
    }

    startArcVerts.push_back(mesh.add_vertex(startPt));
    endArcVerts.push_back(mesh.add_vertex(endPt));
  }

  // Position the wedge apex slightly outside the edge.
  // Use a small outward extension to ensure proper wedge volume.
  // The bisector direction (average of into vectors) points into the solid,
  // so we move opposite to it to place the apex outside.
  const auto   bisector   = normalizeVector(into1 + into2);
  const double outwardExt = radius * 0.05;

  // Base apex points
  auto baseApexStart = edge.start - bisector * outwardExt;
  auto baseApexEnd   = edge.end - bisector * outwardExt;

  // Extended apex points
  auto extendedApexStart = baseApexStart - edgeDir * defaultExtension;
  auto extendedApexEnd   = baseApexEnd + edgeDir * defaultExtension;

  // Apex segment direction
  auto apexSegmentDir = normalizeVector(extendedApexEnd - extendedApexStart);

  Kernel::Point_3 startEdgePt = extendedApexStart;
  Kernel::Point_3 endEdgePt   = extendedApexEnd;

  if (startPlane.has_value()) {
    startEdgePt =
        intersectPointWithPlane(extendedApexStart, apexSegmentDir, *startPlane);
  }

  if (endPlane.has_value()) {
    endEdgePt =
        intersectPointWithPlane(extendedApexStart, apexSegmentDir, *endPlane);
  }

  const auto startEdgeVert = mesh.add_vertex(startEdgePt);
  const auto endEdgeVert   = mesh.add_vertex(endEdgePt);

  // Curved fillet surface
  for (int i = 0; i < segments; ++i) {
    if (edge.isConvex) {
      mesh.add_face(startArcVerts[i], startArcVerts[i + 1], endArcVerts[i + 1]);
      mesh.add_face(startArcVerts[i], endArcVerts[i + 1], endArcVerts[i]);
    } else {
      mesh.add_face(startArcVerts[i], endArcVerts[i], endArcVerts[i + 1]);
      mesh.add_face(startArcVerts[i], endArcVerts[i + 1], startArcVerts[i + 1]);
    }
  }

  // Start cap
  for (int i = 0; i < segments; ++i) {
    if (edge.isConvex) {
      mesh.add_face(startEdgeVert, startArcVerts[i + 1], startArcVerts[i]);
    } else {
      mesh.add_face(startEdgeVert, startArcVerts[i], startArcVerts[i + 1]);
    }
  }

  // End cap
  for (int i = 0; i < segments; ++i) {
    if (edge.isConvex) {
      mesh.add_face(endEdgeVert, endArcVerts[i], endArcVerts[i + 1]);
    } else {
      mesh.add_face(endEdgeVert, endArcVerts[i + 1], endArcVerts[i]);
    }
  }

  // Side faces
  if (edge.isConvex) {
    mesh.add_face(startEdgeVert, startArcVerts[0], endArcVerts[0]);
    mesh.add_face(startEdgeVert, endArcVerts[0], endEdgeVert);
    mesh.add_face(startEdgeVert, endEdgeVert, endArcVerts[segments]);
    mesh.add_face(startEdgeVert, endArcVerts[segments], startArcVerts[segments]);
  } else {
    mesh.add_face(startEdgeVert, endArcVerts[0], startArcVerts[0]);
    mesh.add_face(startEdgeVert, endEdgeVert, endArcVerts[0]);
    mesh.add_face(startEdgeVert, endArcVerts[segments], endEdgeVert);
    mesh.add_face(startEdgeVert, startArcVerts[segments], endArcVerts[segments]);
  }

  CGAL::Polyhedron_3<Kernel> poly;
  CGAL::copy_face_graph(mesh, poly);

  if (poly.is_empty() || !poly.is_closed()) {
#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] Wedge poly invalid: empty=" << poly.is_empty() << ", closed=" << poly.is_closed() << "\n";
#endif
    return {};
  }

#ifdef DEBUG_FILLET
  std::cerr << "[DEBUG] Wedge created successfully with " << poly.size_of_vertices() << " vertices, " << poly.size_of_facets() << " facets\n";
  // Debug: output wedge vertices
  std::cerr << "[DEBUG] Wedge apex: (" << CGAL::to_double(startEdgePt.x()) << ", " << CGAL::to_double(startEdgePt.y()) << ", " << CGAL::to_double(startEdgePt.z()) << ")\n";
  std::cerr << "[DEBUG] Arc center: (" << CGAL::to_double(baseArcCenterStart.x()) << ", " << CGAL::to_double(baseArcCenterStart.y()) << ", " << CGAL::to_double(baseArcCenterStart.z()) << ")\n";
  std::cerr << "[DEBUG] bisector: (" << CGAL::to_double(bisector.x()) << ", " << CGAL::to_double(bisector.y()) << ", " << CGAL::to_double(bisector.z()) << ")\n";
  std::cerr << "[DEBUG] First arc point: (" << CGAL::to_double((baseArcCenterStart + startRadial * radius).x()) << ", " << CGAL::to_double((baseArcCenterStart + startRadial * radius).y()) << ", " << CGAL::to_double((baseArcCenterStart + startRadial * radius).z()) << ")\n";
#endif
  return poly;
}

auto
Fillet3D::toNefPolyhedron() const -> Nef_polyhedron
{
  try {
    if (_geometry->is<Solid>()) {
      const auto &solid = _geometry->as<Solid>();
      auto        poly_ptr =
          solid.exteriorShell().toPolyhedron_3<CGAL::Polyhedron_3<Kernel>>();
      if (poly_ptr && !poly_ptr->is_empty()) {
        return Nef_polyhedron(*poly_ptr);
      }
    } else if (_geometry->is<PolyhedralSurface>()) {
      const auto &ps       = _geometry->as<PolyhedralSurface>();
      auto        poly_ptr = ps.toPolyhedron_3<CGAL::Polyhedron_3<Kernel>>();
      if (poly_ptr && !poly_ptr->is_empty()) {
        return Nef_polyhedron(*poly_ptr);
      }
    }
  } catch (const std::exception &) {
    // NEF conversion failed (e.g., for some triangle prisms)
    return {};
  }

  return {};
}

auto
Fillet3D::fromNefPolyhedron(const Nef_polyhedron &nef) const
    -> std::unique_ptr<Geometry>
{
  if (nef.is_empty()) {
    return _geometry->clone();
  }

  // Helper lambda to build result from a polyhedron
  auto buildResult = [this](const CGAL::Polyhedron_3<Kernel> &poly)
      -> std::unique_ptr<Geometry> {
    if (poly.is_empty()) {
      return _geometry->clone();
    }

    auto result = std::make_unique<PolyhedralSurface>();
#ifdef DEBUG_FILLET
    int skippedFaces = 0;
    int addedFaces = 0;
#endif

    for (auto fit = poly.facets_begin(); fit != poly.facets_end(); ++fit) {
      std::vector<Kernel::Point_3> points;
      auto                         hit = fit->facet_begin();
      do {
        points.push_back(hit->vertex()->point());
        ++hit;
      } while (hit != fit->facet_begin());

      // Filter degenerate faces
      std::vector<Kernel::Point_3> uniquePoints;
      for (const auto &pt : points) {
        bool isDuplicate = false;
        for (const auto &existing : uniquePoints) {
          if (CGAL::squared_distance(pt, existing) < EPSILON * EPSILON) {
            isDuplicate = true;
            break;
          }
        }
        if (!isDuplicate) {
          uniquePoints.push_back(pt);
        }
      }

      if (uniquePoints.size() < 3) {
#ifdef DEBUG_FILLET
        ++skippedFaces;
#endif
        continue;
      }

      if (uniquePoints.size() == 3) {
        auto v1    = uniquePoints[1] - uniquePoints[0];
        auto v2    = uniquePoints[2] - uniquePoints[0];
        auto cross = CGAL::cross_product(v1, v2);
        if (CGAL::to_double(cross.squared_length()) < EPSILON * EPSILON) {
#ifdef DEBUG_FILLET
          ++skippedFaces;
#endif
          continue;
        }
      }

      LineString ring;
      for (const auto &pt : points) {
        ring.addPoint(Point(pt));
      }
      ring.addPoint(ring.pointN(0));
      result->addPolygon(Polygon(ring));
#ifdef DEBUG_FILLET
      ++addedFaces;
#endif
    }

#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] buildResult: added " << addedFaces << " faces, skipped " << skippedFaces << " degenerate faces\n";
#endif

    if (result->numPolygons() == 0) {
      return _geometry->clone();
    }

    if (_inputWasSolid) {
      auto solid             = std::make_unique<Solid>();
      solid->exteriorShell() = *result;
      return solid;
    }

    return result;
  };

  // First, try direct conversion (faster)
  try {
    CGAL::Polyhedron_3<Kernel> poly;
    nef.convert_to_polyhedron(poly);
#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] Direct NEF conversion: " << poly.size_of_vertices() << " vertices, " << poly.size_of_facets() << " facets\n";
#endif
    if (!poly.is_empty()) {
      return buildResult(poly);
    }
#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] Direct NEF conversion produced empty polyhedron, trying regularization\n";
#endif
    // Try regularizing the NEF first
    Nef_polyhedron regularized = nef.regularization();
    regularized.convert_to_polyhedron(poly);
#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] Regularized NEF conversion: " << poly.size_of_vertices() << " vertices, " << poly.size_of_facets() << " facets\n";
#endif
    if (!poly.is_empty()) {
      return buildResult(poly);
    }
  } catch (const std::exception &e) {
#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] Direct NEF conversion failed: " << e.what() << "\n";
#endif
    // Direct conversion failed, try polygon soup approach
  } catch (...) {
#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] Direct NEF conversion failed with unknown exception\n";
#endif
  }

  // Fallback: Extract faces via polygon soup and rebuild
  // This handles cases where direct conversion fails due to topology issues
  try {
    using Point_3 = Kernel::Point_3;
    std::vector<Point_3>             soup_points;
    std::vector<std::vector<size_t>> soup_polygons;
    std::map<Point_3, size_t>        point_index;

    // Helper to add a point and get its index
    auto addPoint = [&](const Point_3 &pt) -> size_t {
      auto it = point_index.find(pt);
      if (it == point_index.end()) {
        size_t idx      = soup_points.size();
        point_index[pt] = idx;
        soup_points.push_back(pt);
        return idx;
      }
      return it->second;
    };

    // Iterate over volumes to get proper boundary faces
    // For each bounded volume (mark() == true), extract its outer shell
    typename Nef_polyhedron::Volume_const_iterator vol;
    CGAL_forall_volumes(vol, nef) {
      // Only process marked volumes (solid parts)
      if (!vol->mark()) {
        continue;
      }

      // Get shells of this volume
      typename Nef_polyhedron::Shell_entry_const_iterator shell_it;
      CGAL_forall_shells_of(shell_it, vol) {
        // Get the halffacets of this shell
        typename Nef_polyhedron::SFace_const_handle sf(shell_it);
        typename Nef_polyhedron::SHalfedge_const_handle she = sf->sface_cycles_begin();

        // Visit all faces reachable from this shell
        // We'll use a workaround: iterate all halffacets and check if they belong to this volume
      }
    }

    // Alternative approach: iterate all halffacets and check incident volume
    typename Nef_polyhedron::Halffacet_const_iterator hf;
    CGAL_forall_halffacets(hf, nef) {
      // Check if this halffacet is on the boundary of a solid volume
      // The incident volume on the positive side should be marked (solid)
      // and on the negative side should be unmarked (void)
      auto vol_pos = hf->incident_volume();
      auto vol_neg = hf->twin()->incident_volume();

      // We want faces where one side is solid and other is void
      bool pos_solid = vol_pos->mark();
      bool neg_solid = vol_neg->mark();

      // Only process if this is a proper boundary (solid on one side, void on other)
      if (pos_solid == neg_solid) {
        continue;
      }

      // Use the halffacet facing outward from the solid (solid on positive side)
      const typename Nef_polyhedron::Halffacet_const_handle *face_to_use = nullptr;
      if (pos_solid && !neg_solid) {
        face_to_use = &hf;
      } else if (!pos_solid && neg_solid) {
        // Skip - we'll process this face via its twin
        continue;
      }

      if (!face_to_use) {
        continue;
      }

      // Extract polygons from this halffacet
      typename Nef_polyhedron::Halffacet_cycle_const_iterator fc;
      CGAL_forall_facet_cycles_of(fc, (*face_to_use)) {
        if (fc.is_shalfedge()) {
          typename Nef_polyhedron::SHalfedge_const_handle she_start(fc);
          typename Nef_polyhedron::SHalfedge_const_handle she = she_start;
          std::vector<size_t>                             polygon;

          do {
            Point_3 pt = she->source()->source()->point();
            polygon.push_back(addPoint(pt));
            she = she->next();
          } while (she != she_start);

          if (polygon.size() >= 3) {
            soup_polygons.push_back(polygon);
          }
        }
      }
    }

    if (soup_polygons.empty()) {
#ifdef DEBUG_FILLET
      std::cerr << "[DEBUG] Polygon soup: no polygons extracted\n";
#endif
      return _geometry->clone();
    }

#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] Polygon soup before repair: " << soup_points.size() << " points, " << soup_polygons.size() << " polygons\n";
#endif

    // Repair and orient the polygon soup
    namespace PMP = CGAL::Polygon_mesh_processing;
    PMP::repair_polygon_soup(soup_points, soup_polygons);
    PMP::orient_polygon_soup(soup_points, soup_polygons);

#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] Polygon soup after repair: " << soup_points.size() << " points, " << soup_polygons.size() << " polygons\n";
#endif

    // Convert to Surface_mesh
    Surface_mesh_3 mesh;
    PMP::polygon_soup_to_polygon_mesh(soup_points, soup_polygons, mesh);

    if (mesh.is_empty()) {
#ifdef DEBUG_FILLET
      std::cerr << "[DEBUG] Polygon soup: mesh conversion failed\n";
#endif
      return _geometry->clone();
    }

#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] Polygon soup: mesh has " << mesh.number_of_vertices() << " vertices, " << mesh.number_of_faces() << " faces\n";
    std::cerr << "[DEBUG] Mesh is_closed: " << mesh.is_valid() << ", closed: " << CGAL::is_closed(mesh) << "\n";
#endif

    // If the mesh isn't closed, try to repair it
    if (!CGAL::is_closed(mesh)) {
#ifdef DEBUG_FILLET
      std::cerr << "[DEBUG] Mesh not closed, trying to repair\n";
#endif
      PMP::stitch_borders(mesh);
#ifdef DEBUG_FILLET
      std::cerr << "[DEBUG] After stitch_borders: closed=" << CGAL::is_closed(mesh) << "\n";
#endif
    }

    // Triangulate faces to help with validity
    PMP::triangulate_faces(mesh);
#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] After triangulation: " << mesh.number_of_vertices() << " vertices, " << mesh.number_of_faces() << " faces\n";
#endif

    // Convert Surface_mesh to Polyhedron
    CGAL::Polyhedron_3<Kernel> poly;
    CGAL::copy_face_graph(mesh, poly);

#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] Final polyhedron: " << poly.size_of_vertices() << " vertices, " << poly.size_of_facets() << " facets, closed=" << poly.is_closed() << "\n";
#endif

    return buildResult(poly);

  } catch (const std::exception &) {
    return _geometry->clone();
  }
}

namespace {
// Check if geometry has triangular faces (which can cause NEF issues)
auto
hasTriangularFaces(const Geometry *geom) -> bool
{
  if (geom->is<Solid>()) {
    const auto &solid = geom->as<Solid>();
    for (size_t i = 0; i < solid.numShells(); ++i) {
      const auto &shell = solid.shellN(i);
      for (size_t j = 0; j < shell.numPolygons(); ++j) {
        if (shell.polygonN(j).exteriorRing().numPoints() <= 4) {
          return true; // Triangle (3 points + closing point = 4)
        }
      }
    }
  } else if (geom->is<PolyhedralSurface>()) {
    const auto &surface = geom->as<PolyhedralSurface>();
    for (size_t i = 0; i < surface.numPolygons(); ++i) {
      if (surface.polygonN(i).exteriorRing().numPoints() <= 4) {
        return true;
      }
    }
  }
  return false;
}
} // namespace

auto
Fillet3D::filletEdges(const EdgeSelector     &selector,
                      const FilletParameters &params) const
    -> std::unique_ptr<Geometry>
{
  // Triangle-based geometries can crash CGAL NEF operations
  if (hasTriangularFaces(_geometry)) {
    return _geometry->clone();
  }

  auto edges = selectEdges(selector);

  if (edges.empty()) {
    return _geometry->clone();
  }

  Nef_polyhedron inputNef = toNefPolyhedron();

  if (inputNef.is_empty()) {
    return _geometry->clone();
  }

  using PointKey = std::tuple<double, double, double>;
  std::map<PointKey, std::vector<size_t>> vertexToEdges;

  auto pointKey = [](const Kernel::Point_3 &p) -> PointKey {
    return {CGAL::to_double(p.x()), CGAL::to_double(p.y()),
            CGAL::to_double(p.z())};
  };

  for (size_t i = 0; i < edges.size(); ++i) {
    vertexToEdges[pointKey(edges[i].start)].push_back(i);
    vertexToEdges[pointKey(edges[i].end)].push_back(i);
  }

#ifdef DEBUG_FILLET
  std::cerr << "[DEBUG] Total edges found: " << edges.size() << "\n";
  for (size_t i = 0; i < edges.size(); ++i) {
    std::cerr << "[DEBUG] Edge " << i << ": ("
              << CGAL::to_double(edges[i].start.x()) << ", "
              << CGAL::to_double(edges[i].start.y()) << ", "
              << CGAL::to_double(edges[i].start.z()) << ") to ("
              << CGAL::to_double(edges[i].end.x()) << ", "
              << CGAL::to_double(edges[i].end.y()) << ", "
              << CGAL::to_double(edges[i].end.z()) << ") dihedral="
              << edges[i].dihedralAngle << "° convex=" << edges[i].isConvex << "\n";
  }
#endif

  // Detect and skip edges at reflex corners to avoid self-intersections
  // Only skip edges that are individually non-convex (concave/reflex)
  std::set<size_t> edgesToSkip;

  for (size_t i = 0; i < edges.size(); ++i) {
    if (!edges[i].isConvex) {
      edgesToSkip.insert(i);
#ifdef DEBUG_FILLET
      std::cerr << "[DEBUG] Skipping reflex edge " << i << ": ("
                << CGAL::to_double(edges[i].start.x()) << ", "
                << CGAL::to_double(edges[i].start.y()) << ", "
                << CGAL::to_double(edges[i].start.z()) << ") to ("
                << CGAL::to_double(edges[i].end.x()) << ", "
                << CGAL::to_double(edges[i].end.y()) << ", "
                << CGAL::to_double(edges[i].end.z()) << ") dihedral="
                << edges[i].dihedralAngle << "°\n";
#endif
    }
  }

  // For inner corner detection, we need ALL edges of the solid, not just selected edges.
  // This ensures proper detection even when using EXPLICIT edge selection.
  auto allEdges = findAllEdgesWithNormals();
  std::map<PointKey, std::vector<size_t>> allVertexToEdges;
  for (size_t i = 0; i < allEdges.size(); ++i) {
    allVertexToEdges[pointKey(allEdges[i].start)].push_back(i);
    allVertexToEdges[pointKey(allEdges[i].end)].push_back(i);
  }

  // Detect inner corner vertices by checking face normal directions relative to centroid.
  // At an outer corner, face normals point away from solid centroid.
  // At an inner corner, some face normals point toward solid centroid.
  std::set<PointKey> innerCornerVertices;
  std::map<PointKey, std::vector<Kernel::Vector_3>> innerCornerFaceNormals;

  // Compute solid centroid using ALL edge endpoints (not just selected)
  Kernel::Vector_3 centroidSum(0, 0, 0);
  int              pointCount = 0;
  for (const auto &e : allEdges) {
    centroidSum = centroidSum + (e.start - CGAL::ORIGIN);
    centroidSum = centroidSum + (e.end - CGAL::ORIGIN);
    pointCount += 2;
  }
  Kernel::Point_3 centroid =
      CGAL::ORIGIN + centroidSum / static_cast<double>(pointCount);

#ifdef DEBUG_FILLET
  std::cerr << "[DEBUG] Solid centroid: (" << CGAL::to_double(centroid.x()) << ", "
            << CGAL::to_double(centroid.y()) << ", " << CGAL::to_double(centroid.z())
            << ")\n";
#endif

  // Use ALL edges of solid for inner corner detection
  for (const auto &[key, edgeIndices] : allVertexToEdges) {
    if (edgeIndices.size() < 3) {
      continue;
    }

    const auto &[vx, vy, vz] = key;
    Kernel::Point_3 vertex(vx, vy, vz);

    // Vector from vertex to centroid
    Kernel::Vector_3 toCentroid = centroid - vertex;
    double           toCentroidLen =
        std::sqrt(CGAL::to_double(toCentroid.squared_length()));
    if (toCentroidLen < EPSILON) {
      continue;
    }
    toCentroid = toCentroid / toCentroidLen;

    // Collect unique face normals at this vertex from ALL edges
    std::vector<Kernel::Vector_3> faceNormals;
    for (size_t idx : edgeIndices) {
      const auto &e = allEdges[idx];
      bool        hasN1 = false, hasN2 = false;
      for (const auto &existing : faceNormals) {
        if (CGAL::to_double((existing - e.normal1).squared_length()) < EPSILON) {
          hasN1 = true;
        }
        if (CGAL::to_double((existing - e.normal2).squared_length()) < EPSILON) {
          hasN2 = true;
        }
      }
      if (!hasN1) {
        faceNormals.push_back(e.normal1);
      }
      if (!hasN2) {
        faceNormals.push_back(e.normal2);
      }
    }

    // At an outer corner, all face normals point away from centroid (dot < 0).
    // At an inner corner, at least one face normal points toward centroid (dot > 0).
    int normalsTowardCentroid = 0;
    for (const auto &n : faceNormals) {
      double dot = CGAL::to_double(n * toCentroid);
      if (dot > 0.1) { // threshold for "pointing toward"
        ++normalsTowardCentroid;
      }
    }

    // Inner corner: majority of normals point toward centroid
    if (normalsTowardCentroid >= 2) {
      innerCornerVertices.insert(key);
      innerCornerFaceNormals[key] = faceNormals;  // Store all face normals at this inner corner
#ifdef DEBUG_FILLET
      std::cerr << "[DEBUG] Detected inner corner at (" << vx << ", " << vy << ", "
                << vz << ") normalsTowardCentroid=" << normalsTowardCentroid << "\n";
#endif
    }
  }

  // At inner corner vertices, we want to keep only the edge that forms the
  // actual inner corner (typically vertical) and skip edges that connect to
  // top/bottom faces. The inner corner edge is the one whose normals are both
  // pointing toward the centroid.
  for (size_t i = 0; i < edges.size(); ++i) {
    if (edgesToSkip.find(i) != edgesToSkip.end()) {
      continue;
    }
    auto startKey = pointKey(edges[i].start);
    auto endKey   = pointKey(edges[i].end);

    // Check if either endpoint is at an inner corner
    bool startIsInner = innerCornerVertices.count(startKey) > 0;
    bool endIsInner   = innerCornerVertices.count(endKey) > 0;

    if (startIsInner || endIsInner) {
      // Check if this edge IS the inner corner edge (both normals point toward centroid)
      const auto &e = edges[i];
      double      dot1 = 0, dot2 = 0;

      if (startIsInner) {
        Kernel::Point_3  startPt(std::get<0>(startKey), std::get<1>(startKey),
                                 std::get<2>(startKey));
        Kernel::Vector_3 toCenter = centroid - startPt;
        toCenter = normalizeVector(toCenter);
        dot1     = CGAL::to_double(e.normal1 * toCenter);
        dot2     = CGAL::to_double(e.normal2 * toCenter);
      } else {
        Kernel::Point_3  endPt(std::get<0>(endKey), std::get<1>(endKey),
                               std::get<2>(endKey));
        Kernel::Vector_3 toCenter = centroid - endPt;
        toCenter = normalizeVector(toCenter);
        dot1     = CGAL::to_double(e.normal1 * toCenter);
        dot2     = CGAL::to_double(e.normal2 * toCenter);
      }

      // If BOTH normals point toward centroid, this is the inner corner edge - SKIP it
      // We'll handle it specially by subtracting a cube and adding a sphere at the corners
      bool isInnerCornerEdge = (dot1 > 0.1 && dot2 > 0.1);

      if (isInnerCornerEdge) {
        edgesToSkip.insert(i);
#ifdef DEBUG_FILLET
        std::cerr << "[DEBUG] Skipping inner corner edge " << i << " (dots=" << dot1
                  << "," << dot2 << ") - will handle with cube+sphere\n";
#endif
      }
    }
  }

  // Filter edges
  std::vector<EdgeInfo> filteredEdges;
  for (size_t i = 0; i < edges.size(); ++i) {
    if (edgesToSkip.find(i) == edgesToSkip.end()) {
      filteredEdges.push_back(edges[i]);
    }
  }

  if (filteredEdges.empty()) {
#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] No edges left after filtering\n";
#endif
    return _geometry->clone();
  }
#ifdef DEBUG_FILLET
  std::cerr << "[DEBUG] Processing " << filteredEdges.size() << " edges\n";
#endif

  // Rebuild vertex mapping
  vertexToEdges.clear();
  for (size_t i = 0; i < filteredEdges.size(); ++i) {
    vertexToEdges[pointKey(filteredEdges[i].start)].push_back(i);
    vertexToEdges[pointKey(filteredEdges[i].end)].push_back(i);
  }

  edges = std::move(filteredEdges);

  // Compute bisector planes at junctions
  std::map<std::pair<PointKey, size_t>, Kernel::Plane_3> edgeBisectorPlanes;

  for (const auto &[key, edgeIndices] : vertexToEdges) {
    if (edgeIndices.size() < 2) {
      continue;
    }

    const auto     &[vx, vy, vz] = key;
    Kernel::Point_3 vertex(vx, vy, vz);

    std::vector<std::pair<size_t, Kernel::Vector_3>> edgeDirs;
    for (size_t idx : edgeIndices) {
      const auto &e   = edges[idx];
      auto        dir = (CGAL::squared_distance(e.start, vertex) < EPSILON)
                            ? normalizeVector(e.end - e.start)
                            : normalizeVector(e.start - e.end);
      edgeDirs.emplace_back(idx, dir);
    }

    if (edgeIndices.size() == 2) {
      // Two edges: bisector is sum of both directions (like buffer3D)
      auto bisectorNormal =
          normalizeVector(edgeDirs[0].second + edgeDirs[1].second);
      auto plane                                = Kernel::Plane_3(vertex, bisectorNormal);
      edgeBisectorPlanes[{key, edgeDirs[0].first}] = plane;
      edgeBisectorPlanes[{key, edgeDirs[1].first}] = plane;
    } else {
      // 3+ edges: for each edge, find the adjacent edge (closest angle) and compute bisector
      for (size_t i = 0; i < edgeDirs.size(); ++i) {
        const auto &[edgeIdx, dir] = edgeDirs[i];

        // Find the edge with smallest angle to this one (most adjacent)
        double     bestDot = -2.0;
        size_t     bestJ   = (i + 1) % edgeDirs.size();
        for (size_t j = 0; j < edgeDirs.size(); ++j) {
          if (j == i)
            continue;
          double dot = CGAL::to_double(dir * edgeDirs[j].second);
          if (dot > bestDot) {
            bestDot = dot;
            bestJ   = j;
          }
        }

        // Compute bisector with the most adjacent edge
        auto bisectorNormal = normalizeVector(dir + edgeDirs[bestJ].second);
        auto plane          = Kernel::Plane_3(vertex, bisectorNormal);
        edgeBisectorPlanes[{key, edgeIdx}] = plane;
      }
    }
  }

  // Create wedges and subtract them one at a time
  // This avoids topology issues when wedges intersect at corners
  Nef_polyhedron resultNef = inputNef;
  int            wedgesProcessed = 0;

  for (size_t i = 0; i < edges.size(); ++i) {
    const auto &edge = edges[i];

    try {
      std::optional<Kernel::Plane_3> startPlane;
      std::optional<Kernel::Plane_3> endPlane;

      auto startKey = pointKey(edge.start);
      auto endKey   = pointKey(edge.end);

      // At inner corner vertices, use a perpendicular face plane (not a 45° bisector).
      // Find a face normal at the inner corner that is NOT part of this edge.
      // This ensures each edge's wedge is cut by the perpendicular face, producing valid geometry.
      if (innerCornerVertices.count(startKey) > 0) {
        Kernel::Point_3 startPt(std::get<0>(startKey), std::get<1>(startKey),
                                std::get<2>(startKey));
        Kernel::Vector_3 toCenter = centroid - startPt;
        toCenter = normalizeVector(toCenter);

        // Find a face normal at this vertex that is different from edge.normal1 and edge.normal2
        const auto &allNormals = innerCornerFaceNormals[startKey];
        Kernel::Vector_3 planeNormal = edge.normal1;  // default
        double bestScore = -2.0;
        for (const auto &n : allNormals) {
          // Check if this normal is different from both edge normals
          bool isDifferent = (CGAL::to_double((n - edge.normal1).squared_length()) > 0.01 &&
                              CGAL::to_double((n - edge.normal2).squared_length()) > 0.01);
          if (!isDifferent) continue;

          // Prefer the normal pointing most toward centroid (the "inner" face)
          double dotCenter = CGAL::to_double(n * toCenter);
          if (dotCenter > bestScore) {
            bestScore = dotCenter;
            planeNormal = n;
          }
        }
        startPlane = Kernel::Plane_3(startPt, planeNormal);
#ifdef DEBUG_FILLET
        std::cerr << "[DEBUG] Inner corner start plane: normal=("
                  << CGAL::to_double(planeNormal.x()) << ", "
                  << CGAL::to_double(planeNormal.y()) << ", "
                  << CGAL::to_double(planeNormal.z()) << ") bestScore=" << bestScore << "\n";
#endif
      } else {
        if (auto it = edgeBisectorPlanes.find({startKey, i});
            it != edgeBisectorPlanes.end()) {
          startPlane = it->second;
        }
      }

      if (innerCornerVertices.count(endKey) > 0) {
        Kernel::Point_3 endPt(std::get<0>(endKey), std::get<1>(endKey),
                              std::get<2>(endKey));
        Kernel::Vector_3 toCenter = centroid - endPt;
        toCenter = normalizeVector(toCenter);

        // Find a face normal at this vertex that is different from edge.normal1 and edge.normal2
        const auto &allNormals = innerCornerFaceNormals[endKey];
        Kernel::Vector_3 planeNormal = edge.normal1;  // default
        double bestScore = -2.0;
        for (const auto &n : allNormals) {
          bool isDifferent = (CGAL::to_double((n - edge.normal1).squared_length()) > 0.01 &&
                              CGAL::to_double((n - edge.normal2).squared_length()) > 0.01);
          if (!isDifferent) continue;

          double dotCenter = CGAL::to_double(n * toCenter);
          if (dotCenter > bestScore) {
            bestScore = dotCenter;
            planeNormal = n;
          }
        }
        endPlane = Kernel::Plane_3(endPt, planeNormal);
#ifdef DEBUG_FILLET
        std::cerr << "[DEBUG] Inner corner end plane: normal=("
                  << CGAL::to_double(planeNormal.x()) << ", "
                  << CGAL::to_double(planeNormal.y()) << ", "
                  << CGAL::to_double(planeNormal.z()) << ") bestScore=" << bestScore << "\n";
#endif
      } else {
        if (auto it = edgeBisectorPlanes.find({endKey, i});
            it != edgeBisectorPlanes.end()) {
          endPlane = it->second;
        }
      }

      auto wedge = createFilletWedgeWithPlanes(edge, params, startPlane, endPlane);

      if (wedge.empty() || !wedge.is_closed()) {
        continue;
      }

      Nef_polyhedron wedgeNef(wedge);

      if (wedgeNef.is_empty()) {
        continue;
      }

      // Subtract this wedge immediately
      // This avoids the topology issues that arise from unioning wedges first
      resultNef = resultNef - wedgeNef;
      ++wedgesProcessed;
#ifdef DEBUG_FILLET
      std::cerr << "[DEBUG] Subtracted wedge " << (i+1) << "/" << edges.size() << "\n";
#endif

    } catch (const std::exception &e) {
#ifdef DEBUG_FILLET
      std::cerr << "[DEBUG] Exception processing edge " << i << ": " << e.what() << "\n";
#endif
      continue;
    }
  }

  // === INNER CORNER VERTEX TREATMENT ===
  // Disabled: The edge fillets naturally meet at inner corners.
  // Adding an extra corner wedge creates an unwanted concave bubble.
  // The inner corner edge should be filleted like any other edge.
  for (const auto &cornerKey : innerCornerVertices) {
    (void)cornerKey;  // Suppress unused variable warning
    continue;  // Skip inner corner treatment - let edge fillets handle it naturally
    double vx = std::get<0>(cornerKey);
    double vy = std::get<1>(cornerKey);
    double vz = std::get<2>(cornerKey);

    // Get face normals at this inner corner
    auto it = innerCornerFaceNormals.find(cornerKey);
    if (it == innerCornerFaceNormals.end() || it->second.size() < 2) {
      continue;
    }

    // Find the two horizontal face normals (perpendicular to Z)
    Kernel::Vector_3 n1(0, 0, 0), n2(0, 0, 0);
    Kernel::Vector_3 nz(0, 0, 0);  // vertical normal (for top/bottom face)
    for (const auto &n : it->second) {
      double nzComp = std::abs(CGAL::to_double(n.z()));
      if (nzComp < 0.1) {
        // Horizontal normal
        if (CGAL::to_double(n1 * n1) < EPSILON) {
          n1 = n;
        } else {
          n2 = n;
        }
      } else {
        nz = n;
      }
    }

    if (CGAL::to_double(n1 * n1) < EPSILON || CGAL::to_double(n2 * n2) < EPSILON) {
      continue;  // Couldn't find two horizontal normals
    }

#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] Inner corner vertex at (" << vx << ", " << vy << ", " << vz << ")\n"
              << "  n1=(" << CGAL::to_double(n1.x()) << ", " << CGAL::to_double(n1.y()) << ")\n"
              << "  n2=(" << CGAL::to_double(n2.x()) << ", " << CGAL::to_double(n2.y()) << ")\n";
#endif

    try {
      // Create corner joint as a curved triangular prism
      // The goal is to connect the two edge fillet endpoints with a smooth arc.
      //
      // Geometry:
      // - Apex at corner vertex (vx, vy)
      // - Arc endpoint A at (vx, vy - r) = where X-edge fillet ends
      // - Arc endpoint B at (vx - r, vy) = where Y-edge fillet ends
      // - Arc from A to B centered at (vx - r, vy - r) with radius r
      //
      // The shape is a triangular prism with one curved edge (the arc).

      // Only process bottom vertices (nz.z() < 0) to avoid duplicate processing
      if (CGAL::to_double(nz.z()) >= 0) {
        continue;
      }

      double cornerHeight = 2.0;  // Height of inner corner edge
      double zBottom = vz;
      double zTop = vz + cornerHeight;

      // n1 and n2 are face normals pointing OUTWARD from the inner corner
      // For L-shape at (1,1): n1=(0,1), n2=(1,0)
      // Arc center is offset INTO the solid by radius in direction opposite to normals
      double arcCenterX = vx - params.radius * CGAL::to_double(n2.x());  // 1 - 0.15 = 0.85
      double arcCenterY = vy - params.radius * CGAL::to_double(n1.y());  // 1 - 0.15 = 0.85

      // Small offset to avoid coincidence with edge wedge end faces
      double eps = params.radius * 0.001;

      // Arc endpoints (where the edge fillets end) - slightly offset inward
      double arcPtAx = vx - eps;                                    // ~0.99985
      double arcPtAy = vy - params.radius * CGAL::to_double(n1.y());  // 0.85
      double arcPtBx = vx - params.radius * CGAL::to_double(n2.x());  // 0.85
      double arcPtBy = vy - eps;                                    // ~0.99985

      // Apex at the arc center (not corner) to properly connect with edge fillets
      double apexX = arcCenterX;
      double apexY = arcCenterY;

#ifdef DEBUG_FILLET
      std::cerr << "[DEBUG] Corner joint: apex=(" << vx << ", " << vy << ")"
                << " arcCenter=(" << arcCenterX << ", " << arcCenterY << ")"
                << " A=(" << arcPtAx << ", " << arcPtAy << ")"
                << " B=(" << arcPtBx << ", " << arcPtBy << ")"
                << " z=[" << zBottom << ", " << zTop << "]\n";
#endif

      // Build the curved triangular prism
      unsigned int numSegs = std::max(8u, static_cast<unsigned int>(params.segments));
      CGAL::Polyhedron_3<Kernel> wedgePoly;

      struct CurvedTriangleBuilder : public CGAL::Modifier_base<CGAL::Polyhedron_3<Kernel>::HalfedgeDS> {
        double apexX, apexY;      // Apex (corner vertex)
        double arcCx, arcCy, r;   // Arc center and radius
        double z0, z1;            // Z extents
        double startAngle, endAngle;  // Arc angles
        unsigned int segs;

        CurvedTriangleBuilder(double ax, double ay, double cx, double cy, double r_,
                              double z0_, double z1_, double sa, double ea, unsigned int s)
            : apexX(ax), apexY(ay), arcCx(cx), arcCy(cy), r(r_),
              z0(z0_), z1(z1_), startAngle(sa), endAngle(ea), segs(s) {}

        void operator()(CGAL::Polyhedron_3<Kernel>::HalfedgeDS &hds) {
          CGAL::Polyhedron_incremental_builder_3<CGAL::Polyhedron_3<Kernel>::HalfedgeDS> B(hds, true);

          // Vertices:
          // - Arc points at z0: indices 0 to segs
          // - Arc points at z1: indices (segs+1) to (2*segs+1)
          // - Apex at z0: index 2*(segs+1)
          // - Apex at z1: index 2*(segs+1)+1
          unsigned int numArcPts = segs + 1;
          unsigned int apex0 = 2 * numArcPts;
          unsigned int apex1 = 2 * numArcPts + 1;
          unsigned int numVerts = 2 * numArcPts + 2;

          // Faces:
          // - Arc surface: segs quads
          // - Bottom cap: triangle fan from apex0 (segs triangles)
          // - Top cap: triangle fan from apex1 (segs triangles)
          // - End face 1: quad from apex to arc[0]
          // - End face 2: quad from apex to arc[segs]
          unsigned int numFaces = segs + segs + segs + 2;

          B.begin_surface(numVerts, numFaces, 0);

          // Add arc vertices at z0
          for (unsigned int i = 0; i <= segs; ++i) {
            double t = static_cast<double>(i) / segs;
            double angle = startAngle + t * (endAngle - startAngle);
            double px = arcCx + r * std::cos(angle);
            double py = arcCy + r * std::sin(angle);
            B.add_vertex(Kernel::Point_3(px, py, z0));
          }
          // Add arc vertices at z1
          for (unsigned int i = 0; i <= segs; ++i) {
            double t = static_cast<double>(i) / segs;
            double angle = startAngle + t * (endAngle - startAngle);
            double px = arcCx + r * std::cos(angle);
            double py = arcCy + r * std::sin(angle);
            B.add_vertex(Kernel::Point_3(px, py, z1));
          }
          // Add apex vertices
          B.add_vertex(Kernel::Point_3(apexX, apexY, z0));
          B.add_vertex(Kernel::Point_3(apexX, apexY, z1));

          // Arc surface quads - normals point OUTWARD from the wedge (toward arc center)
          for (unsigned int i = 0; i < segs; ++i) {
            B.begin_facet();
            B.add_vertex_to_facet(i + 1);                // z0 arc[i+1]
            B.add_vertex_to_facet(i);                    // z0 arc[i]
            B.add_vertex_to_facet(numArcPts + i);        // z1 arc[i]
            B.add_vertex_to_facet(numArcPts + i + 1);    // z1 arc[i+1]
            B.end_facet();
          }

          // Bottom cap: triangle fan from apex0
          for (unsigned int i = 0; i < segs; ++i) {
            B.begin_facet();
            B.add_vertex_to_facet(apex0);
            B.add_vertex_to_facet(i);
            B.add_vertex_to_facet(i + 1);
            B.end_facet();
          }

          // Top cap: triangle fan from apex1 (reversed winding)
          for (unsigned int i = 0; i < segs; ++i) {
            B.begin_facet();
            B.add_vertex_to_facet(apex1);
            B.add_vertex_to_facet(numArcPts + i + 1);
            B.add_vertex_to_facet(numArcPts + i);
            B.end_facet();
          }

          // End face 1: quad from apex to arc[0] (connects to Y-edge fillet)
          B.begin_facet();
          B.add_vertex_to_facet(apex0);
          B.add_vertex_to_facet(apex1);
          B.add_vertex_to_facet(numArcPts);       // z1 arc[0]
          B.add_vertex_to_facet(0);               // z0 arc[0]
          B.end_facet();

          // End face 2: quad from apex to arc[segs] (connects to X-edge fillet)
          B.begin_facet();
          B.add_vertex_to_facet(segs);              // z0 arc[segs]
          B.add_vertex_to_facet(numArcPts + segs);  // z1 arc[segs]
          B.add_vertex_to_facet(apex1);
          B.add_vertex_to_facet(apex0);
          B.end_facet();

          B.end_surface();
        }
      };

      // Arc angles: from n2 direction to n1 direction
      // n2 = (1, 0) -> angle = 0°
      // n1 = (0, 1) -> angle = 90°
      // Arc goes from 0° to 90° (counterclockwise)
      double startAng = std::atan2(CGAL::to_double(n2.y()), CGAL::to_double(n2.x()));  // 0°
      double endAng = std::atan2(CGAL::to_double(n1.y()), CGAL::to_double(n1.x()));    // 90°

      // Handle angle wrap
      if (endAng < startAng) {
        endAng += 2 * M_PI;
      }

#ifdef DEBUG_FILLET
      std::cerr << "[DEBUG] Arc angles: " << (startAng * 180 / M_PI) << "° to "
                << (endAng * 180 / M_PI) << "°\n";
#endif

      CurvedTriangleBuilder builder(apexX, apexY, arcCenterX, arcCenterY, params.radius,
                                    zBottom, zTop, startAng, endAng, numSegs);
      wedgePoly.delegate(builder);

      if (!wedgePoly.is_closed() || wedgePoly.empty()) {
#ifdef DEBUG_FILLET
        std::cerr << "[DEBUG] Corner wedge not closed, skipping\n";
#endif
        continue;
      }

#ifdef DEBUG_FILLET
      std::cerr << "[DEBUG] Corner wedge: " << wedgePoly.size_of_vertices() << " vertices, "
                << wedgePoly.size_of_facets() << " facets, closed=" << wedgePoly.is_closed() << "\n";
#endif

      Nef_polyhedron wedgeNef(wedgePoly);
      if (!wedgeNef.is_empty()) {
        resultNef = resultNef - wedgeNef;
#ifdef DEBUG_FILLET
        std::cerr << "[DEBUG] Subtracted corner wedge at inner corner\n";
#endif
      }

    } catch (const std::exception &ex) {
#ifdef DEBUG_FILLET
      std::cerr << "[DEBUG] Exception at inner corner: " << ex.what() << "\n";
#endif
    }
  }

  // Note: Sphere octant treatment at inner corner vertices was attempted but caused
  // invalid geometry due to overlap with corner wedge. The corner wedge alone is
  // sufficient for handling inner corner edges.

  if (wedgesProcessed == 0) {
#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] No wedges processed\n";
#endif
    return _geometry->clone();
  }

#ifdef DEBUG_FILLET
  std::cerr << "[DEBUG] Processed " << wedgesProcessed << " wedges, converting result\n";
  std::cerr << "[DEBUG] Result NEF: is_empty=" << resultNef.is_empty() << ", number_of_vertices=" << resultNef.number_of_vertices() << ", number_of_facets=" << resultNef.number_of_facets() << ", number_of_volumes=" << resultNef.number_of_volumes() << "\n";
  // Count marked volumes
  int markedVolumes = 0;
  typename Fillet3D::Nef_polyhedron::Volume_const_iterator vol_it;
  CGAL_forall_volumes(vol_it, resultNef) {
    if (vol_it->mark()) ++markedVolumes;
  }
  std::cerr << "[DEBUG] Marked volumes: " << markedVolumes << "\n";
#endif

  if (resultNef.is_empty()) {
    return _geometry->clone();
  }

  return fromNefPolyhedron(resultNef);
}

auto
Fillet3D::fillet(double radius, int segments) const -> std::unique_ptr<Geometry>
{
  return filletEdges(EdgeSelector::allConvex(),
                     FilletParameters::constant(radius, segments));
}

auto
fillet3D(const Geometry &geometry, double radius, int segments)
    -> std::unique_ptr<Geometry>
{
  Fillet3D fillet(geometry);
  return fillet.fillet(radius, segments);
}

} // namespace SFCGAL::algorithm
