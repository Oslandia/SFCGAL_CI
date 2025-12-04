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
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <cmath>
#include <iostream>
#include <map>
#include <set>

// Uncomment to enable debug output
// #define DEBUG_FILLET  // Disabled for production

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
  const auto cross = CGAL::cross_product(n1, n2);
  return CGAL::to_double(cross * edgeDir) > 0;
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
    if (normals.second.has_value()) {
      EdgeInfo info;
      info.start   = edgePoints.at(key).first;
      info.end     = edgePoints.at(key).second;
      info.normal1 = normals.first;
      info.normal2 = *normals.second;

      const auto edgeDir = normalizeVector(info.end - info.start);
      info.dihedralAngle = computeDihedralAngle(info.normal1, info.normal2);
      info.isConvex      = isEdgeConvex(info.normal1, info.normal2, edgeDir);

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

  for (int i = 0; i <= segments; ++i) {
    const double t     = static_cast<double>(i) / segments;
    const double angle = t * arcAngle;
    const auto radial = normalizeVector(rotateVector(startRadial, rotationAxis, angle));

    auto startPt = baseArcCenterStart + radial * radius;
    auto endPt   = baseArcCenterEnd + radial * radius;

    if (startPlane.has_value()) {
      startPt = intersectPointWithPlane(startPt, -edgeDir, *startPlane);
    } else {
      startPt = startPt - edgeDir * defaultExtension;
    }

    if (endPlane.has_value()) {
      endPt = intersectPointWithPlane(endPt, edgeDir, *endPlane);
    } else {
      endPt = endPt + edgeDir * defaultExtension;
    }

    startArcVerts.push_back(mesh.add_vertex(startPt));
    endArcVerts.push_back(mesh.add_vertex(endPt));
  }

  // Position the wedge apex slightly outside the edge.
  // Use a small outward extension to ensure proper wedge volume.
  // The bisector direction (average of into vectors) points into the solid,
  // so we move opposite to it to place the apex outside.
  const auto   bisector    = normalizeVector(into1 + into2);
  const double outwardExt  = radius * 0.05;
  auto         startEdgePt = edge.start - bisector * outwardExt;
  auto         endEdgePt   = edge.end - bisector * outwardExt;

  if (startPlane.has_value()) {
    startEdgePt = intersectPointWithPlane(startEdgePt, -edgeDir, *startPlane);
  } else {
    startEdgePt = startEdgePt - edgeDir * defaultExtension;
  }

  if (endPlane.has_value()) {
    endEdgePt = intersectPointWithPlane(endEdgePt, edgeDir, *endPlane);
  } else {
    endEdgePt = endEdgePt + edgeDir * defaultExtension;
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

  try {
    CGAL::Polyhedron_3<Kernel> poly;
    nef.convert_to_polyhedron(poly);

    if (poly.is_empty()) {
      return _geometry->clone();
    }

    auto result = std::make_unique<PolyhedralSurface>();

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
        continue;
      }

      if (uniquePoints.size() == 3) {
        auto v1    = uniquePoints[1] - uniquePoints[0];
        auto v2    = uniquePoints[2] - uniquePoints[0];
        auto cross = CGAL::cross_product(v1, v2);
        if (CGAL::to_double(cross.squared_length()) < EPSILON * EPSILON) {
          continue;
        }
      }

      LineString ring;
      for (const auto &pt : points) {
        ring.addPoint(Point(pt));
      }
      ring.addPoint(ring.pointN(0));
      result->addPolygon(Polygon(ring));
    }

    if (result->numPolygons() == 0) {
      return _geometry->clone();
    }

    if (_inputWasSolid) {
      auto solid             = std::make_unique<Solid>();
      solid->exteriorShell() = *result;
      return solid;
    }

    return result;
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

  // Detect and skip edges at reflex corners to avoid self-intersections
  // Only skip edges that are individually non-convex (concave/reflex)
  std::set<size_t> edgesToSkip;

  for (size_t i = 0; i < edges.size(); ++i) {
    if (!edges[i].isConvex) {
      edgesToSkip.insert(i);
    }
  }

  // Also skip edges connecting two reflex vertices
  std::set<PointKey> problematicVertices;
  for (size_t idx : edgesToSkip) {
    problematicVertices.insert(pointKey(edges[idx].start));
    problematicVertices.insert(pointKey(edges[idx].end));
  }

  for (size_t i = 0; i < edges.size(); ++i) {
    if (edgesToSkip.find(i) != edgesToSkip.end()) {
      continue;
    }
    auto startKey = pointKey(edges[i].start);
    auto endKey   = pointKey(edges[i].end);
    if (problematicVertices.count(startKey) > 0 &&
        problematicVertices.count(endKey) > 0) {
      edgesToSkip.insert(i);
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
      auto bisectorNormal =
          normalizeVector(edgeDirs[0].second + edgeDirs[1].second);
      auto plane                                = Kernel::Plane_3(vertex, bisectorNormal);
      edgeBisectorPlanes[{key, edgeDirs[0].first}] = plane;
      edgeBisectorPlanes[{key, edgeDirs[1].first}] = plane;
    } else {
      for (size_t i = 0; i < edgeDirs.size(); ++i) {
        const auto &[edgeIdx, dir]             = edgeDirs[i];
        auto plane                             = Kernel::Plane_3(vertex, dir);
        edgeBisectorPlanes[{key, edgeIdx}] = plane;
      }
    }
  }

  // Create wedges
  Nef_polyhedron allWedges;
  bool           firstWedge = true;

  for (size_t i = 0; i < edges.size(); ++i) {
    const auto &edge = edges[i];

    try {
      std::optional<Kernel::Plane_3> startPlane;
      std::optional<Kernel::Plane_3> endPlane;

      auto startKey = pointKey(edge.start);
      auto endKey   = pointKey(edge.end);

      if (auto it = edgeBisectorPlanes.find({startKey, i});
          it != edgeBisectorPlanes.end()) {
        startPlane = it->second;
      }
      if (auto it = edgeBisectorPlanes.find({endKey, i});
          it != edgeBisectorPlanes.end()) {
        endPlane = it->second;
      }

      auto wedge = createFilletWedgeWithPlanes(edge, params, startPlane, endPlane);

      if (wedge.empty() || !wedge.is_closed()) {
        continue;
      }

      Nef_polyhedron wedgeNef(wedge);

      if (wedgeNef.is_empty()) {
        continue;
      }

      if (firstWedge) {
        allWedges  = std::move(wedgeNef);
        firstWedge = false;
      } else {
        allWedges = allWedges + wedgeNef;
      }

    } catch (const std::exception &) {
      continue;
    }
  }

  if (firstWedge || allWedges.is_empty()) {
#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] No wedges created\n";
#endif
    return _geometry->clone();
  }

#ifdef DEBUG_FILLET
  std::cerr << "[DEBUG] Performing NEF subtraction\n";
#endif
  try {
    auto resultNef = inputNef - allWedges;
#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] NEF subtraction completed, result empty=" << resultNef.is_empty() << "\n";
#endif

    if (resultNef.is_empty()) {
      return _geometry->clone();
    }

    return fromNefPolyhedron(resultNef);
  } catch (const std::exception &e) {
#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] NEF subtraction exception: " << e.what() << "\n";
#endif
    // Boolean operation failed
    return _geometry->clone();
  } catch (...) {
#ifdef DEBUG_FILLET
    std::cerr << "[DEBUG] NEF subtraction unknown exception\n";
#endif
    return _geometry->clone();
  }
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
