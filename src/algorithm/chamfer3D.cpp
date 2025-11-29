// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/chamfer3D.h"
#include "SFCGAL/algorithm/difference.h"
#include "SFCGAL/algorithm/union.h"
#include "SFCGAL/io/wkt.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/detail/TypeForDimension.h"
#include "SFCGAL/numeric.h"

#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/merge_border_vertices.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <map>
#include <set>
#include <sstream>

namespace PMP = CGAL::Polygon_mesh_processing;

namespace SFCGAL::algorithm {

// ============================================================================
// Helper functions
// ============================================================================

namespace {

/**
 * @brief Compute the dihedral angle between two faces sharing an edge
 * @param n1 Normal of first face
 * @param n2 Normal of second face
 * @return Dihedral angle in radians [0, PI]
 */
auto
computeDihedralAngleFromNormals(const Kernel::Vector_3 &n1,
                                const Kernel::Vector_3 &n2) -> double
{
  // Normalize vectors
  double len1 = std::sqrt(CGAL::to_double(n1.squared_length()));
  double len2 = std::sqrt(CGAL::to_double(n2.squared_length()));

  if (len1 < EPSILON || len2 < EPSILON) {
    return 0.0;
  }

  Kernel::Vector_3 n1_norm = n1 / len1;
  Kernel::Vector_3 n2_norm = n2 / len2;

  // Dot product gives cos(angle between normals)
  // Dihedral angle = PI - angle between normals (for outward normals)
  double dot =
      CGAL::to_double(n1_norm.x() * n2_norm.x() + n1_norm.y() * n2_norm.y() +
                      n1_norm.z() * n2_norm.z());

  // Clamp to [-1, 1] to avoid numerical issues with acos
  dot = std::max(-1.0, std::min(1.0, dot));

  // The dihedral angle is PI minus the angle between normals
  return M_PI - std::acos(dot);
}

/**
 * @brief Check if two points are approximately equal
 */
auto
pointsApproxEqual(const Kernel::Point_3 &p1, const Kernel::Point_3 &p2,
                  double tolerance = EPSILON) -> bool
{
  return CGAL::to_double(CGAL::squared_distance(p1, p2)) <
         tolerance * tolerance;
}

/**
 * @brief Normalize a vector
 */
auto
normalize(const Kernel::Vector_3 &v) -> Kernel::Vector_3
{
  double len = std::sqrt(CGAL::to_double(v.squared_length()));
  if (len < EPSILON) {
    return Kernel::Vector_3(0, 0, 0);
  }
  return v / len;
}

/**
 * @brief Compute bisector plane for two edge directions meeting at a vertex
 *
 * At an inner corner where two edges meet, the bisector plane is used to create
 * a miter joint. The plane passes through the corner vertex and its normal is
 * the bisector of the two edge directions.
 *
 * @param cornerVertex The vertex where edges meet
 * @param edgeDir1 Direction of first edge (pointing away from corner)
 * @param edgeDir2 Direction of second edge (pointing away from corner)
 * @return Bisector plane
 */
auto
computeBisectorPlane(const Kernel::Point_3 &cornerVertex,
                     const Kernel::Vector_3 &edgeDir1,
                     const Kernel::Vector_3 &edgeDir2) -> Kernel::Plane_3
{
  // Bisector direction is the average of the two edge directions
  Kernel::Vector_3 bisector = normalize(edgeDir1 + edgeDir2);
  return Kernel::Plane_3(cornerVertex, bisector);
}

/**
 * @brief Intersect a line segment with a plane
 *
 * Returns the intersection point of the line (p1, p2) with the plane.
 * The point may be beyond the segment endpoints.
 *
 * @param p1 First point of segment
 * @param p2 Second point of segment
 * @param plane The plane to intersect with
 * @return Intersection point
 */
auto
intersectSegmentPlane(const Kernel::Point_3 &p1,
                      const Kernel::Point_3 &p2,
                      const Kernel::Plane_3 &plane) -> Kernel::Point_3
{
  Kernel::Vector_3 v = p2 - p1;

  // Plane equation: ax + by + cz + d = 0
  // Line: p = p1 + t * v
  // Substitute: a(p1.x + t*v.x) + b(p1.y + t*v.y) + c(p1.z + t*v.z) + d = 0
  // Solve for t: t = -(a*p1.x + b*p1.y + c*p1.z + d) / (a*v.x + b*v.y + c*v.z)

  double denominator = CGAL::to_double(plane.a() * v.x() +
                                        plane.b() * v.y() +
                                        plane.c() * v.z());

  if (std::abs(denominator) < EPSILON) {
    // Line is parallel to plane, return midpoint as fallback
    return Kernel::Point_3(
        (CGAL::to_double(p1.x()) + CGAL::to_double(p2.x())) / 2.0,
        (CGAL::to_double(p1.y()) + CGAL::to_double(p2.y())) / 2.0,
        (CGAL::to_double(p1.z()) + CGAL::to_double(p2.z())) / 2.0);
  }

  double numerator = -CGAL::to_double(plane.a() * p1.x() +
                                       plane.b() * p1.y() +
                                       plane.c() * p1.z() +
                                       plane.d());
  double t = numerator / denominator;

  return Kernel::Point_3(
      CGAL::to_double(p1.x()) + t * CGAL::to_double(v.x()),
      CGAL::to_double(p1.y()) + t * CGAL::to_double(v.y()),
      CGAL::to_double(p1.z()) + t * CGAL::to_double(v.z()));
}

/**
 * @brief Check if a triangle is degenerate (has repeated or nearly-coincident vertices)
 */
auto
isTriangleDegenerate(const Kernel::Point_3 &p0, const Kernel::Point_3 &p1,
                     const Kernel::Point_3 &p2, double tolerance = EPSILON) -> bool
{
  // Check for repeated vertices
  if (pointsApproxEqual(p0, p1, tolerance) ||
      pointsApproxEqual(p1, p2, tolerance) ||
      pointsApproxEqual(p0, p2, tolerance)) {
    return true;
  }
  // Check for collinear vertices (zero-area triangle)
  Kernel::Vector_3 v1 = p1 - p0;
  Kernel::Vector_3 v2 = p2 - p0;
  Kernel::Vector_3 cross = CGAL::cross_product(v1, v2);
  double area2 = CGAL::to_double(cross.squared_length());
  return area2 < tolerance * tolerance;
}

/**
 * @brief Clean a Solid by removing degenerate faces
 */
auto
cleanDegenerateFaces(const Solid &solid) -> std::unique_ptr<Solid>
{
  // Get the exterior shell and create a new one without degenerate faces
  const PolyhedralSurface &shell = solid.exteriorShell();
  auto cleanShell = std::make_unique<PolyhedralSurface>();

  for (size_t i = 0; i < shell.numPolygons(); ++i) {
    const Polygon &poly = shell.polygonN(i);
    const LineString &ring = poly.exteriorRing();

    // For triangulated surfaces, each polygon should have exactly 4 points (3 + closing)
    if (ring.numPoints() >= 4) {
      const Point &pt0 = ring.pointN(0);
      const Point &pt1 = ring.pointN(1);
      const Point &pt2 = ring.pointN(2);

      Kernel::Point_3 p0(pt0.x(), pt0.y(), pt0.z());
      Kernel::Point_3 p1(pt1.x(), pt1.y(), pt1.z());
      Kernel::Point_3 p2(pt2.x(), pt2.y(), pt2.z());

      if (!isTriangleDegenerate(p0, p1, p2)) {
        cleanShell->addPatch(poly);
      }
    }
  }

  return std::make_unique<Solid>(*cleanShell);
}

/**
 * @brief Detect if an edge is concave (reflex) using the perpendicular cross-product test
 *
 * For a truly concave edge (reflex, > 180° exterior angle):
 * - (perp1 × perp2) · edgeDir > 0
 *
 * For a convex edge (< 180° exterior angle):
 * - (perp1 × perp2) · edgeDir < 0
 *
 * The perpendiculars (pointing INTO the solid) give a more reliable indicator
 * of the solid's geometry at the edge than the face normals.
 *
 * @param perp1 Perpendicular vector for first face (points into solid)
 * @param perp2 Perpendicular vector for second face (points into solid)
 * @param edgeDir Edge direction vector (from user-specified edge start to end)
 * @return true if the edge is truly concave (reflex)
 */
auto
isConcaveEdge(const Kernel::Vector_3 &perp1,
              const Kernel::Vector_3 &perp2,
              const Kernel::Vector_3 &edgeDir) -> bool
{
  // Cross product of the two perpendiculars (which point into the solid)
  Kernel::Vector_3 cross = CGAL::cross_product(perp1, perp2);

  // Dot with edge direction determines concavity
  double dot = CGAL::to_double(cross * edgeDir);

  // Positive dot means concave (reflex edge)
  return dot > 0.01;  // Small threshold to handle numerical precision
}

/**
 * @brief Remove degenerate faces from a mesh
 *
 * A face is degenerate if all its vertices are at the same point,
 * or if it has zero area (e.g., collinear vertices in a triangle).
 */
void
removeDegenerateFaces(Surface_mesh_3 &mesh)
{
  std::vector<Surface_mesh_3::Face_index> toRemove;

  for (auto face : mesh.faces()) {
    // Collect unique vertex positions for this face
    std::vector<Kernel::Point_3> points;
    for (auto v : vertices_around_face(mesh.halfedge(face), mesh)) {
      points.push_back(mesh.point(v));
    }

    if (points.size() < 3) {
      toRemove.push_back(face);
      continue;
    }

    // Check if all points are (nearly) the same
    bool allSame = true;
    for (size_t i = 1; i < points.size() && allSame; ++i) {
      if (!pointsApproxEqual(points[0], points[i], EPSILON)) {
        allSame = false;
      }
    }
    if (allSame) {
      toRemove.push_back(face);
      continue;
    }

    // Check for collinear points (zero-area triangle)
    if (points.size() == 3) {
      Kernel::Vector_3 v1 = points[1] - points[0];
      Kernel::Vector_3 v2 = points[2] - points[0];
      Kernel::Vector_3 cross = CGAL::cross_product(v1, v2);
      double           area2 = CGAL::to_double(cross.squared_length());
      if (area2 < EPSILON * EPSILON) {
        toRemove.push_back(face);
      }
    }
  }

  // Remove degenerate faces
  for (auto face : toRemove) {
    CGAL::Euler::remove_face(mesh.halfedge(face), mesh);
  }
}

/**
 * @brief Create a triangular prism (wedge) for filling a concave edge chamfer
 *
 * For concave edges, instead of clipping material, we need to ADD material
 * by creating a triangular prism that fills the void and unioning it.
 *
 * @param edgeStart Start point of the edge
 * @param edgeEnd End point of the edge
 * @param edgeDir Normalized edge direction
 * @param perp1 Perpendicular direction in face 1 plane (pointing into void)
 * @param perp2 Perpendicular direction in face 2 plane (pointing into void)
 * @param d1 Chamfer distance on face 1
 * @param d2 Chamfer distance on face 2
 * @return A Solid representing the fill wedge, or nullptr on failure
 */
auto
createConcaveFillWedge(const Kernel::Point_3  &edgeStart,
                       const Kernel::Point_3  &edgeEnd,
                       const Kernel::Vector_3 &edgeDir,
                       const Kernel::Vector_3 &perp1,
                       const Kernel::Vector_3 &perp2, double d1, double d2)
    -> std::unique_ptr<Solid>
{
  // For a concave chamfer, we create a triangular prism that:
  // - Has its edge along the concave edge of the solid
  // - Extends into the void by d1 along perp1 and d2 along perp2
  // - The chamfer face is the hypotenuse of the triangular cross-section

  // Ensure start/end are ordered correctly relative to edgeDir
  Kernel::Point_3 orderedStart = edgeStart;
  Kernel::Point_3 orderedEnd = edgeEnd;
  Kernel::Vector_3 actualDir = edgeEnd - edgeStart;
  if (CGAL::to_double(actualDir * edgeDir) < 0) {
    // Edge points are reversed relative to edgeDir, swap them
    std::swap(orderedStart, orderedEnd);
  }

  // No extension beyond edge endpoints - keeps fill wedge within solid bounds.
  // The fill wedge is processed BEFORE convex edges, so no extension is needed.
  Kernel::Point_3 extStart = orderedStart;
  Kernel::Point_3 extEnd   = orderedEnd;

  // The wedge cross-section is a right triangle:
  // - Vertex 0: on the edge (slightly inside solid to ensure overlap for union)
  // - Vertex 1: offset by d1 along perp1 (on face 1 cut line)
  // - Vertex 2: offset by d2 along perp2 (on face 2 cut line)

  // Inset toward solid for clean union - must be large enough to ensure overlap
  double          inset = std::min(d1, d2) * 0.3;
  Kernel::Vector_3 insetDir = normalize(perp1 + perp2);

  // Vertices at start of edge
  Kernel::Point_3 s0 = extStart - insetDir * inset; // Slightly into solid
  Kernel::Point_3 s1 = extStart + perp1 * d1;       // On face 1 cut line
  Kernel::Point_3 s2 = extStart + perp2 * d2;       // On face 2 cut line

  // Vertices at end of edge
  Kernel::Point_3 e0 = extEnd - insetDir * inset;
  Kernel::Point_3 e1 = extEnd + perp1 * d1;
  Kernel::Point_3 e2 = extEnd + perp2 * d2;

  // Create the triangular prism as a PolyhedralSurface
  auto ps = std::make_unique<PolyhedralSurface>();

  // Front triangle (at start, viewed from outside looking at -edgeDir)
  {
    auto *ring = new LineString();
    ring->addPoint(Point(s0));
    ring->addPoint(Point(s2));
    ring->addPoint(Point(s1));
    ring->addPoint(Point(s0));
    ps->addPatch(Polygon(ring));
  }

  // Back triangle (at end, viewed from outside looking at +edgeDir)
  {
    auto *ring = new LineString();
    ring->addPoint(Point(e0));
    ring->addPoint(Point(e1));
    ring->addPoint(Point(e2));
    ring->addPoint(Point(e0));
    ps->addPatch(Polygon(ring));
  }

  // Side face 1 (s0-s1-e1-e0): adjacent to face 1 of solid
  {
    auto *ring = new LineString();
    ring->addPoint(Point(s0));
    ring->addPoint(Point(s1));
    ring->addPoint(Point(e1));
    ring->addPoint(Point(e0));
    ring->addPoint(Point(s0));
    ps->addPatch(Polygon(ring));
  }

  // Side face 2 (s0-e0-e2-s2): adjacent to face 2 of solid
  {
    auto *ring = new LineString();
    ring->addPoint(Point(s0));
    ring->addPoint(Point(e0));
    ring->addPoint(Point(e2));
    ring->addPoint(Point(s2));
    ring->addPoint(Point(s0));
    ps->addPatch(Polygon(ring));
  }

  // Chamfer face (s1-s2-e2-e1): the visible chamfer surface
  {
    auto *ring = new LineString();
    ring->addPoint(Point(s1));
    ring->addPoint(Point(s2));
    ring->addPoint(Point(e2));
    ring->addPoint(Point(e1));
    ring->addPoint(Point(s1));
    ps->addPatch(Polygon(ring));
  }

  return std::make_unique<Solid>(*ps);
}

/**
 * @brief Create a bounded subtraction wedge for convex edges using actual mesh endpoints
 *
 * This is used when infinite plane clipping would fail due to complex geometry.
 * The wedge is bounded to the actual edge extent without extension, to avoid
 * overlapping with adjacent edge chamfers at corners.
 *
 * @param edgeStart Actual start point on the mesh edge
 * @param edgeEnd Actual end point on the mesh edge
 * @param edgeDir Normalized edge direction
 * @param perp1 Perpendicular direction in face 1 plane (pointing into solid)
 * @param perp2 Perpendicular direction in face 2 plane (pointing into solid)
 * @param d1 Chamfer distance on face 1
 * @param d2 Chamfer distance on face 2
 * @return A Solid representing the subtraction wedge, or nullptr on failure
 */
auto
createBoundedSubtractionWedge(const Kernel::Point_3  &edgeStart,
                               const Kernel::Point_3  &edgeEnd,
                               const Kernel::Vector_3 &edgeDir,
                               const Kernel::Vector_3 &perp1,
                               const Kernel::Vector_3 &perp2,
                               double d1, double d2)
    -> std::unique_ptr<Solid>
{
  // The wedge cross-section is a right triangle that cuts the corner:
  // - Vertex 0: at the edge (slightly beyond to ensure clean cut)
  // - Vertex 1: offset by d1 along perp1 (on face 1 cut line)
  // - Vertex 2: offset by d2 along perp2 (on face 2 cut line)

  // Compute the chamfer plane normal
  Kernel::Vector_3 chamferNormal = normalize(perp1 * d2 + perp2 * d1);

  // Small extension beyond edge for robust boolean (but bounded to edge endpoints)
  double edgeExtension = std::min(d1, d2) * 0.1;

  // Vertices at start of edge (NO extension beyond actual endpoint)
  Kernel::Point_3 s0 = edgeStart - chamferNormal * edgeExtension;
  Kernel::Point_3 s1 = edgeStart + perp1 * d1;
  Kernel::Point_3 s2 = edgeStart + perp2 * d2;

  // Vertices at end of edge (NO extension beyond actual endpoint)
  Kernel::Point_3 e0 = edgeEnd - chamferNormal * edgeExtension;
  Kernel::Point_3 e1 = edgeEnd + perp1 * d1;
  Kernel::Point_3 e2 = edgeEnd + perp2 * d2;

  // Create triangular prism solid
  auto ps = std::make_unique<PolyhedralSurface>();

  // Front triangle face (at start)
  {
    auto *ring = new LineString();
    ring->addPoint(Point(s0));
    ring->addPoint(Point(s1));
    ring->addPoint(Point(s2));
    ring->addPoint(Point(s0));
    ps->addPatch(Polygon(ring));
  }

  // Back triangle face (at end)
  {
    auto *ring = new LineString();
    ring->addPoint(Point(e0));
    ring->addPoint(Point(e2));
    ring->addPoint(Point(e1));
    ring->addPoint(Point(e0));
    ps->addPatch(Polygon(ring));
  }

  // Face 1 side (s0-e0-e1-s1)
  {
    auto *ring = new LineString();
    ring->addPoint(Point(s0));
    ring->addPoint(Point(e0));
    ring->addPoint(Point(e1));
    ring->addPoint(Point(s1));
    ring->addPoint(Point(s0));
    ps->addPatch(Polygon(ring));
  }

  // Face 2 side (s0-s2-e2-e0)
  {
    auto *ring = new LineString();
    ring->addPoint(Point(s0));
    ring->addPoint(Point(s2));
    ring->addPoint(Point(e2));
    ring->addPoint(Point(e0));
    ring->addPoint(Point(s0));
    ps->addPatch(Polygon(ring));
  }

  // Chamfer face (s1-e1-e2-s2)
  {
    auto *ring = new LineString();
    ring->addPoint(Point(s1));
    ring->addPoint(Point(e1));
    ring->addPoint(Point(e2));
    ring->addPoint(Point(s2));
    ring->addPoint(Point(s1));
    ps->addPatch(Polygon(ring));
  }

  return std::make_unique<Solid>(*ps);
}

/**
 * @brief Create a subtraction wedge with optional bisector plane clipping for miter joints
 *
 * This function creates a chamfer wedge that can be clipped at one or both ends
 * by bisector planes, creating proper miter joints at inner corners.
 *
 * @param edgeStart Actual start point on the mesh edge
 * @param edgeEnd Actual end point on the mesh edge
 * @param edgeDir Normalized edge direction
 * @param perp1 Perpendicular direction in face 1 plane (pointing into solid)
 * @param perp2 Perpendicular direction in face 2 plane (pointing into solid)
 * @param d1 Chamfer distance on face 1
 * @param d2 Chamfer distance on face 2
 * @param startBisectorPlane Optional bisector plane for start end (nullptr if none)
 * @param endBisectorPlane Optional bisector plane for end (nullptr if none)
 * @return A Solid representing the subtraction wedge, or nullptr on failure
 */
auto
createMiteredSubtractionWedge(const Kernel::Point_3  &edgeStart,
                               const Kernel::Point_3  &edgeEnd,
                               const Kernel::Vector_3 &edgeDir,
                               const Kernel::Vector_3 &perp1,
                               const Kernel::Vector_3 &perp2,
                               double d1, double d2,
                               const Kernel::Plane_3 *startBisectorPlane,
                               const Kernel::Plane_3 *endBisectorPlane)
    -> std::unique_ptr<Solid>
{
  // Compute the chamfer plane normal
  Kernel::Vector_3 chamferNormal = normalize(perp1 * d2 + perp2 * d1);

  // Small extension beyond edge for robust boolean
  double edgeExtension = std::min(d1, d2) * 0.1;

  // Extended points for wedge construction (extend past corner for miter)
  double miterExtension = std::max(d1, d2) * 2.0; // Extend significantly for clean miter
  Kernel::Point_3 extStart = edgeStart;
  Kernel::Point_3 extEnd = edgeEnd;

  if (startBisectorPlane != nullptr) {
    extStart = edgeStart - edgeDir * miterExtension;
  }
  if (endBisectorPlane != nullptr) {
    extEnd = edgeEnd + edgeDir * miterExtension;
  }

  // Base vertices (will be clipped if bisector plane exists)
  Kernel::Point_3 s0_base = extStart - chamferNormal * edgeExtension;
  Kernel::Point_3 s1_base = extStart + perp1 * d1;
  Kernel::Point_3 s2_base = extStart + perp2 * d2;

  Kernel::Point_3 e0_base = extEnd - chamferNormal * edgeExtension;
  Kernel::Point_3 e1_base = extEnd + perp1 * d1;
  Kernel::Point_3 e2_base = extEnd + perp2 * d2;

  // Clip start vertices to bisector plane
  Kernel::Point_3 s0 = s0_base, s1 = s1_base, s2 = s2_base;
  if (startBisectorPlane != nullptr) {
    s0 = intersectSegmentPlane(s0_base, e0_base, *startBisectorPlane);
    s1 = intersectSegmentPlane(s1_base, e1_base, *startBisectorPlane);
    s2 = intersectSegmentPlane(s2_base, e2_base, *startBisectorPlane);
  }

  // Clip end vertices to bisector plane
  Kernel::Point_3 e0 = e0_base, e1 = e1_base, e2 = e2_base;
  if (endBisectorPlane != nullptr) {
    e0 = intersectSegmentPlane(s0_base, e0_base, *endBisectorPlane);
    e1 = intersectSegmentPlane(s1_base, e1_base, *endBisectorPlane);
    e2 = intersectSegmentPlane(s2_base, e2_base, *endBisectorPlane);
  }

  // Create triangular prism solid
  auto ps = std::make_unique<PolyhedralSurface>();

  // Front triangle face (at start)
  {
    auto *ring = new LineString();
    ring->addPoint(Point(s0));
    ring->addPoint(Point(s1));
    ring->addPoint(Point(s2));
    ring->addPoint(Point(s0));
    ps->addPatch(Polygon(ring));
  }

  // Back triangle face (at end)
  {
    auto *ring = new LineString();
    ring->addPoint(Point(e0));
    ring->addPoint(Point(e2));
    ring->addPoint(Point(e1));
    ring->addPoint(Point(e0));
    ps->addPatch(Polygon(ring));
  }

  // Face 1 side (s0-e0-e1-s1)
  {
    auto *ring = new LineString();
    ring->addPoint(Point(s0));
    ring->addPoint(Point(e0));
    ring->addPoint(Point(e1));
    ring->addPoint(Point(s1));
    ring->addPoint(Point(s0));
    ps->addPatch(Polygon(ring));
  }

  // Face 2 side (s0-s2-e2-e0)
  {
    auto *ring = new LineString();
    ring->addPoint(Point(s0));
    ring->addPoint(Point(s2));
    ring->addPoint(Point(e2));
    ring->addPoint(Point(e0));
    ring->addPoint(Point(s0));
    ps->addPatch(Polygon(ring));
  }

  // Chamfer face (s1-e1-e2-s2)
  {
    auto *ring = new LineString();
    ring->addPoint(Point(s1));
    ring->addPoint(Point(e1));
    ring->addPoint(Point(e2));
    ring->addPoint(Point(s2));
    ring->addPoint(Point(s1));
    ps->addPatch(Polygon(ring));
  }

  return std::make_unique<Solid>(*ps);
}

/**
 * @brief Create a corner clipping tetrahedron for vertices where 3 edges meet
 *
 * When chamfering multiple edges that meet at a corner vertex, the edge wedges
 * don't fully remove the corner vertex. This function creates a tetrahedron that
 * clips the remaining corner material.
 *
 * @param cornerVertex The corner vertex to clip
 * @param chamferPoints Vector of 3 chamfer points (one from each edge meeting at corner)
 * @return A Solid tetrahedron for subtraction, or nullptr on failure
 */
auto
createCornerClipTetrahedron(const Kernel::Point_3 &cornerVertex,
                             const std::vector<Kernel::Point_3> &chamferPoints)
    -> std::unique_ptr<Solid>
{
  if (chamferPoints.size() != 3) {
    return nullptr; // Only handle 3-edge corners for now
  }

  const auto &p1 = chamferPoints[0];
  const auto &p2 = chamferPoints[1];
  const auto &p3 = chamferPoints[2];

  // Compute centroid of chamfer points
  Kernel::Point_3 centroid(
      (CGAL::to_double(p1.x()) + CGAL::to_double(p2.x()) + CGAL::to_double(p3.x())) / 3.0,
      (CGAL::to_double(p1.y()) + CGAL::to_double(p2.y()) + CGAL::to_double(p3.y())) / 3.0,
      (CGAL::to_double(p1.z()) + CGAL::to_double(p2.z()) + CGAL::to_double(p3.z())) / 3.0);

  // Direction from centroid to corner
  Kernel::Vector_3 toCorner = cornerVertex - centroid;

  // Extend apex significantly beyond corner for clean cut
  // The extension needs to be large enough to fully encompass the corner material
  double cornerDist = std::sqrt(CGAL::to_double(toCorner.squared_length()));
  double extension = cornerDist * 0.5; // 50% extension for robust clip
  Kernel::Point_3 apex = cornerVertex + normalize(toCorner) * extension;

  // Create tetrahedron: 4 triangular faces
  auto ps = std::make_unique<PolyhedralSurface>();

  // Base triangle (p1-p2-p3) - check winding with respect to apex
  // Normal should point away from apex
  Kernel::Vector_3 baseNormal = CGAL::cross_product(p2 - p1, p3 - p1);
  Kernel::Vector_3 toApex = apex - p1;
  bool flipBase = CGAL::to_double(baseNormal * toApex) > 0;

  // Base face
  {
    auto *ring = new LineString();
    if (flipBase) {
      ring->addPoint(Point(p1));
      ring->addPoint(Point(p3));
      ring->addPoint(Point(p2));
    } else {
      ring->addPoint(Point(p1));
      ring->addPoint(Point(p2));
      ring->addPoint(Point(p3));
    }
    ring->addPoint(ring->pointN(0)); // Close the ring
    ps->addPatch(Polygon(ring));
  }

  // Side face 1: apex-p1-p2
  {
    auto *ring = new LineString();
    if (flipBase) {
      ring->addPoint(Point(apex));
      ring->addPoint(Point(p2));
      ring->addPoint(Point(p1));
    } else {
      ring->addPoint(Point(apex));
      ring->addPoint(Point(p1));
      ring->addPoint(Point(p2));
    }
    ring->addPoint(ring->pointN(0));
    ps->addPatch(Polygon(ring));
  }

  // Side face 2: apex-p2-p3
  {
    auto *ring = new LineString();
    if (flipBase) {
      ring->addPoint(Point(apex));
      ring->addPoint(Point(p3));
      ring->addPoint(Point(p2));
    } else {
      ring->addPoint(Point(apex));
      ring->addPoint(Point(p2));
      ring->addPoint(Point(p3));
    }
    ring->addPoint(ring->pointN(0));
    ps->addPatch(Polygon(ring));
  }

  // Side face 3: apex-p3-p1
  {
    auto *ring = new LineString();
    if (flipBase) {
      ring->addPoint(Point(apex));
      ring->addPoint(Point(p1));
      ring->addPoint(Point(p3));
    } else {
      ring->addPoint(Point(apex));
      ring->addPoint(Point(p3));
      ring->addPoint(Point(p1));
    }
    ring->addPoint(ring->pointN(0));
    ps->addPatch(Polygon(ring));
  }

  return std::make_unique<Solid>(*ps);
}

/**
 * @brief Create a corner clip tetrahedron for vertices where 2 convex edges meet
 *
 * When 2 convex edges meet at a vertex (with a concave/reflex edge also present),
 * the chamfer wedges from each edge leave a gap - a tetrahedron-shaped region of
 * material that needs to be removed.
 *
 * The geometry is determined by the perp1/perp2 vectors of each edge's chamfer:
 * - Both edges share one face, so they have the same perp2 direction
 * - Each edge has a different perp1 direction (into its own adjacent face)
 *
 * The tetrahedron has 4 vertices:
 * - V0: corner vertex
 * - V1: corner + perp1_edge1 * d (edge 1's perp1 offset)
 * - V2: corner + perp2_shared * d (shared perp2 offset - same for both edges)
 * - V3: corner + perp1_edge2 * d (edge 2's perp1 offset)
 *
 * @param cornerVertex The corner vertex
 * @param perp1_e1 Perpendicular direction into edge 1's adjacent face
 * @param perp2_shared Shared perpendicular direction (into the common face)
 * @param perp1_e2 Perpendicular direction into edge 2's adjacent face
 * @param d Chamfer distance
 * @return A Solid tetrahedron for subtraction, or nullptr on failure
 */
auto
createTwoEdgeCornerTetrahedron(const Kernel::Point_3 &cornerVertex,
                                const Kernel::Vector_3 &perp1_e1,
                                const Kernel::Vector_3 &perp2_shared,
                                const Kernel::Vector_3 &perp1_e2,
                                double d)
    -> std::unique_ptr<Solid>
{
  // Compute the 4 vertices of the tetrahedron
  // The base vertices are the chamfer face corners
  Kernel::Point_3 v1 = cornerVertex + perp1_e1 * d;    // Edge 1's perp1 offset
  Kernel::Point_3 v2 = cornerVertex + perp2_shared * d; // Shared perp2 offset
  Kernel::Point_3 v3 = cornerVertex + perp1_e2 * d;    // Edge 2's perp1 offset

  // Compute centroid of base triangle
  Kernel::Point_3 centroid(
      (CGAL::to_double(v1.x()) + CGAL::to_double(v2.x()) + CGAL::to_double(v3.x())) / 3.0,
      (CGAL::to_double(v1.y()) + CGAL::to_double(v2.y()) + CGAL::to_double(v3.y())) / 3.0,
      (CGAL::to_double(v1.z()) + CGAL::to_double(v2.z()) + CGAL::to_double(v3.z())) / 3.0);

  // Direction from centroid to corner (outward direction)
  Kernel::Vector_3 toCorner = cornerVertex - centroid;
  double cornerDist = std::sqrt(CGAL::to_double(toCorner.squared_length()));

  // Extend apex significantly beyond corner for clean cut (50% like 3-edge case)
  double ext = cornerDist * 0.5;
  Kernel::Point_3 v0 = cornerVertex + normalize(toCorner) * ext; // Apex extended beyond corner

  // Create tetrahedron: 4 triangular faces
  auto ps = std::make_unique<PolyhedralSurface>();

  // Determine correct winding by checking face normals point outward
  // Base triangle (v1-v2-v3) - check winding with respect to v0
  Kernel::Vector_3 baseNormal = CGAL::cross_product(v2 - v1, v3 - v1);
  Kernel::Vector_3 toApex = v0 - v1;
  bool flipBase = CGAL::to_double(baseNormal * toApex) > 0;

  // Base face: v1-v2-v3
  {
    auto *ring = new LineString();
    if (flipBase) {
      ring->addPoint(Point(v1));
      ring->addPoint(Point(v3));
      ring->addPoint(Point(v2));
    } else {
      ring->addPoint(Point(v1));
      ring->addPoint(Point(v2));
      ring->addPoint(Point(v3));
    }
    ring->addPoint(ring->pointN(0));
    ps->addPatch(Polygon(ring));
  }

  // Side face 1: v0-v1-v2
  {
    auto *ring = new LineString();
    if (flipBase) {
      ring->addPoint(Point(v0));
      ring->addPoint(Point(v2));
      ring->addPoint(Point(v1));
    } else {
      ring->addPoint(Point(v0));
      ring->addPoint(Point(v1));
      ring->addPoint(Point(v2));
    }
    ring->addPoint(ring->pointN(0));
    ps->addPatch(Polygon(ring));
  }

  // Side face 2: v0-v2-v3
  {
    auto *ring = new LineString();
    if (flipBase) {
      ring->addPoint(Point(v0));
      ring->addPoint(Point(v3));
      ring->addPoint(Point(v2));
    } else {
      ring->addPoint(Point(v0));
      ring->addPoint(Point(v2));
      ring->addPoint(Point(v3));
    }
    ring->addPoint(ring->pointN(0));
    ps->addPatch(Polygon(ring));
  }

  // Side face 3: v0-v3-v1
  {
    auto *ring = new LineString();
    if (flipBase) {
      ring->addPoint(Point(v0));
      ring->addPoint(Point(v1));
      ring->addPoint(Point(v3));
    } else {
      ring->addPoint(Point(v0));
      ring->addPoint(Point(v3));
      ring->addPoint(Point(v1));
    }
    ring->addPoint(ring->pointN(0));
    ps->addPatch(Polygon(ring));
  }

  return std::make_unique<Solid>(*ps);
}

} // anonymous namespace

// ============================================================================
// Chamfer3D Implementation
// ============================================================================

Chamfer3D::Chamfer3D(const Geometry &geometry) : _geometry(&geometry)
{
  if (geometry.is<Solid>()) {
    _isSolid = true;
  } else if (geometry.is<PolyhedralSurface>()) {
    _isPolyhedralSurface = true;
  } else {
    throw std::invalid_argument(
        "Chamfer3D: Input geometry must be Solid or PolyhedralSurface");
  }

  if (geometry.isEmpty()) {
    throw std::invalid_argument("Chamfer3D: Input geometry is empty");
  }
}

auto
Chamfer3D::toSurfaceMesh() const -> Surface_mesh_3
{
  if (_isSolid) {
    return _geometry->as<Solid>().toSurfaceMesh();
  }
  return _geometry->as<PolyhedralSurface>().toSurfaceMesh();
}

auto
Chamfer3D::fromSurfaceMesh(const Surface_mesh_3 &mesh) const
    -> std::unique_ptr<Geometry>
{
  auto ps = std::make_unique<PolyhedralSurface>(mesh);

  if (_isSolid) {
    return std::make_unique<Solid>(*ps);
  }
  return ps;
}

/**
 * @brief Extract edges from original geometry patches (not triangulated mesh)
 *
 * This avoids including internal diagonal edges that are artifacts of triangulation.
 * A cube has 12 real edges, but triangulation creates 6 additional diagonal edges.
 */
auto
Chamfer3D::extractOriginalEdges() const -> std::vector<EdgeIdentifier>
{
  std::set<EdgeIdentifier> edgeSet;

  const PolyhedralSurface* ps = nullptr;
  if (_isSolid) {
    ps = &_geometry->as<Solid>().exteriorShell();
  } else {
    ps = &_geometry->as<PolyhedralSurface>();
  }

  // Extract edges from each polygon patch
  for (size_t i = 0; i < ps->numPatches(); ++i) {
    const Polygon& patch = ps->patchN(i);
    const LineString& ring = patch.exteriorRing();

    for (size_t j = 0; j < ring.numPoints() - 1; ++j) {
      const Point& p1 = ring.pointN(j);
      const Point& p2 = ring.pointN(j + 1);
      edgeSet.emplace(p1.toPoint_3(), p2.toPoint_3());
    }
  }

  return std::vector<EdgeIdentifier>(edgeSet.begin(), edgeSet.end());
}

auto
Chamfer3D::resolveEdges(const EdgeSelector &selector) const
    -> std::vector<EdgeIdentifier>
{
  if (selector.mode == EdgeSelector::Mode::EXPLICIT) {
    return selector.edges;
  }

  // BY_ANGLE and ALL modes are not supported - they return multiple edges
  // that can cause corner overlap issues with boolean operations
  throw std::invalid_argument(
      "Chamfer3D: Only explicit edge selection is supported. "
      "Use EdgeSelector::explicit_() with specific edge coordinates.");
}

auto
Chamfer3D::findHalfedge(const Surface_mesh_3  &mesh,
                        const EdgeIdentifier &edge) const
    -> Surface_mesh_3::Halfedge_index
{
  // Use a larger tolerance for matching coordinates
  constexpr double EDGE_TOLERANCE = 1e-6;

  // First pass: try exact endpoint matching
  for (auto he : mesh.halfedges()) {
    auto v1 = mesh.source(he);
    auto v2 = mesh.target(he);

    if ((pointsApproxEqual(mesh.point(v1), edge.start, EDGE_TOLERANCE) &&
         pointsApproxEqual(mesh.point(v2), edge.end, EDGE_TOLERANCE)) ||
        (pointsApproxEqual(mesh.point(v1), edge.end, EDGE_TOLERANCE) &&
         pointsApproxEqual(mesh.point(v2), edge.start, EDGE_TOLERANCE))) {
      return he;
    }
  }

  // Second pass: find edges with same direction that lie on the original edge line
  // This handles cases where previous chamfers have shortened the edge
  Kernel::Vector_3 edgeDir = edge.end - edge.start;
  double           edgeLen = std::sqrt(CGAL::to_double(edgeDir.squared_length()));
  if (edgeLen < EPSILON) {
    return Surface_mesh_3::null_halfedge();
  }
  edgeDir = edgeDir / edgeLen;

  Surface_mesh_3::Halfedge_index bestHe = Surface_mesh_3::null_halfedge();
  double                         bestScore = std::numeric_limits<double>::max();

  for (auto he : mesh.halfedges()) {
    auto            v1  = mesh.source(he);
    auto            v2  = mesh.target(he);
    Kernel::Point_3 p1  = mesh.point(v1);
    Kernel::Point_3 p2  = mesh.point(v2);

    Kernel::Vector_3 heDir = p2 - p1;
    double           heLen = std::sqrt(CGAL::to_double(heDir.squared_length()));
    if (heLen < EPSILON) {
      continue;
    }
    heDir = heDir / heLen;

    // Check if directions are parallel (dot product ~= 1 or -1)
    double dot = std::abs(CGAL::to_double(edgeDir * heDir));
    if (dot < 0.99) {
      continue; // Not parallel
    }

    // Check if midpoint of this edge lies on the original edge line
    Kernel::Point_3  midPt    = CGAL::midpoint(p1, p2);
    Kernel::Vector_3 toMid    = midPt - edge.start;
    double           proj     = CGAL::to_double(toMid * edgeDir);

    // Distance from midpoint to the original edge line
    Kernel::Vector_3 projVec  = edgeDir * proj;
    Kernel::Vector_3 perpVec  = toMid - projVec;
    double           perpDist = std::sqrt(CGAL::to_double(perpVec.squared_length()));

    // Must be close to the original edge line
    if (perpDist > edgeLen * 0.5) {
      continue; // Too far from original edge line
    }

    // Score by perpendicular distance (lower is better)
    if (perpDist < bestScore) {
      bestScore = perpDist;
      bestHe    = he;
    }
  }

  return bestHe;
}

auto
Chamfer3D::findSharpEdges(double angleThreshold) const
    -> std::vector<EdgeIdentifier>
{
  // Get all original edges first
  std::vector<EdgeIdentifier> originalEdges = extractOriginalEdges();
  std::vector<EdgeIdentifier> sharpEdges;

  // For each original edge, compute dihedral angle and filter
  for (const auto& edge : originalEdges) {
    double angle = computeDihedralAngle(edge);
    if (angle < 0) {
      continue; // Boundary edge or not found
    }

    // A "sharp" edge has a dihedral angle that deviates significantly from flat (PI)
    double deviation = std::abs(angle - M_PI);
    if (deviation > angleThreshold) {
      sharpEdges.push_back(edge);
    }
  }
  return sharpEdges;
}

auto
Chamfer3D::computeDihedralAngle(const EdgeIdentifier &edge) const -> double
{
  Surface_mesh_3 mesh = toSurfaceMesh();

  auto he = findHalfedge(mesh, edge);
  if (he == Surface_mesh_3::null_halfedge()) {
    return -1.0;
  }

  auto opp = mesh.opposite(he);
  if (mesh.is_border(he) || mesh.is_border(opp)) {
    return -1.0; // Boundary edge
  }

  auto f1 = mesh.face(he);
  auto f2 = mesh.face(opp);

  // Compute face normals
  auto fnormals =
      mesh.add_property_map<Surface_mesh_3::Face_index, Kernel::Vector_3>(
              "f:normals", CGAL::NULL_VECTOR)
          .first;

  PMP::compute_face_normals(mesh, fnormals);

  return computeDihedralAngleFromNormals(fnormals[f1], fnormals[f2]);
}

auto
Chamfer3D::chamferEdgeDirect(Surface_mesh_3                &mesh,
                             Surface_mesh_3::Halfedge_index he,
                             const ChamferParameters       &params) const -> bool
{
  if (he == Surface_mesh_3::null_halfedge()) {
    return false;
  }

  auto opp = mesh.opposite(he);

  // Skip boundary edges
  if (mesh.is_border(he) || mesh.is_border(opp)) {
    return false;
  }

  // Get vertices of the edge
  auto v1 = mesh.source(he);
  auto v2 = mesh.target(he);

  Kernel::Point_3 p1 = mesh.point(v1);
  Kernel::Point_3 p2 = mesh.point(v2);

  // Get adjacent faces
  auto f1 = mesh.face(he);
  auto f2 = mesh.face(opp);

  // Compute face normals
  auto fnormals =
      mesh.add_property_map<Surface_mesh_3::Face_index, Kernel::Vector_3>(
              "f:normals", CGAL::NULL_VECTOR)
          .first;
  PMP::compute_face_normals(mesh, fnormals);

  Kernel::Vector_3 n1 = fnormals[f1];
  Kernel::Vector_3 n2 = fnormals[f2];

  // Skip essentially flat edges (coplanar faces from triangulation)
  double dihedralAngle = computeDihedralAngleFromNormals(n1, n2);
  double deviation     = std::abs(dihedralAngle - M_PI);
  if (deviation < EPSILON) {
    return true; // Skip flat edges but report success (nothing to chamfer)
  }

  // Edge direction
  Kernel::Vector_3 edgeDir = normalize(p2 - p1);

  // Compute perpendicular directions in each face plane
  // p1 is inward direction perpendicular to edge in face 1
  Kernel::Vector_3 perp1 = normalize(CGAL::cross_product(edgeDir, n1));
  // p2 is inward direction perpendicular to edge in face 2
  Kernel::Vector_3 perp2 = normalize(CGAL::cross_product(n2, edgeDir));

  double d1 = params.distance1;
  double d2 = (params.type == ChamferParameters::Type::SYMMETRIC) ? params.distance1
                                                                   : params.distance2;

  // Compute chamfer vertices
  // c1_start and c1_end are on the face 1 side
  // c2_start and c2_end are on the face 2 side
  Kernel::Point_3 c1_start = p1 + perp1 * d1;
  Kernel::Point_3 c1_end   = p2 + perp1 * d1;
  Kernel::Point_3 c2_start = p1 + perp2 * d2;
  Kernel::Point_3 c2_end   = p2 + perp2 * d2;

  // Add new vertices for chamfer
  auto vc1_start = mesh.add_vertex(c1_start);
  auto vc1_end   = mesh.add_vertex(c1_end);
  auto vc2_start = mesh.add_vertex(c2_start);
  auto vc2_end   = mesh.add_vertex(c2_end);

  // Create the chamfer face (quadrilateral)
  // Note: The face orientation should match the mesh convention
  std::vector<Surface_mesh_3::Vertex_index> chamferFaceVerts = {
      vc1_start, vc1_end, vc2_end, vc2_start};

  auto chamferFace = mesh.add_face(chamferFaceVerts);

  if (chamferFace == Surface_mesh_3::null_face()) {
    // Failed to add face - remove added vertices
    mesh.remove_vertex(vc1_start);
    mesh.remove_vertex(vc1_end);
    mesh.remove_vertex(vc2_start);
    mesh.remove_vertex(vc2_end);
    return false;
  }

  // The complex part: modifying existing faces to incorporate chamfer vertices
  // This requires topology changes that are challenging with Surface_mesh
  // For now, we'll mark this as requiring the boolean fallback for complex
  // topology changes

  // Clean up - we successfully added the chamfer face but didn't complete
  // the topology modification. Mark for fallback.
  // In a complete implementation, we would:
  // 1. Split face f1 to include edge from vc1_start to vc1_end
  // 2. Split face f2 to include edge from vc2_start to vc2_end
  // 3. Remove the original edge and reconnect

  // For this implementation, return false to trigger boolean fallback
  // The chamfer face vertices are kept as they will be useful for the
  // boolean approach

  mesh.remove_face(chamferFace);
  mesh.remove_vertex(vc1_start);
  mesh.remove_vertex(vc1_end);
  mesh.remove_vertex(vc2_start);
  mesh.remove_vertex(vc2_end);

  return false; // Trigger boolean fallback
}

auto
Chamfer3D::createChamferWedge(const EdgeIdentifier    &edge,
                              const ChamferParameters &params,
                              const Kernel::Vector_3  &n1,
                              const Kernel::Vector_3  &n2) const
    -> std::unique_ptr<Solid>
{
  // Edge direction (normalized)
  Kernel::Vector_3 edgeDir = normalize(edge.end - edge.start);

  // Compute inward-facing perpendicular directions in each face plane
  // n1 and n2 are OUTWARD normals, so we need to negate the cross products
  // to get vectors pointing INTO the solid (away from the edge toward solid interior)
  Kernel::Vector_3 perp1 = -normalize(CGAL::cross_product(edgeDir, n1));
  Kernel::Vector_3 perp2 = -normalize(CGAL::cross_product(n2, edgeDir));

  double d1 = params.distance1;
  double d2 = (params.type == ChamferParameters::Type::SYMMETRIC) ? params.distance1
                                                                   : params.distance2;

  // Extend edge slightly beyond endpoints for robust boolean
  double           edgeLen   = std::sqrt(CGAL::to_double(CGAL::squared_distance(edge.start, edge.end)));
  double           extension = std::max(edgeLen * 0.05, d1 * 0.5); // 5% or half d1
  Kernel::Point_3  extStart  = edge.start - edgeDir * extension;
  Kernel::Point_3  extEnd    = edge.end + edgeDir * extension;

  // The wedge is a triangular prism that cuts the corner
  // Its cross-section is a triangle with:
  //   - One vertex at the edge (or slightly beyond to ensure clean cut)
  //   - One vertex offset by d1 in perp1 direction (on face 1)
  //   - One vertex offset by d2 in perp2 direction (on face 2)

  // Small extension beyond edge for robust boolean operation
  double edgeExtension = std::min(d1, d2) * 0.1;

  // Compute the chamfer plane normal (perpendicular to chamfer face)
  // The chamfer face connects the two cut lines on the adjacent faces
  Kernel::Vector_3 chamferNormal = normalize(perp1 * d2 + perp2 * d1);

  // Wedge vertices at start of edge
  // s0 extends slightly beyond the edge in the direction opposite to chamfer
  Kernel::Point_3 s0 = extStart - chamferNormal * edgeExtension;
  Kernel::Point_3 s1 = extStart + perp1 * d1;       // Offset by d1 on face 1
  Kernel::Point_3 s2 = extStart + perp2 * d2;       // Offset by d2 on face 2

  // Wedge vertices at end of edge
  Kernel::Point_3 e0 = extEnd - chamferNormal * edgeExtension;
  Kernel::Point_3 e1 = extEnd + perp1 * d1;
  Kernel::Point_3 e2 = extEnd + perp2 * d2;

  // Create triangular prism solid
  // All faces must have counterclockwise winding when viewed from outside (outward normals)
  auto ps = std::make_unique<PolyhedralSurface>();

  // Front triangle face (at start, normal toward -edgeDir)
  // Viewed from start looking toward end, counterclockwise: s0 -> s1 -> s2
  {
    auto *ring = new LineString();
    ring->addPoint(Point(s0));
    ring->addPoint(Point(s1));
    ring->addPoint(Point(s2));
    ring->addPoint(Point(s0));
    ps->addPatch(Polygon(ring));
  }

  // Back triangle face (at end, normal toward +edgeDir)
  // Viewed from end looking toward start, counterclockwise: e0 -> e2 -> e1
  {
    auto *ring = new LineString();
    ring->addPoint(Point(e0));
    ring->addPoint(Point(e2));
    ring->addPoint(Point(e1));
    ring->addPoint(Point(e0));
    ps->addPatch(Polygon(ring));
  }

  // Face 1 side (s0-s1-e1-e0): normal toward -perp2 direction
  {
    auto *ring = new LineString();
    ring->addPoint(Point(s0));
    ring->addPoint(Point(e0));
    ring->addPoint(Point(e1));
    ring->addPoint(Point(s1));
    ring->addPoint(Point(s0));
    ps->addPatch(Polygon(ring));
  }

  // Face 2 side (s0-s2-e2-e0): normal toward -perp1 direction
  {
    auto *ring = new LineString();
    ring->addPoint(Point(s0));
    ring->addPoint(Point(s2));
    ring->addPoint(Point(e2));
    ring->addPoint(Point(e0));
    ring->addPoint(Point(s0));
    ps->addPatch(Polygon(ring));
  }

  // Chamfer face (s1-s2-e2-e1): normal toward chamfer direction (outward from wedge)
  {
    auto *ring = new LineString();
    ring->addPoint(Point(s1));
    ring->addPoint(Point(e1));
    ring->addPoint(Point(e2));
    ring->addPoint(Point(s2));
    ring->addPoint(Point(s1));
    ps->addPatch(Polygon(ring));
  }

  return std::make_unique<Solid>(*ps);
}

auto
Chamfer3D::chamferEdgeBoolean(const Geometry          &geom,
                              const EdgeIdentifier    &edge,
                              const ChamferParameters &params) const
    -> std::unique_ptr<Geometry>
{
  // Handle GeometryCollection by extracting the first solid-like geometry
  const Geometry *workGeom = &geom;
  if (geom.is<GeometryCollection>()) {
    const auto &gc = geom.as<GeometryCollection>();
    for (size_t i = 0; i < gc.numGeometries(); ++i) {
      const auto &g = gc.geometryN(i);
      if (g.is<Solid>() || g.is<PolyhedralSurface>() ||
          g.is<TriangulatedSurface>()) {
        workGeom = &g;
        break;
      }
    }
    if (workGeom == &geom) {
      // No solid-like geometry found in collection
      return geom.clone();
    }
  }

  // Convert input geometry to Surface_mesh for plane clipping
  Surface_mesh_3 mesh;
  if (workGeom->is<Solid>()) {
    mesh = workGeom->as<Solid>().toSurfaceMesh();
  } else if (workGeom->is<PolyhedralSurface>()) {
    mesh = workGeom->as<PolyhedralSurface>().toSurfaceMesh();
  } else if (workGeom->is<TriangulatedSurface>()) {
    mesh = workGeom->as<TriangulatedSurface>().toSurfaceMesh();
  } else {
    throw std::runtime_error(
        "Chamfer3D: Input geometry must be Solid or PolyhedralSurface, got: " +
        std::to_string(static_cast<int>(geom.geometryTypeId())));
  }

  // Repair mesh immediately after conversion
  // The toSurfaceMesh conversion may create boundary edges if faces don't
  // properly share vertices (e.g., due to floating point precision issues)
  PMP::stitch_borders(mesh);
  PMP::merge_duplicated_vertices_in_boundary_cycles(mesh);
  mesh.collect_garbage();

  // Get face normals at the edge
  auto he = findHalfedge(mesh, edge);
  if (he == Surface_mesh_3::null_halfedge()) {
    // Edge not found - skip and continue with current geometry
    return workGeom->clone();
  }

  auto opp = mesh.opposite(he);
  if (mesh.is_border(he) || mesh.is_border(opp)) {
    // Boundary edge - nothing to chamfer (no two faces meet)
    return workGeom->clone();
  }

  auto fnormals =
      mesh.add_property_map<Surface_mesh_3::Face_index, Kernel::Vector_3>(
              "f:normals", CGAL::NULL_VECTOR)
          .first;
  PMP::compute_face_normals(mesh, fnormals);

  Kernel::Vector_3 n1 = fnormals[mesh.face(he)];
  Kernel::Vector_3 n2 = fnormals[mesh.face(opp)];

  // Check if the found halfedge direction matches the original edge direction
  // If they point in opposite directions, swap the face normals
  Kernel::Vector_3 heDir = mesh.point(mesh.target(he)) - mesh.point(mesh.source(he));
  Kernel::Vector_3 reqDir = edge.end - edge.start;
  if (CGAL::to_double(heDir * reqDir) < 0) {
    // Halfedge is in opposite direction - swap normals so they match the expected orientation
    std::swap(n1, n2);
  }

  // Skip essentially flat edges (coplanar faces from triangulation)
  double dihedralAngle = computeDihedralAngleFromNormals(n1, n2);
  double deviation     = std::abs(dihedralAngle - M_PI);
  if (deviation < EPSILON) {
    // Skip flat edges (nothing to chamfer)
    return workGeom->clone();
  }


  // Chamfer distances
  double d1 = params.distance1;
  double d2 = (params.type == ChamferParameters::Type::SYMMETRIC)
                  ? params.distance1
                  : params.distance2;

  // Get ACTUAL edge coordinates from the mesh
  // This is critical for multi-edge chamfers where previous chamfers have modified vertices
  Kernel::Point_3 heStart = mesh.point(mesh.source(he));
  Kernel::Point_3 heEnd = mesh.point(mesh.target(he));

  // Use ORIGINAL edge direction for consistency in perpendicular calculations
  // This ensures the chamfer is computed with the intended orientation
  Kernel::Vector_3 originalDir = edge.end - edge.start;
  double           originalLen = std::sqrt(CGAL::to_double(originalDir.squared_length()));
  if (originalLen < EPSILON) {
    return workGeom->clone();
  }
  Kernel::Vector_3 edgeDir = originalDir / originalLen;

  // Pick actual mesh point for chamfer positioning
  // Use heStart (which lies on the current mesh edge) as the reference point
  Kernel::Point_3 actualStart = heStart;

  // Compute perpendicular directions in each face plane
  // cross(edgeDir, n) gives a vector perpendicular to edge, lying in the face plane
  Kernel::Vector_3 perp1_raw = CGAL::cross_product(edgeDir, n1);
  Kernel::Vector_3 perp2_raw = CGAL::cross_product(n2, edgeDir);

  // Ensure perpendiculars point INTO the solid (toward the dihedral angle interior)
  // For a convex edge, perp1 should point somewhat toward face 2's interior (same direction as -n2)
  // For a convex edge, perp2 should point somewhat toward face 1's interior (same direction as -n1)
  Kernel::Vector_3 perp1 = (CGAL::to_double(perp1_raw * (-n2)) >= 0) ? perp1_raw : -perp1_raw;
  Kernel::Vector_3 perp2 = (CGAL::to_double(perp2_raw * (-n1)) >= 0) ? perp2_raw : -perp2_raw;

  perp1 = normalize(perp1);
  perp2 = normalize(perp2);

  // Chamfer cut points on the edge (using actual start point, not original)
  // Point P1: on face F1, at distance d1 from edge along perp1
  // Point P2: on face F2, at distance d2 from edge along perp2
  Kernel::Point_3 chamferPt1 = actualStart + perp1 * d1;
  Kernel::Point_3 chamferPt2 = actualStart + perp2 * d2;

  // The chamfer plane passes through the two chamfer lines
  // Normal is computed from edge direction and the line connecting chamfer points
  Kernel::Vector_3 chamferDir = chamferPt2 - chamferPt1;
  Kernel::Vector_3 chamferNormal = normalize(CGAL::cross_product(edgeDir, chamferDir));

  // PMP::clip keeps the NEGATIVE half-space (where normal points away from)
  // We want to REMOVE the corner (at actualStart) and KEEP the solid body
  // So the normal should point TOWARD the corner to remove
  Kernel::Vector_3 toCorner = actualStart - chamferPt1;
  if (CGAL::to_double(chamferNormal * toCorner) < 0) {
    // Normal points away from corner, flip it to point toward corner
    chamferNormal = -chamferNormal;
  }

  // Create chamfer plane through chamferPt1 with computed normal
  Kernel::Plane_3 chamferPlane(chamferPt1, chamferNormal);

  // SAFETY CHECK: Detect if chamfer plane would clip non-adjacent vertices
  // This happens for internal/concave edges where the infinite plane intersects
  // other parts of the geometry. In such cases, skip this edge.
  {
    // Collect vertices adjacent to this edge (in both faces)
    std::set<Surface_mesh_3::Vertex_index> adjacentVerts;

    // Add vertices from face containing he
    for (auto v : vertices_around_face(he, mesh)) {
      adjacentVerts.insert(v);
    }
    // Add vertices from face containing opp (opposite halfedge)
    for (auto v : vertices_around_face(opp, mesh)) {
      adjacentVerts.insert(v);
    }

    // Check all mesh vertices not adjacent to the edge
    // If any are on the "removed" side of the plane (positive side where normal points),
    // then this chamfer would damage distant parts of the geometry
    bool wouldClipNonAdjacent = false;
    for (auto v : mesh.vertices()) {
      if (adjacentVerts.count(v) > 0) {
        continue; // Skip adjacent vertices
      }

      Kernel::Point_3  pt     = mesh.point(v);
      Kernel::Vector_3 toVert = pt - chamferPt1;
      double           signedDist = CGAL::to_double(chamferNormal * toVert);

      // If vertex is on the positive side (will be clipped) by more than a small tolerance
      // relative to chamfer distance, this edge should be skipped
      // Use a larger threshold to avoid false positives from nearby chamfered vertices
      double threshold = std::max(d1, d2) * 2.0;
      if (signedDist > threshold) {
        wouldClipNonAdjacent = true;
        break;
      }
    }

    if (wouldClipNonAdjacent) {
      // Determine if edge is truly concave using the perpendicular cross-product test
      // Use edgeDir (from original user-specified edge) for consistent direction
      bool isTrueConcave = isConcaveEdge(perp1, perp2, edgeDir);

      if (isTrueConcave) {
        // This is a CONCAVE edge - ADD material by creating a fill wedge
        // For concave edges, perp1/perp2 point INTO the solid
        // For the fill wedge, we need them pointing INTO THE VOID (negate)
        auto fillWedge = createConcaveFillWedge(
            heStart, heEnd, edgeDir, -perp1, -perp2, d1, d2);

        if (fillWedge && !fillWedge->isEmpty()) {
          try {
            auto unionResult = union3D(*workGeom, *fillWedge);
            if (unionResult && !unionResult->isEmpty()) {
              return unionResult;
            }
          } catch (...) {
            // Union failed
          }
        }
      } else {
        // This is a CONVEX edge where infinite plane clipping would fail
        // Use BOUNDED WEDGE SUBTRACTION with actual mesh endpoints
        // (no extension beyond endpoints to avoid overlapping at corners)
        auto subtractWedge = createBoundedSubtractionWedge(
            heStart, heEnd, edgeDir, perp1, perp2, d1, d2);
        if (subtractWedge && !subtractWedge->isEmpty()) {
          try {
            auto diffResult = difference3D(*workGeom, *subtractWedge);
            if (diffResult && !diffResult->isEmpty()) {
              return diffResult;
            }
          } catch (...) {
            // Difference failed
          }
        }
      }
      // Fallback: skip this edge
      return workGeom->clone();
    }
  }

  // Triangulate mesh for clipping (required by PMP::clip)
  PMP::triangulate_faces(mesh);

  // IMPORTANT: Repair any existing boundary edges BEFORE clipping
  // This prevents accumulation of boundaries that create holes
  PMP::stitch_borders(mesh);
  PMP::merge_duplicated_vertices_in_boundary_cycles(mesh);

  // Clip mesh with chamfer plane - removes the corner material
  bool clip_ok =
      PMP::clip(mesh, chamferPlane, CGAL::parameters::clip_volume(true));
  if (!clip_ok || mesh.is_empty()) {
    return workGeom->clone();
  }

  // IMPORTANT: Remove degenerate faces created by clipping
  // PMP::clip can create zero-area triangles at intersection boundaries
  removeDegenerateFaces(mesh);

  // Repair mesh topology - stitch borders and merge duplicate vertices
  // This is crucial for multi-edge chamfers where successive clips can create
  // small gaps that need to be stitched together
  PMP::stitch_borders(mesh);
  PMP::merge_duplicated_vertices_in_boundary_cycles(mesh);

  // Fill any remaining holes in the mesh
  // This is critical for multi-edge chamfers where the clip operation may leave holes
  std::vector<Surface_mesh_3::Halfedge_index> borderCycles;
  PMP::extract_boundary_cycles(mesh, std::back_inserter(borderCycles));
  for (auto borderHe : borderCycles) {
    try {
      PMP::triangulate_hole(mesh, borderHe);
    } catch (...) {
      // Hole filling failed (e.g., self-intersecting boundary)
      // Continue with other holes - the mesh may have some unfilled holes
    }
  }

  // Try to remove degenerate faces and edges
  PMP::remove_degenerate_faces(mesh);
  PMP::remove_degenerate_edges(mesh);

  // Collect garbage to remove orphaned vertices/faces
  mesh.collect_garbage();

  // Convert result back to SFCGAL geometry
  auto ps = std::make_unique<PolyhedralSurface>(mesh);

  if (_isSolid) {
    return std::make_unique<Solid>(*ps);
  }
  return ps;
}

auto
Chamfer3D::chamferVertexDirect(Surface_mesh_3                &mesh,
                               Surface_mesh_3::Vertex_index   v,
                               const VertexChamferParameters &params) const
    -> bool
{
  if (v == Surface_mesh_3::null_vertex()) {
    return false;
  }

  Kernel::Point_3 vertexPos = mesh.point(v);

  // Collect incident halfedges and compute cut points
  std::vector<Surface_mesh_3::Halfedge_index> incidentHalfedges;
  std::vector<Kernel::Point_3>                cutPoints;

  auto startHe = mesh.halfedge(v);
  if (startHe == Surface_mesh_3::null_halfedge()) {
    return false;
  }

  for (auto he : mesh.halfedges_around_target(startHe)) {
    incidentHalfedges.push_back(he);

    auto           source    = mesh.source(he);
    Kernel::Point_3 sourcePos = mesh.point(source);

    // Direction from vertex to source (along the edge)
    Kernel::Vector_3 dir = sourcePos - vertexPos;
    double           len = std::sqrt(CGAL::to_double(dir.squared_length()));

    if (len < EPSILON) {
      continue;
    }

    // Cut point at distance from vertex
    double          cutDist  = std::min(params.distance, len * 0.9); // Don't exceed 90%
    Kernel::Point_3 cutPoint = vertexPos + normalize(dir) * cutDist;
    cutPoints.push_back(cutPoint);
  }

  if (cutPoints.size() < 3) {
    return false; // Need at least 3 edges to form a chamfer face
  }

  // Add cut point vertices
  std::vector<Surface_mesh_3::Vertex_index> cutVertices;
  for (const auto &pt : cutPoints) {
    cutVertices.push_back(mesh.add_vertex(pt));
  }

  // Create chamfer face from cut points
  auto chamferFace = mesh.add_face(cutVertices);

  if (chamferFace == Surface_mesh_3::null_face()) {
    // Failed - clean up
    for (auto cv : cutVertices) {
      mesh.remove_vertex(cv);
    }
    return false;
  }

  // Similar to edge chamfer, complete topology modification is complex
  // Return false to use fallback approach
  mesh.remove_face(chamferFace);
  for (auto cv : cutVertices) {
    mesh.remove_vertex(cv);
  }

  return false;
}

auto
Chamfer3D::chamferEdges(const EdgeSelector      &selector,
                        const ChamferParameters &params) const
    -> std::unique_ptr<Geometry>
{
  if (params.distance1 <= 0 ||
      (params.type == ChamferParameters::Type::ASYMMETRIC &&
       params.distance2 <= 0)) {
    throw std::invalid_argument("Chamfer distance must be positive");
  }

  std::vector<EdgeIdentifier> edges = resolveEdges(selector);

  if (edges.empty()) {
    return _geometry->clone();
  }

  // Classify edges as concave or convex based on the perpendicular cross-product test
  // Process CONCAVE edges FIRST (fill wedges), then CONVEX edges (subtraction)
  // This ensures fill wedges establish corner geometry before subtractions
  std::vector<EdgeIdentifier> concaveEdges;
  std::vector<EdgeIdentifier> convexEdges;
  {
    Surface_mesh_3 mesh = toSurfaceMesh();
    PMP::stitch_borders(mesh);

    auto fnormals =
        mesh.add_property_map<Surface_mesh_3::Face_index, Kernel::Vector_3>(
                "f:normals", CGAL::NULL_VECTOR)
            .first;
    PMP::compute_face_normals(mesh, fnormals);

    for (const auto &edge : edges) {
      auto he = findHalfedge(mesh, edge);
      if (he == Surface_mesh_3::null_halfedge() ||
          mesh.is_border(he) || mesh.is_border(mesh.opposite(he))) {
        convexEdges.push_back(edge);  // Default to convex for unmatched edges
        continue;
      }

      Kernel::Vector_3 n1 = fnormals[mesh.face(he)];
      Kernel::Vector_3 n2 = fnormals[mesh.face(mesh.opposite(he))];

      // Check if halfedge direction matches original edge direction
      Kernel::Vector_3 heDir = mesh.point(mesh.target(he)) - mesh.point(mesh.source(he));
      Kernel::Vector_3 reqDir = edge.end - edge.start;
      if (CGAL::to_double(heDir * reqDir) < 0) {
        std::swap(n1, n2);
      }

      Kernel::Vector_3 edgeDir = normalize(edge.end - edge.start);

      // Concavity test using cross product of face normals
      // Cross product of normals gives the edge direction (for convex) or opposite (for concave)
      // For convex edge: cross(n1,n2) · edgeDir > 0 (right-hand rule from n1 to n2)
      // For concave edge: cross(n1,n2) · edgeDir < 0
      Kernel::Vector_3 normalCross = CGAL::cross_product(n1, n2);
      double crossDotEdge = CGAL::to_double(normalCross * edgeDir);
      bool isConcave = crossDotEdge < -0.01;

      if (isConcave) {
        concaveEdges.push_back(edge);
      } else {
        convexEdges.push_back(edge);
      }
    }
  }

  // Helper to create string key for vertex coordinates (with tolerance)
  auto vertexKey = [](const Kernel::Point_3 &p) {
    auto round = [](double v) { return std::round(v * 10000.0) / 10000.0; };
    std::ostringstream oss;
    oss << round(CGAL::to_double(p.x())) << ","
        << round(CGAL::to_double(p.y())) << ","
        << round(CGAL::to_double(p.z()));
    return oss.str();
  };

  // Helper to check if edge is vertical (Z-aligned)
  auto isVerticalEdge = [](const EdgeIdentifier &e) {
    double dx = std::abs(CGAL::to_double(e.end.x() - e.start.x()));
    double dy = std::abs(CGAL::to_double(e.end.y() - e.start.y()));
    double dz = std::abs(CGAL::to_double(e.end.z() - e.start.z()));
    return dz > dx && dz > dy;
  };

  // Helper to get edge Z level (for horizontal edges)
  auto getEdgeZ = [](const EdgeIdentifier &e) {
    return std::min(CGAL::to_double(e.start.z()), CGAL::to_double(e.end.z()));
  };

  // Sort convex edges: vertical first, then horizontal (bottom before top)
  // This improves chamfer quality at corners where vertical and horizontal edges meet
  std::stable_sort(convexEdges.begin(), convexEdges.end(),
                   [&isVerticalEdge, &getEdgeZ](const EdgeIdentifier &a, const EdgeIdentifier &b) {
                     bool aVert = isVerticalEdge(a);
                     bool bVert = isVerticalEdge(b);
                     if (aVert != bVert) return aVert; // Vertical first
                     if (!aVert && !bVert) {
                       // Both horizontal: bottom (lower Z) before top
                       return getEdgeZ(a) < getEdgeZ(b);
                     }
                     return false;
                   });

  // BATCH PROCESSING: Create all wedges on ORIGINAL geometry, then single boolean
  // This avoids accumulated numerical errors from sequential boolean operations

  // Get mesh and normals from ORIGINAL geometry (only once)
  Surface_mesh_3 originalMesh = toSurfaceMesh();
  PMP::stitch_borders(originalMesh);

  auto fnormals =
      originalMesh.add_property_map<Surface_mesh_3::Face_index, Kernel::Vector_3>(
              "f:normals", CGAL::NULL_VECTOR)
          .first;
  PMP::compute_face_normals(originalMesh, fnormals);

  // Collect all fill wedges (for concave edges) and subtraction wedges (for convex edges)
  std::vector<std::unique_ptr<Solid>> fillWedges;
  std::vector<std::unique_ptr<Solid>> subtractionWedges;

  double d1 = params.distance1;
  double d2 = (params.type == ChamferParameters::Type::SYMMETRIC) ? params.distance1 : params.distance2;

  // NOTE: Concave (reflex) edges are automatically SKIPPED
  // Chamfering concave edges requires adding material (fill wedges) which doesn't
  // join cleanly with the subtraction wedges of adjacent convex edges at corners.
  // For clean chamfer results, only convex edges are processed.

  // =========================================================================
  // PRE-COMPUTE BISECTOR PLANES FOR INNER CORNERS (2 convex + concave)
  // This is done BEFORE wedge creation so we can apply miter joints
  // =========================================================================

  // Map vertex key -> bisector plane (for inner corners only)
  std::map<std::string, Kernel::Plane_3> innerCornerBisectorPlanes;
  {
    // Build vertex -> edges maps
    std::map<std::string, std::vector<size_t>> vertexToConvexEdges;
    for (size_t i = 0; i < convexEdges.size(); ++i) {
      vertexToConvexEdges[vertexKey(convexEdges[i].start)].push_back(i);
      vertexToConvexEdges[vertexKey(convexEdges[i].end)].push_back(i);
    }

    std::map<std::string, std::vector<size_t>> vertexToConcaveEdges;
    for (size_t i = 0; i < concaveEdges.size(); ++i) {
      vertexToConcaveEdges[vertexKey(concaveEdges[i].start)].push_back(i);
      vertexToConcaveEdges[vertexKey(concaveEdges[i].end)].push_back(i);
    }

    // Find inner corners: exactly 2 convex edges + at least 1 concave edge
    for (const auto &[vkey, convexEdgeIndices] : vertexToConvexEdges) {
      if (convexEdgeIndices.size() != 2) continue;

      // Check for concave edge at this vertex
      auto concaveIt = vertexToConcaveEdges.find(vkey);
      if (concaveIt == vertexToConcaveEdges.end()) continue;
      if (concaveIt->second.empty()) continue;

      // This is an inner corner - compute bisector plane
      const EdgeIdentifier &edge1 = convexEdges[convexEdgeIndices[0]];
      const EdgeIdentifier &edge2 = convexEdges[convexEdgeIndices[1]];

      // Find the corner vertex
      Kernel::Point_3 cornerVertex;
      if (vertexKey(edge1.start) == vkey) {
        cornerVertex = edge1.start;
      } else {
        cornerVertex = edge1.end;
      }

      // Compute edge directions (pointing AWAY from corner)
      Kernel::Vector_3 edgeDir1, edgeDir2;
      if (vertexKey(edge1.start) == vkey) {
        edgeDir1 = normalize(edge1.end - edge1.start);
      } else {
        edgeDir1 = normalize(edge1.start - edge1.end);
      }
      if (vertexKey(edge2.start) == vkey) {
        edgeDir2 = normalize(edge2.end - edge2.start);
      } else {
        edgeDir2 = normalize(edge2.start - edge2.end);
      }

      // Compute and store bisector plane
      Kernel::Plane_3 bisectorPlane = computeBisectorPlane(cornerVertex, edgeDir1, edgeDir2);
      innerCornerBisectorPlanes[vkey] = bisectorPlane;
    }
  }

  // Process only CONVEX edges on the ORIGINAL mesh
  for (const auto &edge : convexEdges) {
    auto he = findHalfedge(originalMesh, edge);
    if (he == Surface_mesh_3::null_halfedge()) continue;
    if (originalMesh.is_border(he) || originalMesh.is_border(originalMesh.opposite(he))) continue;

    Kernel::Vector_3 n1 = fnormals[originalMesh.face(he)];
    Kernel::Vector_3 n2 = fnormals[originalMesh.face(originalMesh.opposite(he))];

    // Check halfedge direction
    Kernel::Vector_3 heDir = originalMesh.point(originalMesh.target(he)) - originalMesh.point(originalMesh.source(he));
    if (CGAL::to_double(heDir * (edge.end - edge.start)) < 0) std::swap(n1, n2);

    // Skip flat edges
    double dihedralAngle = computeDihedralAngleFromNormals(n1, n2);
    if (std::abs(dihedralAngle - M_PI) < EPSILON) continue;

    Kernel::Point_3 heStart = originalMesh.point(originalMesh.source(he));
    Kernel::Point_3 heEnd = originalMesh.point(originalMesh.target(he));
    Kernel::Vector_3 edgeDir = normalize(edge.end - edge.start);

    Kernel::Vector_3 perp1 = normalize(CGAL::cross_product(edgeDir, n1));
    Kernel::Vector_3 perp2 = normalize(CGAL::cross_product(n2, edgeDir));
    if (CGAL::to_double(perp1 * (-n2)) < 0) perp1 = -perp1;
    if (CGAL::to_double(perp2 * (-n1)) < 0) perp2 = -perp2;

    // Check if this edge has inner corner endpoints that need miter joints
    std::string startKey = vertexKey(edge.start);
    std::string endKey = vertexKey(edge.end);

    auto startBisectorIt = innerCornerBisectorPlanes.find(startKey);
    auto endBisectorIt = innerCornerBisectorPlanes.find(endKey);

    const Kernel::Plane_3 *startBisector = (startBisectorIt != innerCornerBisectorPlanes.end())
                                            ? &startBisectorIt->second : nullptr;
    const Kernel::Plane_3 *endBisector = (endBisectorIt != innerCornerBisectorPlanes.end())
                                          ? &endBisectorIt->second : nullptr;

    // Create subtraction wedge (with optional miter joint clipping)
    std::unique_ptr<Solid> subWedge;
    if (startBisector != nullptr || endBisector != nullptr) {
      // Use mitered wedge for inner corner joints
      subWedge = createMiteredSubtractionWedge(heStart, heEnd, edgeDir, perp1, perp2, d1, d2,
                                                startBisector, endBisector);
    } else {
      // Use regular bounded wedge
      subWedge = createBoundedSubtractionWedge(heStart, heEnd, edgeDir, perp1, perp2, d1, d2);
    }

    if (subWedge && !subWedge->isEmpty()) {
      subtractionWedges.push_back(std::move(subWedge));
    }
  }

  // =========================================================================
  // CORNER VERTEX HANDLING
  // Detect vertices where convex edges meet and create corner clips
  // =========================================================================
  {
    // Map vertex -> list of convex edge indices sharing that vertex
    std::map<std::string, std::vector<size_t>> vertexToConvexEdges;
    for (size_t i = 0; i < convexEdges.size(); ++i) {
      vertexToConvexEdges[vertexKey(convexEdges[i].start)].push_back(i);
      vertexToConvexEdges[vertexKey(convexEdges[i].end)].push_back(i);
    }

    // Map vertex -> list of concave edge indices sharing that vertex
    std::map<std::string, std::vector<size_t>> vertexToConcaveEdges;
    for (size_t i = 0; i < concaveEdges.size(); ++i) {
      vertexToConcaveEdges[vertexKey(concaveEdges[i].start)].push_back(i);
      vertexToConcaveEdges[vertexKey(concaveEdges[i].end)].push_back(i);
    }

    // Process each vertex based on number of convex edges
    for (const auto &[vkey, convexEdgeIndices] : vertexToConvexEdges) {
      // Find the corner vertex
      Kernel::Point_3 cornerVertex;
      bool foundCorner = false;
      for (size_t idx : convexEdgeIndices) {
        const auto &edge = convexEdges[idx];
        if (vertexKey(edge.start) == vkey) {
          cornerVertex = edge.start;
          foundCorner = true;
          break;
        }
        if (vertexKey(edge.end) == vkey) {
          cornerVertex = edge.end;
          foundCorner = true;
          break;
        }
      }
      if (!foundCorner) continue;

      // Compute chamfer data (edgeDir, perp1, perp2) for each convex edge at this corner
      // This uses the same computation as the main edge processing loop
      struct EdgeChamferData {
        Kernel::Vector_3 edgeDir;
        Kernel::Vector_3 perp1;
        Kernel::Vector_3 perp2;
      };
      std::vector<EdgeChamferData> edgeChamferData;
      std::vector<Kernel::Point_3> chamferPoints;
      double dist = (d1 + d2) / 2.0;

      for (size_t idx : convexEdgeIndices) {
        const auto &edge = convexEdges[idx];

        // Find the halfedge in the mesh
        auto he = findHalfedge(originalMesh, edge);
        if (he == Surface_mesh_3::null_halfedge()) continue;
        if (originalMesh.is_border(he) || originalMesh.is_border(originalMesh.opposite(he))) continue;

        // Get face normals
        Kernel::Vector_3 n1 = fnormals[originalMesh.face(he)];
        Kernel::Vector_3 n2 = fnormals[originalMesh.face(originalMesh.opposite(he))];

        // Check halfedge direction matches edge direction
        Kernel::Vector_3 heDir = originalMesh.point(originalMesh.target(he)) - originalMesh.point(originalMesh.source(he));
        if (CGAL::to_double(heDir * (edge.end - edge.start)) < 0) std::swap(n1, n2);

        // Compute edge direction (pointing away from corner)
        Kernel::Vector_3 edgeVec;
        if (vertexKey(edge.start) == vkey) {
          edgeVec = edge.end - edge.start;
        } else {
          edgeVec = edge.start - edge.end;
        }
        double edgeLen = std::sqrt(CGAL::to_double(edgeVec.squared_length()));
        if (edgeLen < 1e-10) continue;
        Kernel::Vector_3 edgeDir = edgeVec * (1.0 / edgeLen);

        // Compute perp1 and perp2 (same as in main edge processing loop)
        Kernel::Vector_3 perp1 = normalize(CGAL::cross_product(edgeDir, n1));
        Kernel::Vector_3 perp2 = normalize(CGAL::cross_product(n2, edgeDir));
        if (CGAL::to_double(perp1 * (-n2)) < 0) perp1 = -perp1;
        if (CGAL::to_double(perp2 * (-n1)) < 0) perp2 = -perp2;

        edgeChamferData.push_back({edgeDir, perp1, perp2});
        chamferPoints.push_back(cornerVertex + edgeDir * dist);
      }

      if (convexEdgeIndices.size() == 3 && chamferPoints.size() == 3) {
        // 3-edge corner (external corner like cube vertex): create tetrahedron clip
        auto cornerClip = createCornerClipTetrahedron(cornerVertex, chamferPoints);
        if (cornerClip && !cornerClip->isEmpty()) {
          subtractionWedges.push_back(std::move(cornerClip));
        }
      }
      // NOTE: 2-edge inner corners (with concave edge) are now handled by
      // miter joints in the wedge creation, no tetrahedron clip needed.
    }
  }

  // Now perform batch boolean operations
  // Ensure we start with a Solid (not PolyhedralSurface) for proper boolean operations
  std::unique_ptr<Geometry> result;
  if (_isSolid) {
    result = _geometry->clone();
  } else {
    // Convert PolyhedralSurface to Solid
    const auto &ps = _geometry->as<PolyhedralSurface>();
    result = std::make_unique<Solid>(ps);
  }

  // Step 1: Union all fill wedges with original
  if (!fillWedges.empty()) {
    for (auto &fillWedge : fillWedges) {
      try {
        auto unionResult = union3D(*result, *fillWedge);
        if (unionResult && !unionResult->isEmpty()) {
          result = std::move(unionResult);
        }
      } catch (...) {
        // Individual fill wedge failed, continue with others
      }
    }
  }

  // Step 2: Subtract all subtraction wedges
  if (!subtractionWedges.empty()) {
    // Union all subtraction wedges into one geometry for single subtraction
    std::unique_ptr<Geometry> allSubtractions = subtractionWedges[0]->clone();
    for (size_t i = 1; i < subtractionWedges.size(); ++i) {
      try {
        auto unionResult = union3D(*allSubtractions, *subtractionWedges[i]);
        if (unionResult && !unionResult->isEmpty()) {
          allSubtractions = std::move(unionResult);
        }
      } catch (...) {
        // Individual union failed
      }
    }

    // Single difference operation
    try {
      auto diffResult = difference3D(*result, *allSubtractions);
      if (diffResult && !diffResult->isEmpty()) {
        result = std::move(diffResult);
      }
    } catch (...) {
      // Batch subtraction failed, fall back to individual subtractions
      for (auto &subWedge : subtractionWedges) {
        try {
          auto diffResult = difference3D(*result, *subWedge);
          if (diffResult && !diffResult->isEmpty()) {
            result = std::move(diffResult);
          }
        } catch (...) {
          // Individual subtraction failed
        }
      }
    }
  }

  // Clean degenerate faces from the result
  if (result && result->geometryTypeId() == TYPE_SOLID) {
    auto cleaned = cleanDegenerateFaces(static_cast<const Solid &>(*result));
    if (cleaned && !cleaned->isEmpty()) {
      return cleaned;
    }
  }

  return result;
}

auto
Chamfer3D::chamferVertices(const std::vector<Kernel::Point_3> &vertices,
                           const VertexChamferParameters      &params) const
    -> std::unique_ptr<Geometry>
{
  if (params.distance <= 0) {
    throw std::invalid_argument("Vertex chamfer distance must be positive");
  }

  if (vertices.empty()) {
    return _geometry->clone();
  }

  Surface_mesh_3 mesh = toSurfaceMesh();

  // Find vertices in mesh
  std::vector<Surface_mesh_3::Vertex_index> targetVertices;

  for (const auto &pt : vertices) {
    for (auto v : mesh.vertices()) {
      if (pointsApproxEqual(mesh.point(v), pt)) {
        targetVertices.push_back(v);
        break;
      }
    }
  }

  if (targetVertices.empty()) {
    return _geometry->clone();
  }

  // Try direct vertex chamfer
  bool allSucceeded = true;
  for (auto v : targetVertices) {
    if (!chamferVertexDirect(mesh, v, params)) {
      allSucceeded = false;
      break;
    }
  }

  if (allSucceeded) {
    return fromSurfaceMesh(mesh);
  }

  // Boolean fallback for vertex chamfer
  // Convert each vertex chamfer to edge chamfers on incident edges
  std::vector<EdgeIdentifier> incidentEdges;

  mesh = toSurfaceMesh(); // Reset mesh

  for (auto v : targetVertices) {
    auto startHe = mesh.halfedge(v);
    if (startHe == Surface_mesh_3::null_halfedge()) {
      continue;
    }
    for (auto he : mesh.halfedges_around_target(startHe)) {
      auto source = mesh.source(he);
      incidentEdges.emplace_back(mesh.point(v), mesh.point(source));
    }
  }

  // Apply small chamfer to all incident edges
  auto edgeParams = ChamferParameters::symmetric(params.distance * 0.5);
  auto edgeSelector = EdgeSelector::explicit_(incidentEdges);

  return chamferEdges(edgeSelector, edgeParams);
}

// ============================================================================
// Convenience Functions
// ============================================================================

auto
chamfer3D(const Geometry                    &geometry,
          const std::vector<EdgeIdentifier> &edges,
          const ChamferParameters           &params) -> std::unique_ptr<Geometry>
{
  Chamfer3D op(geometry);
  auto      selector = EdgeSelector::explicit_(edges);
  return op.chamferEdges(selector, params);
}

auto
chamfer3DVertices(const Geometry                     &geometry,
                  const std::vector<Kernel::Point_3> &vertices,
                  const VertexChamferParameters      &params)
    -> std::unique_ptr<Geometry>
{
  Chamfer3D op(geometry);
  return op.chamferVertices(vertices, params);
}

} // namespace SFCGAL::algorithm
