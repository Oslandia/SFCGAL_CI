// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_CHAMFER3D_ALGORITHM
#define SFCGAL_CHAMFER3D_ALGORITHM

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Surface_mesh.h>

#include <cmath>
#include <memory>
#include <vector>

namespace SFCGAL::algorithm {

/**
 * @brief Represents an edge identified by its two endpoint coordinates.
 *
 * Used to specify which edges should be chamfered in a 3D geometry.
 */
struct SFCGAL_API EdgeIdentifier {
  Kernel::Point_3 start; ///< Start point of the edge
  Kernel::Point_3 end;   ///< End point of the edge

  /**
   * @brief Default constructor
   */
  EdgeIdentifier() = default;

  /**
   * @brief Construct from two points
   * @param s Start point
   * @param e End point
   */
  EdgeIdentifier(const Kernel::Point_3 &s, const Kernel::Point_3 &e)
      : start(s), end(e)
  {
  }

  /**
   * @brief Construct from SFCGAL points
   * @param s Start point
   * @param e End point
   */
  EdgeIdentifier(const Point &s, const Point &e)
      : start(s.toPoint_3()), end(e.toPoint_3())
  {
  }

  /**
   * @brief Equality comparison (order-independent)
   * @param other Edge to compare with
   * @return true if edges represent the same geometric edge
   */
  auto
  operator==(const EdgeIdentifier &other) const -> bool
  {
    return (start == other.start && end == other.end) ||
           (start == other.end && end == other.start);
  }

  /**
   * @brief Less-than comparison for use in ordered containers
   * @param other Edge to compare with
   * @return true if this edge is less than other
   */
  auto
  operator<(const EdgeIdentifier &other) const -> bool
  {
    // Helper to compare two points lexicographically
    auto pointLess = [](const Kernel::Point_3 &a,
                        const Kernel::Point_3 &b) -> bool {
      if (a.x() != b.x()) {
        return a.x() < b.x();
      }
      if (a.y() != b.y()) {
        return a.y() < b.y();
      }
      return a.z() < b.z();
    };

    // Normalize so smaller point comes first
    auto [minThis, maxThis]   = std::minmax(start, end, pointLess);
    auto [minOther, maxOther] = std::minmax(other.start, other.end, pointLess);

    if (minThis != minOther) {
      return pointLess(minThis, minOther);
    }
    return pointLess(maxThis, maxOther);
  }
};

/**
 * @brief Edge selection criteria for chamfering operations.
 *
 * Supports explicit edge selection, angle-based selection, or all edges.
 */
struct SFCGAL_API EdgeSelector {
  /**
   * @brief Selection mode enumeration
   */
  enum class Mode {
    EXPLICIT, ///< Select specific edges by coordinates
    BY_ANGLE, ///< Select edges where dihedral angle exceeds threshold
    ALL       ///< Select all edges
  };

  Mode                        mode{Mode::ALL};      ///< Selection mode
  std::vector<EdgeIdentifier> edges;                ///< For EXPLICIT mode
  double angleThreshold{M_PI / 4}; ///< For BY_ANGLE mode (radians, default 45
                                   ///< degrees)

  /**
   * @brief Create selector for explicit edge list
   * @param edgeList List of edges to select
   * @return EdgeSelector configured for explicit selection
   */
  static auto
  explicit_(std::vector<EdgeIdentifier> edgeList) -> EdgeSelector
  {
    EdgeSelector sel;
    sel.mode  = Mode::EXPLICIT;
    sel.edges = std::move(edgeList);
    return sel;
  }

  /**
   * @brief Create selector for angle-based selection
   * @param threshold Minimum dihedral angle in radians to select an edge
   * @return EdgeSelector configured for angle-based selection
   */
  static auto
  byAngle(double threshold) -> EdgeSelector
  {
    EdgeSelector sel;
    sel.mode           = Mode::BY_ANGLE;
    sel.angleThreshold = threshold;
    return sel;
  }

  /**
   * @brief Create selector for all edges
   * @return EdgeSelector configured to select all edges
   */
  static auto
  all() -> EdgeSelector
  {
    EdgeSelector sel;
    sel.mode = Mode::ALL;
    return sel;
  }
};

/**
 * @brief Chamfer parameters for edge chamfering.
 *
 * Supports symmetric (single distance) and asymmetric (two distances) chamfers.
 */
struct SFCGAL_API ChamferParameters {
  /**
   * @brief Chamfer type enumeration
   */
  enum class Type {
    SYMMETRIC,  ///< Equal distance from both adjacent faces
    ASYMMETRIC  ///< Different distances from each adjacent face
  };

  Type   type{Type::SYMMETRIC}; ///< Chamfer type
  double distance1{0.0};        ///< Primary distance (used for both symmetric
                                ///< and asymmetric)
  double distance2{0.0};        ///< Secondary distance (only for asymmetric)

  /**
   * @brief Create symmetric chamfer parameters
   * @param distance Distance from edge to chamfer cut
   * @return ChamferParameters configured for symmetric chamfer
   */
  static auto
  symmetric(double distance) -> ChamferParameters
  {
    ChamferParameters params;
    params.type      = Type::SYMMETRIC;
    params.distance1 = distance;
    params.distance2 = distance;
    return params;
  }

  /**
   * @brief Create asymmetric chamfer parameters
   * @param dist1 Distance from first adjacent face
   * @param dist2 Distance from second adjacent face
   * @return ChamferParameters configured for asymmetric chamfer
   */
  static auto
  asymmetric(double dist1, double dist2) -> ChamferParameters
  {
    ChamferParameters params;
    params.type      = Type::ASYMMETRIC;
    params.distance1 = dist1;
    params.distance2 = dist2;
    return params;
  }
};

/**
 * @brief Vertex chamfer parameters.
 *
 * Controls how vertices are chamfered by specifying the distance along
 * each incident edge where the cut should be made.
 */
struct SFCGAL_API VertexChamferParameters {
  double distance{0.0}; ///< Distance from vertex along each incident edge

  /**
   * @brief Create vertex chamfer parameters
   * @param dist Distance from vertex
   * @return VertexChamferParameters
   */
  static auto
  create(double dist) -> VertexChamferParameters
  {
    VertexChamferParameters params;
    params.distance = dist;
    return params;
  }
};

/**
 * @brief 3D chamfer operation on edges and vertices of solids or polyhedral
 * surfaces.
 *
 * Uses a hybrid approach: attempts direct mesh modification first, then
 * falls back to boolean operations if direct modification fails.
 */
class SFCGAL_API Chamfer3D {
public:
  /**
   * @brief Construct a Chamfer3D operation
   * @param geometry Input geometry (must be Solid or PolyhedralSurface)
   * @throws std::invalid_argument if geometry is not Solid or PolyhedralSurface
   */
  explicit Chamfer3D(const Geometry &geometry);

  /**
   * @brief Apply edge chamfer to selected edges
   * @param selector Edge selection criteria
   * @param params Chamfer parameters
   * @return New geometry with chamfered edges
   * @throws std::runtime_error if chamfering fails
   */
  auto
  chamferEdges(const EdgeSelector      &selector,
               const ChamferParameters &params) const
      -> std::unique_ptr<Geometry>;

  /**
   * @brief Apply chamfer to all edges meeting at specified vertices
   * @param vertices List of vertex coordinates to chamfer
   * @param params Vertex chamfer parameters
   * @return New geometry with chamfered vertices
   * @throws std::runtime_error if chamfering fails
   */
  auto
  chamferVertices(const std::vector<Kernel::Point_3> &vertices,
                  const VertexChamferParameters      &params) const
      -> std::unique_ptr<Geometry>;

  /**
   * @brief Find all sharp edges based on dihedral angle
   * @param angleThreshold Minimum dihedral angle (radians) to consider "sharp"
   * @return List of edge identifiers for edges exceeding the threshold
   */
  auto
  findSharpEdges(double angleThreshold) const -> std::vector<EdgeIdentifier>;

  /**
   * @brief Compute dihedral angle at an edge
   * @param edge The edge to compute angle for
   * @return Dihedral angle in radians, or -1.0 if edge not found
   */
  auto
  computeDihedralAngle(const EdgeIdentifier &edge) const -> double;

private:
  const Geometry *_geometry;           ///< Input geometry (non-owning)
  bool            _isSolid{false};     ///< True if input is Solid
  bool            _isPolyhedralSurface{false}; ///< True if input is
                                               ///< PolyhedralSurface

  /**
   * @brief Convert input geometry to Surface_mesh
   * @return CGAL Surface_mesh representation
   */
  auto
  toSurfaceMesh() const -> Surface_mesh_3;

  /**
   * @brief Convert Surface_mesh back to SFCGAL geometry
   * @param mesh The mesh to convert
   * @return Appropriate SFCGAL geometry type
   */
  auto
  fromSurfaceMesh(const Surface_mesh_3 &mesh) const -> std::unique_ptr<Geometry>;

  /**
   * @brief Resolve edge selector to list of edges
   * @param selector Edge selection criteria
   * @return List of resolved edge identifiers
   */
  auto
  resolveEdges(const EdgeSelector &selector) const
      -> std::vector<EdgeIdentifier>;

  /**
   * @brief Find halfedge in mesh corresponding to edge identifier
   * @param mesh The mesh to search
   * @param edge The edge to find
   * @return Halfedge index, or null_halfedge() if not found
   */
  auto
  findHalfedge(const Surface_mesh_3 &mesh, const EdgeIdentifier &edge) const
      -> Surface_mesh_3::Halfedge_index;

  /**
   * @brief Extract edges from original geometry patches (not triangulated mesh)
   * @return List of edges from the original polygon patches
   */
  auto
  extractOriginalEdges() const -> std::vector<EdgeIdentifier>;

  /**
   * @brief Apply chamfer to single edge using direct mesh modification
   * @param mesh The mesh to modify
   * @param he Halfedge representing the edge
   * @param params Chamfer parameters
   * @return true if successful, false if modification failed
   */
  auto
  chamferEdgeDirect(Surface_mesh_3                   &mesh,
                    Surface_mesh_3::Halfedge_index    he,
                    const ChamferParameters          &params) const -> bool;

  /**
   * @brief Apply chamfer to single edge using boolean subtraction (fallback)
   * @param geom Current geometry state
   * @param edge Edge to chamfer
   * @param params Chamfer parameters
   * @return New geometry with edge chamfered
   */
  auto
  chamferEdgeBoolean(const Geometry          &geom,
                     const EdgeIdentifier    &edge,
                     const ChamferParameters &params) const
      -> std::unique_ptr<Geometry>;

  /**
   * @brief Create wedge solid for boolean chamfer operation
   * @param edge Edge to create wedge for
   * @param params Chamfer parameters
   * @param n1 Normal of first adjacent face
   * @param n2 Normal of second adjacent face
   * @return Solid representing material to remove
   */
  auto
  createChamferWedge(const EdgeIdentifier    &edge,
                     const ChamferParameters &params,
                     const Kernel::Vector_3  &n1,
                     const Kernel::Vector_3  &n2) const
      -> std::unique_ptr<Solid>;

  /**
   * @brief Apply chamfer to single vertex using direct mesh modification
   * @param mesh The mesh to modify
   * @param v Vertex to chamfer
   * @param params Vertex chamfer parameters
   * @return true if successful, false if modification failed
   */
  auto
  chamferVertexDirect(Surface_mesh_3                     &mesh,
                      Surface_mesh_3::Vertex_index        v,
                      const VertexChamferParameters      &params) const -> bool;
};

// ============================================================================
// Convenience Functions
// ============================================================================

/**
 * @brief Apply chamfer to specific edges
 *
 * @param geometry Input geometry (Solid or PolyhedralSurface)
 * @param edges List of edges to chamfer
 * @param params Chamfer parameters
 * @return Chamfered geometry
 * @throws std::invalid_argument if geometry type is not supported
 */
SFCGAL_API auto
chamfer3D(const Geometry                    &geometry,
          const std::vector<EdgeIdentifier> &edges,
          const ChamferParameters           &params) -> std::unique_ptr<Geometry>;

/**
 * @brief Apply vertex chamfer to specific vertices
 *
 * @param geometry Input geometry (Solid or PolyhedralSurface)
 * @param vertices List of vertex coordinates to chamfer
 * @param params Vertex chamfer parameters
 * @return Chamfered geometry
 * @throws std::invalid_argument if geometry type is not supported
 */
SFCGAL_API auto
chamfer3DVertices(const Geometry                     &geometry,
                  const std::vector<Kernel::Point_3> &vertices,
                  const VertexChamferParameters      &params)
    -> std::unique_ptr<Geometry>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_CHAMFER3D_ALGORITHM
