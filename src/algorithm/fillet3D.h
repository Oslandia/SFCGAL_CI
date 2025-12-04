// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_FILLET3D_ALGORITHM
#define SFCGAL_FILLET3D_ALGORITHM

#include "SFCGAL/Geometry.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/export.h"
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <optional>

namespace SFCGAL::algorithm {

/**
 * @brief Parameters for fillet operation
 */
struct SFCGAL_API FilletParameters {
  double radius{0.0};   ///< Fillet radius
  int    segments{16};  ///< Number of segments for arc approximation

  /**
   * @brief Create constant fillet parameters
   * @param r Fillet radius
   * @param segs Number of segments (default: 16)
   * @return FilletParameters with specified values
   */
  static auto
  constant(double r, int segs = 16) -> FilletParameters
  {
    return FilletParameters{r, segs};
  }
};

/**
 * @brief Edge identifier using two endpoints
 */
struct SFCGAL_API EdgeIdentifier {
  Kernel::Point_3 start;
  Kernel::Point_3 end;

  EdgeIdentifier(const Kernel::Point_3 &s, const Kernel::Point_3 &e)
      : start(s), end(e)
  {
  }

  auto
  operator==(const EdgeIdentifier &other) const -> bool
  {
    return (start == other.start && end == other.end) ||
           (start == other.end && end == other.start);
  }
};

/**
 * @brief Edge selector for specifying which edges to fillet
 */
class SFCGAL_API EdgeSelector {
public:
  enum class Type { ALL_CONVEX, EXPLICIT, ANGLE_THRESHOLD };

  /**
   * @brief Select all convex edges
   */
  static auto
  allConvex() -> EdgeSelector
  {
    EdgeSelector sel;
    sel._type = Type::ALL_CONVEX;
    return sel;
  }

  /**
   * @brief Select explicit edges by their identifiers
   */
  static auto
  explicit_(const std::vector<EdgeIdentifier> &edges) -> EdgeSelector
  {
    EdgeSelector sel;
    sel._type  = Type::EXPLICIT;
    sel._edges = edges;
    return sel;
  }

  /**
   * @brief Select edges with dihedral angle above threshold
   */
  static auto
  byAngleThreshold(double angleThresholdDegrees) -> EdgeSelector
  {
    EdgeSelector sel;
    sel._type           = Type::ANGLE_THRESHOLD;
    sel._angleThreshold = angleThresholdDegrees;
    return sel;
  }

  [[nodiscard]] auto
  type() const -> Type
  {
    return _type;
  }
  [[nodiscard]] auto
  edges() const -> const std::vector<EdgeIdentifier> &
  {
    return _edges;
  }
  [[nodiscard]] auto
  angleThreshold() const -> double
  {
    return _angleThreshold;
  }

private:
  Type                        _type{Type::ALL_CONVEX};
  std::vector<EdgeIdentifier> _edges;
  double                      _angleThreshold{90.0};
};

/**
 * @brief Computes 3D fillet (rounded edges) on solids
 *
 * The fillet operation creates rounded edges by constructing wedge-shaped
 * volumes along selected edges and subtracting them from the input geometry
 * using NEF polyhedron boolean operations.
 *
 * @note Known limitation: Edges at reflex (concave) corners are automatically
 * skipped to ensure valid output geometry. This means inner corners of L-shaped,
 * U-shaped, or similar geometries will remain sharp. This is a trade-off to
 * guarantee valid, non-self-intersecting output.
 */
class SFCGAL_API Fillet3D {
public:
  using Nef_polyhedron = CGAL::Nef_polyhedron_3<Kernel>;

  /**
   * @brief Constructs a Fillet3D object
   * @param inputGeometry The input geometry (must be a Solid or
   * PolyhedralSurface)
   * @throws std::invalid_argument if the input geometry is not suitable
   */
  explicit Fillet3D(const Geometry &inputGeometry);

  /**
   * @brief Apply fillet to selected edges
   * @param selector Edge selection criteria
   * @param params Fillet parameters
   * @return Filleted geometry
   */
  [[nodiscard]] auto
  filletEdges(const EdgeSelector     &selector,
              const FilletParameters &params) const -> std::unique_ptr<Geometry>;

  /**
   * @brief Apply fillet to all convex edges
   * @param radius Fillet radius
   * @param segments Number of arc segments (default: 16)
   * @return Filleted geometry
   */
  [[nodiscard]] auto
  fillet(double radius, int segments = 16) const -> std::unique_ptr<Geometry>;

private:
  const Geometry *_geometry;
  bool            _inputWasSolid;

  // Internal structures for edge processing
  struct EdgeInfo {
    Kernel::Point_3  start;
    Kernel::Point_3  end;
    Kernel::Vector_3 normal1;  // Normal of first adjacent face
    Kernel::Vector_3 normal2;  // Normal of second adjacent face
    double           dihedralAngle;
    bool             isConvex;
  };

  // Edge selection
  [[nodiscard]] auto
  selectEdges(const EdgeSelector &selector) const -> std::vector<EdgeInfo>;

  [[nodiscard]] auto
  findAllEdgesWithNormals() const -> std::vector<EdgeInfo>;

  [[nodiscard]] auto
  computeDihedralAngle(const Kernel::Vector_3 &n1,
                       const Kernel::Vector_3 &n2) const -> double;

  [[nodiscard]] auto
  isEdgeConvex(const Kernel::Vector_3 &n1, const Kernel::Vector_3 &n2,
               const Kernel::Vector_3 &edgeDir) const -> bool;

  // Wedge creation
  [[nodiscard]] auto
  createFilletWedgeWithPlanes(const EdgeInfo                       &edge,
                              const FilletParameters               &params,
                              const std::optional<Kernel::Plane_3> &startPlane,
                              const std::optional<Kernel::Plane_3> &endPlane) const
      -> CGAL::Polyhedron_3<Kernel>;

  // Geometry conversion
  [[nodiscard]] auto
  toNefPolyhedron() const -> Nef_polyhedron;

  [[nodiscard]] auto
  fromNefPolyhedron(const Nef_polyhedron &nef) const
      -> std::unique_ptr<Geometry>;

  // Vector utilities
  [[nodiscard]] static auto
  normalizeVector(const Kernel::Vector_3 &v) -> Kernel::Vector_3;

  [[nodiscard]] static auto
  rotateVector(const Kernel::Vector_3 &v, const Kernel::Vector_3 &axis,
               double angle) -> Kernel::Vector_3;
};

/**
 * @brief Convenience function to compute 3D fillet
 * @param geometry Input geometry (Solid or PolyhedralSurface)
 * @param radius Fillet radius
 * @param segments Number of arc segments (default: 16)
 * @return Filleted geometry
 *
 * @note Edges at reflex corners are skipped to ensure valid output.
 */
SFCGAL_API auto
fillet3D(const Geometry &geometry, double radius,
         int segments = 16) -> std::unique_ptr<Geometry>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_FILLET3D_ALGORITHM
