// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TRIANGULATE_CONSTRAINTDELAUNAYTRIANGULATION_H_
#define SFCGAL_TRIANGULATE_CONSTRAINTDELAUNAYTRIANGULATION_H_

#include <optional>

#include "SFCGAL/Coordinate.h"
#include "SFCGAL/config.h"

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

namespace SFCGAL {
class TriangulatedSurface;
}

namespace SFCGAL {
namespace triangulate {

/**
 * @brief 2DZ constraint Delaunay triangulation
 */
class SFCGAL_API ConstraintDelaunayTriangulation {
public:
  /**
   * @brief vertex info in triangulation
   */
  struct VertexInfo {
    VertexInfo() : original() {}
    Coordinate original; ///< Original coordinate of the vertex
  };

  /**
   * face information (depth)
   */
  struct FaceInfo {
    FaceInfo() : nestingLevel(-1) {}
    int nestingLevel; ///< Nesting level for domain determination

    /**
     * @brief Check if face is in the domain
     * @return True if face is in the domain (odd nesting level)
     */
    auto
    in_domain() -> bool
    {
      return nestingLevel % 2 == 1;
    }
  };

  typedef CGAL::Triangulation_vertex_base_with_info_2<VertexInfo, Kernel>
      Triangulation_vertex_base; ///< CGAL vertex base type with info
  typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, Kernel>
      Triangulation_face_base; ///< CGAL face base type with info
  typedef CGAL::Constrained_triangulation_face_base_2<Kernel,
                                                      Triangulation_face_base>
      Constrained_triangulation_face_base; ///< CGAL constrained triangulation
                                           ///< face base
  typedef CGAL::Triangulation_data_structure_2<
      Triangulation_vertex_base, Constrained_triangulation_face_base>
      Triangulation_data_structure; ///< CGAL triangulation data structure type

  typedef CGAL::Constrained_Delaunay_triangulation_2<
      Kernel, Triangulation_data_structure, CGAL::Exact_predicates_tag>
      CDT; ///< CGAL constrained Delaunay triangulation type

  typedef CDT::Vertex_handle
      Vertex_handle;                    ///< Handle to vertex in triangulation
  using Face_handle = CDT::Face_handle; ///< Handle to face in triangulation
  typedef CDT::All_faces_iterator
      All_faces_iterator; ///< Iterator over all faces
  typedef CDT::Finite_faces_iterator
      Finite_faces_iterator; ///< Iterator over finite faces

  /**
   * @brief default constructor
   */
  ConstraintDelaunayTriangulation();

  /**
   * @brief Add a vertex to the triangulation
   * @param position The coordinate position of the vertex
   * @return Handle to the added vertex
   */
  auto
  addVertex(const Coordinate &position) -> Vertex_handle;
  /**
   * @brief Add a constraint edge to the triangulation
   * @param source The source vertex handle
   * @param target The target vertex handle
   */
  auto
  addConstraint(Vertex_handle source, Vertex_handle target) -> void;

  /**
   * @brief clear the triangulation
   */
  auto
  clear() -> void;

  /**
   * @brief Returns the number of finite vertices
   * @return The number of finite vertices in the triangulation
   */
  auto
  numVertices() const -> size_t;
  /**
   * @brief Returns the number of finite faces
   * @return The number of finite triangles in the triangulation
   */
  auto
  numTriangles() const -> size_t;

  /**
   * @brief test if a projection plane is defined
   * @return True if a projection plane has been set
   */
  inline auto
  hasProjectionPlane() const -> bool
  {
    return _projectionPlane.has_value();
  }
  /**
   * @brief define projection plane
   * @param projectionPlane The 3D plane to use for projection
   */
  auto
  setProjectionPlane(const Kernel::Plane_3 &projectionPlane) -> void;
  /**
   * @brief get the projection plane (OXY if not defined)
   * @return The projection plane (defaults to OXY plane if not set)
   */
  auto
  projectionPlane() const -> Kernel::Plane_3;

  /**
   * @brief test if the vertex is infinite
   * @param vertex The vertex handle to test
   * @return True if the vertex is infinite
   */
  inline auto
  isInfinite(Vertex_handle vertex) const -> bool
  {
    return _cdt.is_infinite(vertex);
  }
  /**
   * @brief test if the face has infinite vertex
   * @param face The face handle to test
   * @return True if the face has an infinite vertex
   */
  inline auto
  isInfinite(Face_handle face) const -> bool
  {
    return _cdt.is_infinite(face);
  }

  /**
   * @brief Append Triangles to a TriangulatedSurface
   * @param triangulatedSurface The surface to append triangles to
   * @param filterExteriorParts If true, exclude exterior triangles
   */
  auto
  getTriangles(TriangulatedSurface &triangulatedSurface,
               bool                 filterExteriorParts = false) const -> void;
  /**
   * get the resulting TriangulatedSurface
   * @return Unique pointer to the triangulated surface
   */
  auto
  getTriangulatedSurface() const -> std::unique_ptr<TriangulatedSurface>;

  /**
   * @brief get finite face iterator
   * @return Iterator to beginning of finite faces
   */
  inline auto
  finite_faces_begin() const -> Finite_faces_iterator
  {
    return _cdt.finite_faces_begin();
  }
  /**
   * @brief get finite face iterator
   * @return Iterator to end of finite faces
   */
  inline auto
  finite_faces_end() const -> Finite_faces_iterator
  {
    return _cdt.finite_faces_end();
  }

  /**
   * @brief get all face iterator
   * @return Iterator to beginning of all faces
   */
  inline auto
  all_faces_begin() const -> All_faces_iterator
  {
    return _cdt.all_faces_begin();
  }
  /**
   * @brief get all face iterator
   * @return Iterator to end of all faces
   */
  inline auto
  all_faces_end() const -> All_faces_iterator
  {
    return _cdt.all_faces_end();
  }

  /**
   * @brief fill nesting_level info in face info
   */
  auto
  markDomains() -> void;

  /**
   * @brief [advanced]get the CGAL object
   * @return Reference to the underlying CGAL triangulation
   */
  inline auto
  cdt() -> CDT &
  {
    return _cdt;
  }
  /**
   * @brief [advanced]get the CGAL object
   * @return Const reference to the underlying CGAL triangulation
   */
  inline auto
  cdt() const -> const CDT &
  {
    return _cdt;
  }

private:
  /**
   * @brief wrapped triangulation
   */
  CDT _cdt;
  /**
   * @brief plan in which the triangulation is done
   */
  std::optional<Kernel::Plane_3> _projectionPlane;
};

} // namespace triangulate
} // namespace SFCGAL

#endif
