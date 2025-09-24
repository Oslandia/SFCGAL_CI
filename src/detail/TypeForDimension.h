// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_DETAIL_TYPE_FOR_DIMENSION_H
#define SFCGAL_DETAIL_TYPE_FOR_DIMENSION_H

#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Triangle_3.h>

namespace SFCGAL {
namespace detail {

/*
/// Type traits for CGAL types.
/// CGAL types cannot be directly parametrized by their dimension.
/// For instance, there are both a Triangle_2<K> and a Triangle_3<K> type
TypeForKernel<K, 2>::Point is equivalent to CGAL::Point_2<K> for instance
TypeForDimension<2>::Bbox is equivalent to CGAL::Bbox_2
*/

struct SFCGAL_API NoVolume{};

/// Generic traits, default dimension is 2
template <int Dim>
struct TypeForDimension {
  typedef CGAL::Bbox_2                       Bbox; ///< Bounding box type for 2D
  typedef Kernel::Point_2                    Point;    ///< Point type for 2D
  typedef Kernel::Segment_2                  Segment;  ///< Segment type for 2D
  typedef Kernel::Triangle_2                 Triangle; ///< Triangle type for 2D
  typedef CGAL::Polygon_with_holes_2<Kernel> Surface;  ///< Surface type for 2D
  typedef NoVolume Volume; ///< Volume type for 2D (none)
};

/// Extended polyhedron_3 type with a boolean marker on halfedges
/// This is used internally for polyhedra boolean operations
template <class Refs>
struct Halfedge_with_mark : public CGAL::HalfedgeDS_halfedge_base<Refs> {
  Halfedge_with_mark() : CGAL::HalfedgeDS_halfedge_base<Refs>(), mark(false) {}

  bool mark; ///< A boundary marker for faces with different status
  /**
   * @brief Set the boundary marker to true
   *
   * This will be called by Intersection_of_Polyhedra_3
   */
  void
  set_mark()
  {
    mark = true;
  }
};

/// An items type using my halfedge.
struct SFCGAL_API Items_with_mark_on_hedge : public CGAL::Polyhedron_items_3 {
  template <class Refs, class Traits>
  /**
   * @brief Wrapper for halfedges with mark functionality
   */
  struct Halfedge_wrapper {
    typedef Halfedge_with_mark<Refs>
        Halfedge; ///< Halfedge type with marking capability
  };
};

/// @brief CGAL 3D polyhedron with marked half-edges
typedef CGAL::Polyhedron_3<Kernel, Items_with_mark_on_hedge> MarkedPolyhedron;

/// Specialization for dimension = 3
template <>
struct TypeForDimension<3> {
  typedef CGAL::Bbox_3       Bbox;     ///< Bounding box type for 3D
  typedef Kernel::Point_3    Point;    ///< Point type for 3D
  typedef Kernel::Segment_3  Segment;  ///< Segment type for 3D
  typedef Kernel::Triangle_3 Triangle; ///< Triangle type for 3D
  typedef Kernel::Triangle_3 Surface;  ///< Surface type for 3D
  typedef MarkedPolyhedron   Volume; ///< Volume type for 3D (marked polyhedron)
};

/// Another way of looking at TypeForDimension<Dim>::Point
template <int Dim>
struct Point_d {
  typedef typename TypeForDimension<Dim>::Point
      Type; ///< Point type for given dimension
};
/// Another way of looking at TypeForDimension<Dim>::Segment
template <int Dim>
struct Segment_d {
  typedef typename TypeForDimension<Dim>::Segment
      Type; ///< Segment type for given dimension
};
/// Another way of looking at TypeForDimension<Dim>::Surface
template <int Dim>
struct Surface_d {
  typedef typename TypeForDimension<Dim>::Surface
      Type; ///< Surface type for given dimension
};
/// Another way of looking at TypeForDimension<Dim>::Volume
template <int Dim>
struct Volume_d {
  typedef typename TypeForDimension<Dim>::Volume
      Type; ///< Volume type for given dimension
};

/// Create a distinct type for each dimension
template <int N>
struct dim_t {
  enum { v = N }; ///< Dimension value
};

/// Get a primitive dimension (0: point, 1: line, 2: surface, 3: volume) from a
/// type
template <class T>
struct PrimitiveDimension {
  static const int value = 0; ///< Default dimension value for points
};

/**
 * @brief Specialization for 2D segments
 */
template <>
struct PrimitiveDimension<TypeForDimension<2>::Segment> {
  static const int value = 1; ///< 2D segments have dimension 1
};
/**
 * @brief Specialization for 3D segments
 */
template <>
struct PrimitiveDimension<TypeForDimension<3>::Segment> {
  static const int value = 1; ///< 3D segments have dimension 1
};
/**
 * @brief Specialization for 2D surfaces
 */
template <>
struct PrimitiveDimension<TypeForDimension<2>::Surface> {
  static const int value = 2; ///< 2D surfaces have dimension 2
};
/**
 * @brief Specialization for 3D surfaces
 */
template <>
struct PrimitiveDimension<TypeForDimension<3>::Surface> {
  static const int value = 2; ///< 3D surfaces have dimension 2
};
/**
 * @brief Specialization for 2D volumes (should not exist but for completeness)
 */
template <>
struct PrimitiveDimension<TypeForDimension<2>::Volume> {
  static const int value = 3; ///< 2D volumes have dimension 3
};
/**
 * @brief Specialization for 3D volumes
 */
template <>
struct PrimitiveDimension<TypeForDimension<3>::Volume> {
  static const int value = 3; ///< 3D volumes have dimension 3
};

/// Tests if a primitive type has a larger dimension than another one
template <class X, class Y>
struct IsPrimitiveLarger {
  static const bool value =
      PrimitiveDimension<X>::value >
      PrimitiveDimension<Y>::value; ///< True if X has larger dimension than Y
};

} // namespace detail
} // namespace SFCGAL

#endif
