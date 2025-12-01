// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_DETAIL_GEOMETRY_SET_H_
#define SFCGAL_DETAIL_GEOMETRY_SET_H_

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/variant.hpp>

#include <cstdint>

#include "SFCGAL/config.h"

#include "SFCGAL/Kernel.h"
#include "SFCGAL/detail/TypeForDimension.h"

#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Box_intersection_d/Box_with_handle_d.h>

/// comparison operator on 2D segments, for use in a std::set
auto
operator<(const CGAL::Segment_2<SFCGAL::Kernel> &segmentA,
          const CGAL::Segment_2<SFCGAL::Kernel> &segmentB) -> bool;

/// comparison operator on 3D segments, for use in a std::set
auto
operator<(const CGAL::Segment_3<SFCGAL::Kernel> &segmentA,
          const CGAL::Segment_3<SFCGAL::Kernel> &segmentB) -> bool;

namespace SFCGAL {
class Geometry;
}

namespace SFCGAL::detail {

/// Primitive type enumeration. Note that the value is the dimension !
enum PrimitiveType : std::uint8_t {
  PrimitivePoint   = 0,
  PrimitiveSegment = 1,
  PrimitiveSurface = 2,
  PrimitiveVolume  = 3
};

/// Primitive handle. Holds a pointer to a primitive, through the 'handle'
/// member
template <int Dim>
struct PrimitiveHandle {
  //
  // We use boost::variant here for convenience, whereas it is needed
  using Type = boost::variant<const typename Point_d<Dim>::Type *,
                              const typename Segment_d<Dim>::Type *,
                              const typename Surface_d<Dim>::Type *,
                              const typename Volume_d<Dim>::Type *>;
  Type handle; ///< The stored primitive handle

  /**
   * @brief Constructor from primitive pointer
   * @param primitivePointer Pointer to the primitive to store
   */
  template <class T>
  PrimitiveHandle(const T *primitivePointer) : handle(primitivePointer)
  {
  }

  /**
   * @brief Cast handle to specific type
   * @return Pointer to primitive of specified type
   */
  template <class T>
  [[nodiscard]] auto
  as() const -> const T *
  {
    return boost::get<const T *>(handle);
  }
};

/// PrimitiveBox. Type used for CGAL::Box_intersection_d
template <int Dim>
struct PrimitiveBox {
  using Type =
      CGAL::Box_intersection_d::Box_with_handle_d<double, Dim,
                                                  PrimitiveHandle<Dim> *>;
};

/// BoxCollection for use with CGAL::Box_intersection_d
template <int Dim>
struct BoxCollection {
  using Type = std::vector<typename PrimitiveBox<Dim>::Type>;
};

/// HandleCollection. Used to store PrimitiveHandle
template <int Dim>
struct HandleCollection {
  using Type = std::list<PrimitiveHandle<Dim>>;
};

/// Flags available for each type of Geometry type.
/// Primitives can be 'flagged' in order to speed up recomposition
enum ElementFlag : std::uint8_t {
  // the polyhedron is planar => build a triangle or a polygon
  FLAG_IS_PLANAR = 1
};

/// CollectionElement, a Primitive with flags
/// Primitive : Point_d, Segment_d, Surface_d, Volume_d
template <class Primitive>
class CollectionElement {
public:
  /**
   * @brief Get element flags
   * @return The flags associated with this element
   */
  [[nodiscard]] auto
  flags() const -> int
  {
    return _flags;
  }
  /**
   * @brief Set element flags
   * @param flags The flags to set for this element
   */
  auto
  setFlags(int flags) -> void
  {
    _flags = flags;
  }

  /**
   * @brief Get mutable reference to primitive
   * @return Reference to the underlying primitive
   */
  [[nodiscard]] auto
  primitive() -> Primitive &
  {
    return _primitive;
  }
  /**
   * @brief Get const reference to primitive
   * @return Const reference to the underlying primitive
   */
  [[nodiscard]] auto
  primitive() const -> const Primitive &
  {
    return _primitive;
  }

  /**
   * @brief Default constructor
   */
  CollectionElement() : _flags(0) {}

  /**
   * @brief Constructor from primitive
   * @param primitive The primitive to wrap
   */
  CollectionElement(const Primitive &primitive)
      : _primitive(primitive), _flags(0)
  {
  }

  /**
   * @brief Constructor from primitive with flags
   * @param primitive The primitive to wrap
   * @param flags The initial flags
   */
  CollectionElement(const Primitive &primitive, int flags)
      : _primitive(primitive), _flags(flags)
  {
  }

  /**
   * @brief Copy constructor
   * @param other The element to copy from
   */
  CollectionElement(const CollectionElement &other)
      : _primitive(other._primitive), _flags(other._flags)
  {
  }

  /**
   * @brief Copy assignment operator
   * @param other The element to copy from
   * @return Reference to this
   */
  auto
  operator=(const CollectionElement &other) -> CollectionElement & = default;

  /**
   * @brief Less-than comparison operator
   * @param other The element to compare with
   * @return True if this element is less than other
   */
  auto
  operator<(const CollectionElement &other) const -> bool
  {
    return _primitive < other._primitive;
  }

private:
  Primitive _primitive;
  int       _flags;
};

/**
 * Stream output operator for a CollectionElement.
 *
 * @tparam Primitive The type of the underlying primitive stored in the
 * CollectionElement.
 * @param ostr The output stream to write to.
 * @param collection The CollectionElement to display.
 * @return The output stream after writing the CollectionElement.
 */
template <class Primitive>
auto
operator<<(std::ostream &ostr, const CollectionElement<Primitive> &collection)
    -> std::ostream &
{
  ostr << collection.primitive() << " flags: " << collection.flags();
  return ostr;
}

/// A GeometrySet represents a set of CGAL primitives.
/// Primitive are either of dimension 0 (points),
/// dimension 1 (segments), dimension 2 (surfaces, a.k.a. polygon or triangles)
/// or dimension 3 (polyhedron)
template <int Dim>
class GeometrySet {
public:
  // Points are stored in an ordered set
  using PointCollection =
      std::set<CollectionElement<typename Point_d<Dim>::Type>>;
  // Segments are stored in an ordered set
  using SegmentCollection =
      std::set<CollectionElement<typename Segment_d<Dim>::Type>>;
  using SurfaceCollection =
      std::list<CollectionElement<typename Surface_d<Dim>::Type>>;
  using VolumeCollection =
      std::list<CollectionElement<typename Volume_d<Dim>::Type>>;

  GeometrySet();

  /**
   * Construct a GeometrySet from a SFCGAL::Geometry
   * @param g The geometry to convert to GeometrySet
   */
  GeometrySet(const Geometry &g);

  /**
   * Construct a GeometrySet from a Point
   * @param g The point to add to the set
   * @param flags Optional flags for the point
   */
  GeometrySet(const typename TypeForDimension<Dim>::Point &g, int flags = 0);

  /**
   * Construct a GeometrySet from a Segment
   * @param g The segment to add to the set
   * @param flags Optional flags for the segment
   */
  GeometrySet(const typename TypeForDimension<Dim>::Segment &g, int flags = 0);

  /**
   * Construct a GeometrySet from a Surface
   * @param g The surface to add to the set
   * @param flags Optional flags for the surface
   */
  GeometrySet(const typename TypeForDimension<Dim>::Surface &g, int flags = 0);

  /**
   * Construct a GeometrySet from a Volume
   * @param g The volume to add to the set
   * @param flags Optional flags for the volume
   */
  GeometrySet(const typename TypeForDimension<Dim>::Volume &g, int flags = 0);

  /**
   * Add primitives from another set
   * @param g The geometry set to merge from
   */
  auto
  merge(const GeometrySet<Dim> &g) -> void;

  /**
   * Add a geometry by decomposing it into CGAL primitives
   * @param g The geometry to decompose and add
   */
  void
  addGeometry(const Geometry &g);

  /**
   * add a primitive from a PrimitiveHandle  to the set
   * @param handle The primitive handle to add
   */
  void
  addPrimitive(const PrimitiveHandle<Dim> &handle);

  /**
   * add a primitive from a CGAL::Object to the set
   * pointsAsRing : if set to true, build a polygon if obj is a vector of points
   * @param obj The CGAL object to add as primitive
   * @param pointsAsRing If true, build polygon from point vector
   */
  void
  addPrimitive(const CGAL::Object &obj, bool pointsAsRing = false);

  /**
   * add a point to the set
   * @param g The point to add
   * @param flags Optional flags for the point
   */
  void
  addPrimitive(const typename TypeForDimension<Dim>::Point &g, int flags = 0);
  /**
   * @brief Add multiple points from iterator range
   * @param ibegin Iterator to first point
   * @param iend Iterator to end of points
   */
  template <class IT>
  void
  addPoints(IT ibegin, IT iend)
  {
    std::copy(ibegin, iend, std::inserter(_points, _points.end()));
  }

  /**
   * collect all points of b and add them to the point list
   * @param b The primitive handle to extract points from
   */
  void
  collectPoints(const PrimitiveHandle<Dim> &b);

  /**
   * add a segment to the set
   * @param g The segment to add
   * @param flags Optional flags for the segment
   */
  void
  addPrimitive(const typename TypeForDimension<Dim>::Segment &g, int flags = 0);
  /**
   * @brief Add multiple segments from iterator range
   * @param ibegin Iterator to first segment
   * @param iend Iterator to end of segments
   */
  template <class IT>
  void
  addSegments(IT ibegin, IT iend)
  {
    std::copy(ibegin, iend, std::inserter(_segments, _segments.end()));
  }

  /**
   * add a surface to the set
   * @param g The surface to add
   * @param flags Optional flags for the surface
   */
  void
  addPrimitive(const typename TypeForDimension<Dim>::Surface &g, int flags = 0);
  /**
   * @brief Add multiple surfaces from iterator range
   * @param ibegin Iterator to first surface
   * @param iend Iterator to end of surfaces
   */
  template <class IT>
  void
  addSurfaces(IT ibegin, IT iend)
  {
    std::copy(ibegin, iend, std::back_inserter(_surfaces));
  }

  /**
   * add a volume to the set
   * @param g The volume to add
   * @param flags Optional flags for the volume
   */
  void
  addPrimitive(const typename TypeForDimension<Dim>::Volume &g, int flags = 0);
  /**
   * @brief Add multiple volumes from iterator range
   * @param ibegin Iterator to first volume
   * @param iend Iterator to end of volumes
   */
  template <class IT>
  void
  addVolumes(IT ibegin, IT iend)
  {
    std::copy(ibegin, iend, std::back_inserter(_volumes));
  }

  /**
   * Get the maximum geometry dimension of the set
   * -1 : empty
   * 0 : there are points
   * 1 : there are segments
   * 2 : there are surfaces
   * 3 : there are volumes
   * @return The maximum dimension of geometries in the set
   */
  [[nodiscard]] auto
  dimension() const -> int;

  /**
   * Add the boundary (segments) of a surface
   * @param surface The surface whose boundary to add
   */
  void
  addBoundary(const typename TypeForDimension<Dim>::Surface &surface);

  /**
   * Add the boundary (surfaces) of a volume
   * @param volume The volume whose boundary to add
   */
  void
  addBoundary(const typename TypeForDimension<Dim>::Volume &volume);

  /**
   * Compute all bounding boxes and handles of the set
   * @param handles Output collection to store primitive handles
   * @param boxes Output collection to store bounding boxes
   */
  void
  computeBoundingBoxes(typename HandleCollection<Dim>::Type &handles,
                       typename BoxCollection<Dim>::Type    &boxes) const;

  /**
   * @brief Get mutable reference to point collection
   * @return Reference to the point collection
   */
  [[nodiscard]] auto
  points() -> PointCollection &
  {
    return _points;
  }
  /**
   * @brief Get const reference to point collection
   * @return Const reference to the point collection
   */
  [[nodiscard]] auto
  points() const -> const PointCollection &
  {
    return _points;
  }

  /**
   * @brief Get mutable reference to segment collection
   * @return Reference to the segment collection
   */
  [[nodiscard]] auto
  segments() -> SegmentCollection &
  {
    return _segments;
  }
  /**
   * @brief Get const reference to segment collection
   * @return Const reference to the segment collection
   */
  [[nodiscard]] auto
  segments() const -> const SegmentCollection &
  {
    return _segments;
  }

  /**
   * @brief Get mutable reference to surface collection
   * @return Reference to the surface collection
   */
  [[nodiscard]] auto
  surfaces() -> SurfaceCollection &
  {
    return _surfaces;
  }
  /**
   * @brief Get const reference to surface collection
   * @return Const reference to the surface collection
   */
  [[nodiscard]] auto
  surfaces() const -> const SurfaceCollection &
  {
    return _surfaces;
  }

  /**
   * @brief Get mutable reference to volume collection
   * @return Reference to the volume collection
   */
  [[nodiscard]] auto
  volumes() -> VolumeCollection &
  {
    return _volumes;
  }
  /**
   * @brief Get const reference to volume collection
   * @return Const reference to the volume collection
   */
  [[nodiscard]] auto
  volumes() const -> const VolumeCollection &
  {
    return _volumes;
  }

  /**
   * Returns true if the set holds points
   * @return True if the set contains points
   */
  [[nodiscard]] auto
  hasPoints() const -> bool;
  /**
   * Returns true if the set holds segments
   * @return True if the set contains segments
   */
  [[nodiscard]] auto
  hasSegments() const -> bool;
  /**
   * Returns true if the set holds surfaces
   * @return True if the set contains surfaces
   */
  [[nodiscard]] auto
  hasSurfaces() const -> bool;
  /**
   * Returns true if the set holds volumes
   * @return True if the set contains volumes
   */
  [[nodiscard]] auto
  hasVolumes() const -> bool;

  /**
   * convert the set to a SFCGAL::Geometry
   * @return Unique pointer to the recomposed geometry
   */
  [[nodiscard]] auto
  recompose() const -> std::unique_ptr<Geometry>;

  /**
   * Filter (remove) primitives that are already covered by others
   * @param output The output geometry set with filtered primitives
   */
  void
  filterCovered(GeometrySet<Dim> &output) const;

private:
  ///
  /// Given an input SFCGAL::Geometry, decompose it into CGAL primitives
  void
  _decompose(const Geometry &g);

  PointCollection   _points;
  SegmentCollection _segments;
  SurfaceCollection _surfaces;
  VolumeCollection  _volumes;
};

/**
 * @brief Display operator for 2D GeometrySet.
 * @param ostr The output stream to write to.
 * @param geomSet The 2D GeometrySet to display.
 * @return The output stream after writing.
 */
SFCGAL_API auto
operator<<(std::ostream &ostr, const GeometrySet<2> &geomSet) -> std::ostream &;

/**
 * @brief Display operator for 3D GeometrySet.
 * @param ostr The output stream to write to.
 * @param geomSet The 3D GeometrySet to display.
 * @return The output stream after writing.
 */
SFCGAL_API auto
operator<<(std::ostream &ostr, const GeometrySet<3> &geomSet) -> std::ostream &;

/**
 * Compute the bounding box of a volume in 2D.
 *
 * This overload is for the `NoVolume` type and will never be called in
 * practice.
 *
 * @return An empty 2D bounding box.
 */
inline auto
compute_solid_bbox(const NoVolume & /*unused*/, dim_t<2> /*unused*/)
    -> CGAL::Bbox_2
{
  return {};
}

/**
 * Compute the bounding box of a 3D volume.
 *
 * @param volume The volume whose bounding box is to be computed.
 * @return The bounding box of the volume as a CGAL::Bbox_3.
 *
 */
inline auto
compute_solid_bbox(const TypeForDimension<3>::Volume &volume,
                   dim_t<3> /*unused*/) -> CGAL::Bbox_3
{
  BOOST_ASSERT(volume.size_of_vertices());
  MarkedPolyhedron::Point_const_iterator pit = volume.points_begin();
  CGAL::Bbox_3                           ret(pit->bbox());
  ++pit;

  for (; pit != volume.points_end(); ++pit) {
    ret = ret + pit->bbox();
  }

  return ret;
}
} // namespace SFCGAL::detail

#endif
