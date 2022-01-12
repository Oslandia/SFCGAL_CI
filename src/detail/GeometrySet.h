// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _SFCGAL_DETAIL_GEOMETRY_SET_H_
#define _SFCGAL_DETAIL_GEOMETRY_SET_H_

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/variant.hpp>

#include <SFCGAL/config.h>

#include <SFCGAL/Kernel.h>
#include <SFCGAL/detail/TypeForDimension.h>

#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Box_intersection_d/Box_with_handle_d.h>

// comparison operator on segments, for use in a std::set
bool
operator<(const CGAL::Segment_2<SFCGAL::Kernel> &sega,
          const CGAL::Segment_2<SFCGAL::Kernel> &segb);
bool
operator<(const CGAL::Segment_3<SFCGAL::Kernel> &sega,
          const CGAL::Segment_3<SFCGAL::Kernel> &segb);

namespace SFCGAL {
class Geometry;
namespace detail {

///
/// Primitive type enumeration. Note that the value is the dimension !
enum PrimitiveType {
  PrimitivePoint   = 0,
  PrimitiveSegment = 1,
  PrimitiveSurface = 2,
  PrimitiveVolume  = 3
};

///
/// Primitive handle. Holds a pointer to a primitive, through the 'handle'
/// member
template <int Dim>
struct PrimitiveHandle {
  //
  // We use boost::variant here for convenience, whereas it is needed
  typedef boost::variant<const typename Point_d<Dim>::Type *,
                         const typename Segment_d<Dim>::Type *,
                         const typename Surface_d<Dim>::Type *,
                         const typename Volume_d<Dim>::Type *>
       Type;
  Type handle;

  template <class T>
  PrimitiveHandle(const T *p) : handle(p)
  {
  }

  template <class T>
  inline const T *
  as() const
  {
    return boost::get<const T *>(handle);
  }
};

///
/// PrimitiveBox. Type used for CGAL::Box_intersection_d
template <int Dim>
struct PrimitiveBox {
  typedef CGAL::Box_intersection_d::Box_with_handle_d<double, Dim,
                                                      PrimitiveHandle<Dim> *>
      Type;
};

///
/// BoxCollection for use with CGAL::Box_intersection_d
template <int Dim>
struct BoxCollection {
  typedef std::vector<typename PrimitiveBox<Dim>::Type> Type;
};

///
/// HandleCollection. Used to store PrimitiveHandle
template <int Dim>
struct HandleCollection {
  typedef std::list<PrimitiveHandle<Dim>> Type;
};

///
/// Flags available for each type of Geometry type.
/// Primitives can be 'flagged' in order to speed up recomposition
enum ElementFlag {
  // the polyhedron is planar => build a triangle or a polygon
  FLAG_IS_PLANAR = 1
};

///
/// CollectionElement, a Primitive with flags
/// Primitive : Point_d, Segment_d, Surface_d, Volume_d
template <class Primitive>
class CollectionElement {
public:
  int
  flags() const
  {
    return _flags;
  }
  void
  setFlags(int flags)
  {
    _flags = flags;
  }

  Primitive &
  primitive()
  {
    return _primitive;
  }
  const Primitive &
  primitive() const
  {
    return _primitive;
  }

  // constructor from Primitive
  CollectionElement() : _flags(0) {}
  CollectionElement(const Primitive &p) : _primitive(p), _flags(0) {}
  CollectionElement(const Primitive &p, int f) : _primitive(p), _flags(f) {}

  CollectionElement(const CollectionElement &other)
      : _primitive(other._primitive), _flags(other._flags)
  {
  }
  bool
  operator<(const CollectionElement &other) const
  {
    return _primitive < other._primitive;
  }

private:
  Primitive _primitive;
  int       _flags;
};

template <class Primitive>
std::ostream &
operator<<(std::ostream &ostr, const CollectionElement<Primitive> &p)
{
  ostr << p.primitive() << " flags: " << p.flags();
  return ostr;
}

///
/// A GeometrySet represents a set of CGAL primitives.
/// Primitive are either of dimension 0 (points),
/// dimension 1 (segments), dimension 2 (surfaces, a.k.a. polygon or triangles)
/// or dimension 3 (polyhedron)
template <int Dim>
class GeometrySet {
public:
  // Points are stored in an ordered set
  typedef std::set<CollectionElement<typename Point_d<Dim>::Type>>
      PointCollection;
  // Segments are stored in an ordered set
  typedef std::set<CollectionElement<typename Segment_d<Dim>::Type>>
      SegmentCollection;
  typedef std::list<CollectionElement<typename Surface_d<Dim>::Type>>
      SurfaceCollection;
  typedef std::list<CollectionElement<typename Volume_d<Dim>::Type>>
      VolumeCollection;

  GeometrySet();

  /**
   * Construct a GeometrySet from a SFCGAL::Geometry
   */
  GeometrySet(const Geometry &g);

  /**
   * Construct a GeometrySet from a Point
   */
  GeometrySet(const typename TypeForDimension<Dim>::Point &g, int flags = 0);

  /**
   * Construct a GeometrySet from a Segment
   */
  GeometrySet(const typename TypeForDimension<Dim>::Segment &g, int flags = 0);

  /**
   * Construct a GeometrySet from a Surface
   */
  GeometrySet(const typename TypeForDimension<Dim>::Surface &g, int flags = 0);

  /**
   * Construct a GeometrySet from a Volume
   */
  GeometrySet(const typename TypeForDimension<Dim>::Volume &g, int flags = 0);

  /**
   * Add primitives from another set
   */
  void
  merge(const GeometrySet<Dim> &g);

  /**
   * Add a geometry by decomposing it into CGAL primitives
   */
  void
  addGeometry(const Geometry &g);

  /**
   * add a primitive from a PrimitiveHandle  to the set
   */
  void
  addPrimitive(const PrimitiveHandle<Dim> &p);

  /**
   * add a primitive from a CGAL::Object to the set
   * pointsAsRing : if set to true, build a polygon if o is a vector of points
   */
  void
  addPrimitive(const CGAL::Object &o, bool pointsAsRing = false);

  /**
   * add a point to the set
   */
  void
  addPrimitive(const typename TypeForDimension<Dim>::Point &g, int flags = 0);
  template <class IT>
  void
  addPoints(IT ibegin, IT iend)
  {
    std::copy(ibegin, iend, std::inserter(_points, _points.end()));
  }

  /**
   * collect all points of b and add them to the point list
   */
  void
  collectPoints(const PrimitiveHandle<Dim> &b);

  /**
   * add a segment to the set
   */
  void
  addPrimitive(const typename TypeForDimension<Dim>::Segment &g, int flags = 0);
  template <class IT>
  void
  addSegments(IT ibegin, IT iend)
  {
    std::copy(ibegin, iend, std::inserter(_segments, _segments.end()));
  }

  /**
   * add a surface to the set
   */
  void
  addPrimitive(const typename TypeForDimension<Dim>::Surface &g, int flags = 0);
  template <class IT>
  void
  addSurfaces(IT ibegin, IT iend)
  {
    std::copy(ibegin, iend, std::back_inserter(_surfaces));
  }

  /**
   * add a volume to the set
   */
  void
  addPrimitive(const typename TypeForDimension<Dim>::Volume &g, int flags = 0);
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
   */
  int
  dimension() const;

  /**
   * Add the boundary (segments) of a surface
   */
  void
  addBoundary(const typename TypeForDimension<Dim>::Surface &surface);

  /**
   * Add the boundary (surfaces) of a volume
   */
  void
  addBoundary(const typename TypeForDimension<Dim>::Volume &volume);

  /**
   * Compute all bounding boxes and handles of the set
   */
  void
  computeBoundingBoxes(typename HandleCollection<Dim>::Type &handles,
                       typename BoxCollection<Dim>::Type    &boxes) const;

  inline PointCollection &
  points()
  {
    return _points;
  }
  inline const PointCollection &
  points() const
  {
    return _points;
  }

  inline SegmentCollection &
  segments()
  {
    return _segments;
  }
  inline const SegmentCollection &
  segments() const
  {
    return _segments;
  }

  inline SurfaceCollection &
  surfaces()
  {
    return _surfaces;
  }
  inline const SurfaceCollection &
  surfaces() const
  {
    return _surfaces;
  }

  inline VolumeCollection &
  volumes()
  {
    return _volumes;
  }
  inline const VolumeCollection &
  volumes() const
  {
    return _volumes;
  }

  /**
   * Returns true if the set holds points
   */
  bool
  hasPoints() const;
  /**
   * Returns true if the set holds segments
   */
  bool
  hasSegments() const;
  /**
   * Returns true if the set holds surfaces
   */
  bool
  hasSurfaces() const;
  /**
   * Returns true if the set holds volumes
   */
  bool
  hasVolumes() const;

  /**
   * convert the set to a SFCGAL::Geometry
   */
  std::unique_ptr<Geometry>
  recompose() const;

  /**
   * Filter (remove) primitives that are already covered by others
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

///
/// Display operator
SFCGAL_API std::ostream            &
operator<<(std::ostream &, const GeometrySet<2> &g);
///
/// Display operator
SFCGAL_API std::ostream            &
operator<<(std::ostream &, const GeometrySet<3> &g);

// bbox of a 'volume' for 2D, will never be called
inline CGAL::Bbox_2
compute_solid_bbox(const NoVolume &, dim_t<2>)
{
  return CGAL::Bbox_2();
}

inline CGAL::Bbox_3
compute_solid_bbox(const TypeForDimension<3>::Volume &vol, dim_t<3>)
{
  BOOST_ASSERT(vol.size_of_vertices());
  MarkedPolyhedron::Point_const_iterator pit = vol.points_begin();
  CGAL::Bbox_3                           ret(pit->bbox());
  ++pit;

  for (; pit != vol.points_end(); ++pit) {
    ret = ret + pit->bbox();
  }

  return ret;
}
} // namespace detail
} // namespace SFCGAL

#endif
