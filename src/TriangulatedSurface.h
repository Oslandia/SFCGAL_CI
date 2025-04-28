// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TRIANGULATED_SURFACE_H_
#define SFCGAL_TRIANGULATED_SURFACE_H_

#include <boost/assert.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/serialize_ptr_vector.hpp>
#include <boost/serialization/base_object.hpp>
#include <set>
#include <vector>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include "SFCGAL/Exception.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Triangle.h"

namespace SFCGAL {

/**
 * A TriangulatedSurface in SFA modeled as a Triangle soup
 * @todo do better than a "triangle soup" or add topological view?
 */
class SFCGAL_API TriangulatedSurface : public Surface {
public:
  typedef boost::ptr_vector<Triangle>::iterator       iterator;
  typedef boost::ptr_vector<Triangle>::const_iterator const_iterator;

  /**
   * Empty TriangulatedSurface constructor
   */
  TriangulatedSurface();
  /**
   * Constructor with a vector of triangles
   */
  TriangulatedSurface(const std::vector<Triangle> &triangle);
  /**
   * Copy constructor
   */
  TriangulatedSurface(const TriangulatedSurface &other);
  /**
   * assign operator
   */
  TriangulatedSurface &
  operator=(TriangulatedSurface other);
  /**
   * destructor
   */
  ~TriangulatedSurface();

  //-- SFCGAL::Geometry
  TriangulatedSurface *
  clone() const override;

  //-- SFCGAL::Geometry
  std::string
  geometryType() const override;
  //-- SFCGAL::Geometry
  GeometryType
  geometryTypeId() const override;
  //-- SFCGAL::Geometry
  int
  dimension() const override;
  //-- SFCGAL::Geometry
  int
  coordinateDimension() const override;
  //-- SFCGAL::Geometry
  bool
  isEmpty() const override;
  //-- SFCGAL::Geometry
  bool
  is3D() const override;
  //-- SFCGAL::Geometry
  bool
  isMeasured() const override;

  auto
  dropZ() -> bool override;

  auto
  dropM() -> bool override;

  auto
  swapXY() -> void override;

  /**
   * [SFA/OGC]Returns the number of patchs
   * @deprecated see numGeometries()
   */
  inline auto
  numPatches() const -> size_t
  {
    return _triangles.size();
  }
  /**
   * [SFA/OGC]Returns the n-th patch
   */
  auto
  patchN(size_t const &n) const -> const Triangle &;
  /**
   * [SFA/OGC]Returns the n-th patch
   */
  auto
  patchN(size_t const &n) -> Triangle &;
  /**
   * add a Patch to the TriangulatedSurface
   */
  inline void
  addPatch(const Triangle &patch)
  {
    addPatch(patch.clone());
  }
  /**
   * add a Patch to the TriangulatedSurface
   */
  inline void
  addPatch(Triangle *patch)
  {
    _triangles.push_back(patch);
  }
  /**
   * add patchs from an other TriangulatedSurface
   */
  void
  addPatchs(const TriangulatedSurface &other);

  /**
   * [SFA/OGC]Returns the number of points
   * @deprecated see numPatchs()
   */
  inline size_t
  numTriangles() const
  {
    return _triangles.size();
  }
  /**
   * [SFA/OGC]Returns the n-th point
   * @deprecated see patchN()
   */
  inline const Triangle &
  triangleN(size_t const &n) const
  {
    BOOST_ASSERT(n < _triangles.size());
    return _triangles[n];
  }
  /**
   * [SFA/OGC]Returns the n-th point
   * @deprecated see patchN()
   */
  inline Triangle &
  triangleN(size_t const &n)
  {
    BOOST_ASSERT(n < _triangles.size());
    return _triangles[n];
  }
  /**
   * add a Triangle to the TriangulatedSurface
   * @deprecated see addPatch()
   */
  inline void
  addTriangle(const Triangle &triangle)
  {
    addTriangle(triangle.clone());
  }
  /**
   * add a Triangle to the TriangulatedSurface
   * @deprecated see addPatch()
   */
  inline void
  addTriangle(Triangle *triangle)
  {
    _triangles.push_back(triangle);
  }
  /**
   * add triangles from an other TriangulatedSurface
   * @deprecated see addPatchs()
   */
  void
  addTriangles(const TriangulatedSurface &other);

  /**
   * Sets the n-th Geometry, starting at zero
   * It needs to be a triangle.
   */
  void
  setPatchN(const Geometry &geometry, size_t const &n);

  /**
   * Sets the n-th Geometry, starting at zero
   * It needs to be a triangle.
   * The ownership of the polygon is taken. The caller is not responsible
   * anymore of its deallocation.
   */
  void
  setPatchN(Geometry *geometry, size_t const &n);

  /**
   * Sets the n-th Patch, starting at zero
   */
  void
  setPatchN(const Triangle &triangle, size_t const &n);

  /**
   * Sets the n-th Patch, starting at zero
   * The ownership of the polygon is taken. The caller is not responsible
   * anymore of its deallocation.
   */
  void
  setPatchN(Triangle *triangle, size_t const &n);

  //-- optimization

  void
  reserve(const size_t &n);

  //-- iterators

  inline iterator
  begin()
  {
    return _triangles.begin();
  }
  inline const_iterator
  begin() const
  {
    return _triangles.begin();
  }

  inline iterator
  end()
  {
    return _triangles.end();
  }
  inline const_iterator
  end() const
  {
    return _triangles.end();
  }

  //-- visitors

  //-- SFCGAL::Geometry
  void
  accept(GeometryVisitor &visitor) override;
  //-- SFCGAL::Geometry
  void
  accept(ConstGeometryVisitor &visitor) const override;

  //-- helpers

  /**
   * @brief Converts a TriangulatedSurface to a CGAL::Polyhedron_3
   */
  template <typename K, typename Polyhedron>
  std::unique_ptr<Polyhedron>
  toPolyhedron_3() const;

  /**
   * Serializer
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int /*version*/)
  {
    ar &boost::serialization::base_object<Geometry>(*this);
    ar & _triangles;
  }

private:
  boost::ptr_vector<Triangle> _triangles;

  void
  swap(TriangulatedSurface &other)
  {
    std::swap(_triangles, other._triangles);
  }
};
} // namespace SFCGAL

#endif
