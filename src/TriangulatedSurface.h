// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TRIANGULATED_SURFACE_H_
#define SFCGAL_TRIANGULATED_SURFACE_H_

#include <boost/assert.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/unique_ptr.hpp>
#include <set>
#include <vector>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include "SFCGAL/DereferenceIterator.h"
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
  using iterator =
      DereferenceIterator<std::vector<std::unique_ptr<Triangle>>::iterator>;
  using const_iterator = DereferenceIterator<
      std::vector<std::unique_ptr<Triangle>>::const_iterator>;

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
   * [SFA/OGC]Returns the number of patches
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
   * @brief Adds a patch to the TriangulatedSurface.
   *
   * @param patch The Triangle object representing the new patch to add.
   */
  void
  addPatch(const Triangle &patch)
  {
    addPatch(std::unique_ptr<Triangle>(patch.clone()));
  }
  /**
   * @brief Adds a patch to the TriangulatedSurface.
   *
   * @param patch A raw pointer to the Triangle object representing the new
   * patch to add. Ownership of this patch is transferred to the
   * TriangulatedSurface.
   *
   * @deprecated The unique_ptr version should be used instead
   */
  void
  addPatch(Triangle *patch)
  {
    addPatch(std::unique_ptr<Triangle>(patch));
  }
  /**
   * @brief Adds a patch to the TriangulatedSurface.
   *
   * @param patch A unique pointer to the Triangle object representing the new
   * patch to add. Ownership of this patch is moved into the
   * TriangulatedSurface.
   */
  void
  addPatch(std::unique_ptr<Triangle> patch)
  {
    _triangles.push_back(std::move(patch));
  }
  /**
   * add patchs from an other TriangulatedSurface
   */
  void
  addPatchs(const TriangulatedSurface &other);

  /**
   * [SFA/OGC]Returns the number of points
   * @deprecated see numPatches()
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
    return *_triangles[n];
  }
  /**
   * [SFA/OGC]Returns the n-th point
   * @deprecated see patchN()
   */
  inline Triangle &
  triangleN(size_t const &n)
  {
    BOOST_ASSERT(n < _triangles.size());
    return *_triangles[n];
  }
  /**
   * add a Triangle to the TriangulatedSurface
   * @deprecated see addPatch()
   */
  inline void
  addTriangle(const Triangle &triangle)
  {
    addTriangle(std::unique_ptr<Triangle>(triangle.clone()));
  }
  /**
   * add a Triangle to the TriangulatedSurface
   * @deprecated see addPatch()
   */
  inline void
  addTriangle(Triangle *triangle)
  {
    addTriangle(std::unique_ptr<Triangle>(triangle));
  }
  /**
   * add a Triangle to the TriangulatedSurface
   * @deprecated see addPatch()
   */
  void
  addTriangle(std::unique_ptr<Triangle> triangle)
  {
    _triangles.push_back(std::move(triangle));
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
  setPatchN(const Geometry &geometry, size_t const &idx);

  /**
   * Sets the n-th Geometry, starting at zero
   * It needs to be a triangle.
   * The ownership of the polygon is taken. The caller is not responsible
   * anymore of its deallocation.
   */
  void
  setPatchN(Geometry *geometry, size_t const &idx);

  /**
   * @brief Sets the n-th patch in the TriangulatedSurface.
   *
   * Replaces the patch at the specified index with the provided one.
   * Only Triangle geometries are allowed. Indexing starts at zero.
   *
   * @param geometry A unique pointer to the Geometry object to set. Ownership
   * is transferred to this class.
   * @param idx The zero-based index of the patch to replace.
   *
   */
  void
  setPatchN(std::unique_ptr<Geometry> geometry, size_t const &idx);

  /**
   * Sets the n-th Patch, starting at zero
   */
  void
  setPatchN(const Triangle &triangle, size_t const &idx);

  /**
   * Sets the n-th Patch, starting at zero
   * The ownership of the polygon is taken. The caller is not responsible
   * anymore of its deallocation.
   *
   * @deprecated The unique_ptr version should be used instead
   */
  void
  setPatchN(Triangle *triangle, size_t const &idx);

  /**
   * @brief Sets the n-th patch in the TriangulatedSurface.
   *
   * Replaces the patch at the specified index with the provided one.
   * Indexing starts at zero.
   *
   * @param triangle A unique pointer to the Triangle object to set. Ownership
   * is transferred to this class.
   * @param idx The zero-based index of the patch to replace.
   *
   */
  void
  setPatchN(std::unique_ptr<Triangle> triangle, size_t const &idx);

  //-- optimization

  void
  reserve(const size_t &n);

  //-- iterators

  inline iterator
  begin()
  {
    return dereference_iterator(_triangles.begin());
  }
  inline const_iterator
  begin() const
  {
    return dereference_iterator(_triangles.begin());
  }

  inline iterator
  end()
  {
    return dereference_iterator(_triangles.end());
  }
  inline const_iterator
  end() const
  {
    return dereference_iterator(_triangles.end());
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
  std::vector<std::unique_ptr<Triangle>> _triangles;

  void
  swap(TriangulatedSurface &other)
  {
    std::swap(_triangles, other._triangles);
  }
};
} // namespace SFCGAL

#endif
