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
class SFCGAL_API TriangulatedSurface
    : public GeometryImpl<TriangulatedSurface, Surface> {
public:
  /// @brief Iterator type for triangulated surface patches
  using iterator =
      DereferenceIterator<std::vector<std::unique_ptr<Triangle>>::iterator>;
  /// @brief Const iterator type for triangulated surface patches
  using const_iterator = DereferenceIterator<
      std::vector<std::unique_ptr<Triangle>>::const_iterator>;

  /**
   * Empty TriangulatedSurface constructor
   */
  TriangulatedSurface();
  /**
   * Constructor with a vector of triangles
   * @param triangle Vector of triangles to initialize the surface
   */
  TriangulatedSurface(const std::vector<Triangle> &triangle);
  /**
   * Copy constructor
   * @param other The triangulated surface to copy from
   */
  TriangulatedSurface(const TriangulatedSurface &other);
  /**
   * assign operator
   * @param other The triangulated surface to assign from
   * @return Reference to this triangulated surface
   */
  auto
  operator=(TriangulatedSurface other) -> TriangulatedSurface &;
  /**
   * destructor
   */
  ~TriangulatedSurface() override;

  //-- SFCGAL::Geometry
  /// @brief Get the geometry type as string
  /// @return "TriangulatedSurface"
  [[nodiscard]] auto
  geometryType() const -> std::string override;
  //-- SFCGAL::Geometry
  /// @brief Get the geometry type identifier
  /// @return TYPE_TRIANGULATEDSURFACE
  [[nodiscard]] auto
  geometryTypeId() const -> GeometryType override;
  //-- SFCGAL::Geometry
  /// @brief Get the dimension of the surface
  /// @return 2 (surfaces are 2-dimensional)
  [[nodiscard]] auto
  dimension() const -> int override;
  //-- SFCGAL::Geometry
  /// @brief Get the coordinate dimension
  /// @return Number of coordinates per point
  [[nodiscard]] auto
  coordinateDimension() const -> int override;
  //-- SFCGAL::Geometry
  /// @brief Check if the surface is empty
  /// @return true if empty, false otherwise
  [[nodiscard]] auto
  isEmpty() const -> bool override;
  //-- SFCGAL::Geometry
  /// @brief Check if the surface has 3D coordinates
  /// @return true if 3D, false otherwise
  [[nodiscard]] auto
  is3D() const -> bool override;
  //-- SFCGAL::Geometry
  /// @brief Check if the surface has measured coordinates
  /// @return true if measured, false otherwise
  [[nodiscard]] auto
  isMeasured() const -> bool override;

  /// @brief Drop Z coordinate from all points
  /// @return true if Z was dropped, false otherwise
  auto
  dropZ() -> bool override;

  /// @brief Drop M coordinate from all points
  /// @return true if M was dropped, false otherwise
  auto
  dropM() -> bool override;

  /// @brief Swap X and Y coordinates of all points
  auto
  swapXY() -> void override;

  /**
   * [SFA/OGC]Returns the number of patches
   * @return Number of triangular patches
   * @deprecated see numGeometries()
   */
  [[nodiscard]] auto
  numPatches() const -> size_t
  {
    return _triangles.size();
  }
  /**
   * [SFA/OGC]Returns the n-th patch
   * @param n The index of the patch to get
   * @return Const reference to the nth triangle patch
   */
  [[nodiscard]] auto
  patchN(size_t const &n) const -> const Triangle &;
  /**
   * [SFA/OGC]Returns the n-th patch
   * @param n The index of the patch to get
   * @return Reference to the nth triangle patch
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
    addPatch(patch.clone());
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
    _triangles.emplace_back(std::move(patch));
  }
  /**
   * add patchs from an other TriangulatedSurface
   * @param other The triangulated surface to add patches from
   */
  void
  addPatchs(const TriangulatedSurface &other);

  /**
   * [SFA/OGC]Returns the number of triangles
   * @return Number of triangles in the surface
   * @deprecated see numPatches()
   */
  [[nodiscard]] auto
  numTriangles() const -> size_t
  {
    return _triangles.size();
  }
  /**
   * [SFA/OGC]Returns the n-th triangle
   * @param n The triangle index
   * @return Const reference to the nth triangle
   * @deprecated see patchN()
   */
  [[nodiscard]] auto
  triangleN(size_t const &n) const -> const Triangle &
  {
    BOOST_ASSERT(n < _triangles.size());
    return *_triangles[n];
  }
  /**
   * [SFA/OGC]Returns the n-th triangle
   * @param n The triangle index
   * @return Reference to the nth triangle
   * @deprecated see patchN()
   */
  auto
  triangleN(size_t const &n) -> Triangle &
  {
    BOOST_ASSERT(n < _triangles.size());
    return *_triangles[n];
  }
  /**
   * add a Triangle to the TriangulatedSurface
   * @param triangle The triangle to add
   * @deprecated see addPatch()
   */
  void
  addTriangle(const Triangle &triangle)
  {
    addTriangle(triangle.clone());
  }
  /**
   * add a Triangle to the TriangulatedSurface
   * @param triangle Pointer to the triangle to add
   * @deprecated see addPatch()
   */
  void
  addTriangle(Triangle *triangle)
  {
    addTriangle(std::unique_ptr<Triangle>(triangle));
  }
  /**
   * add a Triangle to the TriangulatedSurface
   * @param triangle Unique pointer to the triangle to add
   * @deprecated see addPatch()
   */
  void
  addTriangle(std::unique_ptr<Triangle> triangle)
  {
    _triangles.emplace_back(std::move(triangle));
  }
  /**
   * add triangles from an other TriangulatedSurface
   * @param other The triangulated surface to copy triangles from
   * @deprecated see addPatchs()
   */
  void
  addTriangles(const TriangulatedSurface &other);

  /**
   * Sets the n-th Geometry, starting at zero
   * It needs to be a triangle.
   * @param geometry The geometry to set (must be a triangle)
   * @param idx The index of the patch to set
   */
  void
  setPatchN(const Geometry &geometry, size_t const &idx);

  /**
   * Sets the n-th Geometry, starting at zero
   * It needs to be a triangle.
   * The ownership of the polygon is taken. The caller is not responsible
   * anymore of its deallocation.
   * @param geometry Pointer to the geometry to set (must be a triangle)
   * @param idx The index of the patch to set
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
   * @param triangle The triangle to set
   * @param idx The index of the patch to set
   */
  void
  setPatchN(const Triangle &triangle, size_t const &idx);

  /**
   * Sets the n-th Patch, starting at zero
   * The ownership of the polygon is taken. The caller is not responsible
   * anymore of its deallocation.
   * @param triangle Pointer to the triangle to set
   * @param idx The index of the patch to set
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

  /**
   * @brief Reserve space for triangles
   * @param n Number of triangles to reserve space for
   */
  void
  reserve(const size_t &n);

  //-- iterators

  /**
   * @brief Get iterator to beginning of patches
   * @return Iterator to first patch
   */
  auto
  begin() -> iterator
  {
    return dereference_iterator(_triangles.begin());
  }
  /**
   * @brief Get const iterator to beginning of patches
   * @return Const iterator to first patch
   */
  [[nodiscard]] auto
  begin() const -> const_iterator
  {
    return dereference_iterator(_triangles.begin());
  }

  /**
   * @brief Get iterator to end of patches
   * @return Iterator to past-the-end
   */
  auto
  end() -> iterator
  {
    return dereference_iterator(_triangles.end());
  }
  /**
   * @brief Get const iterator to end of patches
   * @return Const iterator to past-the-end
   */
  [[nodiscard]] auto
  end() const -> const_iterator
  {
    return dereference_iterator(_triangles.end());
  }

  //-- visitors

  //-- SFCGAL::Geometry
  /// @brief Accept a geometry visitor
  /// @param visitor Visitor to accept
  void
  accept(GeometryVisitor &visitor) override;
  //-- SFCGAL::Geometry
  /// @brief Accept a const geometry visitor
  /// @param visitor Const visitor to accept
  void
  accept(ConstGeometryVisitor &visitor) const override;

  //-- helpers

  /**
   * @brief Converts the current TriangulatedSurface into a CGAL::Polyhedron_3.
   *
   * This function creates a new `Polyhedron` instance that represents the
   * geometry of the TriangulatedSurface.
   *
   * @tparam Polyhedron The CGAL Polyhedron_3 type.
   *
   * @return std::unique_ptr<Polyhedron> A unique pointer to the newly created
   * Polyhedron.
   */
  template <typename Polyhedron>
  auto
  toPolyhedron_3() const -> std::unique_ptr<Polyhedron>;

  /**
   * Serializer
   * @param ar Archive for serialization
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
  swap(TriangulatedSurface &other) noexcept
  {
    std::swap(_triangles, other._triangles);
  }
};
} // namespace SFCGAL

#endif
