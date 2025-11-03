// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_POLYHEDRALSURFACE_H_
#define SFCGAL_POLYHEDRALSURFACE_H_

#include <boost/assert.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/unique_ptr.hpp>
#include <memory>
#include <vector>

#include "SFCGAL/DereferenceIterator.h"
#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/triangulate/triangulatePolygon.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

using InexactKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh          = CGAL::Surface_mesh<SFCGAL::Kernel::Point_3>;
using InexactMesh   = CGAL::Surface_mesh<InexactKernel::Point_3>;

namespace SFCGAL {

/**
 * A PolyhedralSurface in SFA modeled as a Polygon soup
 * @todo do better than a "polygon soup" or add topological view?
 */
class SFCGAL_API PolyhedralSurface
    : public GeometryImpl<PolyhedralSurface, Surface> {
public:
  /// @brief Iterator type for polyhedral surface patches
  using iterator =
      DereferenceIterator<std::vector<std::unique_ptr<Polygon>>::iterator>;
  /// @brief Const iterator type for polyhedral surface patches
  using const_iterator = DereferenceIterator<
      std::vector<std::unique_ptr<Polygon>>::const_iterator>;

  /**
   * Empty PolyhedralSurface constructor
   */
  PolyhedralSurface();
  /**
   * Constructor with a vector of polygons
   * @param polygons Vector of polygons to initialize the surface
   */
  PolyhedralSurface(const std::vector<Polygon> &polygons);

  /**
   * Constructor with a Geometry
   * @param geometry The geometry to convert to polyhedral surface
   */
  PolyhedralSurface(const std::unique_ptr<Geometry> &geometry);

  /**
   * Constructor from a Polyhedron (detail::MarkedPolyhedron or
   * CGAL::Polyhedron_3)
   * @tparam Polyhedron The CGAL polyhedron type
   * @param poly The polyhedron to convert
   */
  template <typename Polyhedron>
  PolyhedralSurface(const Polyhedron &poly)
  {
    for (typename Polyhedron::Facet_const_iterator fit = poly.facets_begin();
         fit != poly.facets_end(); ++fit) {
      auto *face = new LineString();
      typename Polyhedron::Halfedge_around_facet_const_circulator hit =
          fit->facet_begin();
      do {
        face->addPoint(hit->vertex()->point());
        ++hit;
      } while (hit != fit->facet_begin());
      // close the ring
      face->addPoint(hit->vertex()->point());
      _polygons.push_back(std::make_unique<Polygon>(face));
    }
  }

  /**
   * Constructor from a CGAL::Surface_mesh
   * @param sm The surface mesh to convert
   */
  PolyhedralSurface(const Mesh &sm);

  /**
   * Constructor from an inexact CGAL::Surface_mesh
   * @param inexactMesh The inexact surface mesh to convert
   */
  PolyhedralSurface(const InexactMesh &inexactMesh);
  /**
   * Copy constructor
   * @param other The polyhedral surface to copy from
   */
  PolyhedralSurface(const PolyhedralSurface &other);
  /**
   * assign operator
   * @param other The polyhedral surface to assign from
   * @return Reference to this polyhedral surface
   */
  auto
  operator=(PolyhedralSurface other) -> PolyhedralSurface &;
  /**
   * destructor
   */
  ~PolyhedralSurface() override;

  //-- SFCGAL::Geometry
  /// @brief Get the geometry type as string
  /// @return "PolyhedralSurface"
  [[nodiscard]] auto
  geometryType() const -> std::string override;
  //-- SFCGAL::Geometry
  /// @brief Get the geometry type identifier
  /// @return TYPE_POLYHEDRALSURFACE
  [[nodiscard]] auto
  geometryTypeId() const -> GeometryType override;
  //-- SFCGAL::Geometry
  /// @brief Get the dimension of the polyhedral surface
  /// @return 2 (surfaces are 2-dimensional)
  [[nodiscard]] auto
  dimension() const -> int override;
  //-- SFCGAL::Geometry
  /// @brief Get the coordinate dimension
  /// @return Number of coordinates per point
  [[nodiscard]] auto
  coordinateDimension() const -> int override;
  //-- SFCGAL::Geometry
  /// @brief Check if the polyhedral surface is empty
  /// @return true if empty, false otherwise
  [[nodiscard]] auto
  isEmpty() const -> bool override;
  //-- SFCGAL::Geometry
  /// @brief Check if the polyhedral surface has 3D coordinates
  /// @return true if 3D, false otherwise
  [[nodiscard]] auto
  is3D() const -> bool override;
  //-- SFCGAL::Geometry
  /// @brief Check if the polyhedral surface has measured coordinates
  /// @return true if measured, false otherwise
  [[nodiscard]] auto
  isMeasured() const -> bool override;

  /// @brief Drop Z coordinate from all polygons
  /// @return true if Z was dropped, false otherwise
  auto
  dropZ() -> bool override;

  /// @brief Drop M coordinate from all polygons
  /// @return true if M was dropped, false otherwise
  auto
  dropM() -> bool override;

  /// @brief Swap X and Y coordinates of all polygons
  auto
  swapXY() -> void override;

  /**
   * Convert PolyhedralSurface to TriangulatedSurface
   * @return A triangulated surface representation
   */
  [[nodiscard]] auto
  toTriangulatedSurface() const -> TriangulatedSurface;

  /**
   * [SFA/OGC]Returns the number of patches
   * @return Number of polygonal patches in the surface
   * @warning PolyhedralSurface is treated as one geometry. numGeometries
   * returns 1 or 0 for empty PolyhedralSurface
   */
  [[nodiscard]] auto
  numPatches() const -> size_t
  {
    return _polygons.size();
  }

  /**
   * @brief [SFA/OGC]Returns the number of polygons
   * @return Number of polygons in the surface
   * @warning PolyhedralSurface is treated as one geometry. numGeometries
   * returns 1 or 0 for empty PolyhedralSurface
   * @deprecated see numPatches
   * @see numGeometries()
   */
  [[nodiscard]] auto
  numPolygons() const -> size_t
  {
    return numPatches();
  }

  /**
   * [SFA/OGC]Returns the n-th patch
   * @param n The index of the patch to get
   * @return Const reference to the nth patch
   */
  [[nodiscard]] auto
  patchN(size_t const &n) const -> const Polygon &
  {
    BOOST_ASSERT(n < _polygons.size());
    return *_polygons[n];
  }
  /**
   * [SFA/OGC]Returns the n-th patch
   * @param n The index of the patch to get
   * @return Reference to the nth patch
   */
  auto
  patchN(size_t const &n) -> Polygon &
  {
    BOOST_ASSERT(n < _polygons.size());
    return *_polygons[n];
  }
  /**
   * @brief Adds a polygonal patch to the PolyhedralSurface.
   *
   * @param patch The Polygon object representing the new patch to add.
   */
  void
  addPatch(const Polygon &patch);
  /**
   * @brief Adds a polygonal patch to the PolyhedralSurface.
   *
   * @param patch A raw pointer to the Polygon object representing the new
   * patch to add. Ownership of this patch is moved into the PolyhedralSurface.
   *
   * @deprecated The unique_ptr version should be used instead
   */
  void
  addPatch(Polygon *patch);
  /**
   * @brief Adds a polygonal patch to the PolyhedralSurface.
   *
   * @param patch A unique pointer to the Polygon object representing the new
   * patch to add. Ownership of this patch is moved into the PolyhedralSurface.
   */
  void
  addPatch(std::unique_ptr<Polygon> patch);
  /**
   * add patches from an other PolyhedralSurface
   * @param polyhedralSurface The polyhedral surface to add patches from
   */
  void
  addPatchs(const PolyhedralSurface &polyhedralSurface);

  /**
   * @brief [SFA/OGC]Returns the n-th polygon
   * @param n Index of the polygon to get
   * @return Const reference to the nth polygon
   * @deprecated see patchN()
   */
  [[nodiscard]] auto
  polygonN(size_t const &n) const -> const Polygon &
  {
    return patchN(n);
  }
  /**
   * @brief [SFA/OGC]Returns the n-th polygon
   * @param n Index of the polygon to get
   * @return Reference to the nth polygon
   * @deprecated see patchN()
   */
  auto
  polygonN(size_t const &n) -> Polygon &
  {
    return patchN(n);
  }
  /**
   * add a polygon to the PolyhedralSurface
   * @param polygon The polygon to add
   * @deprecated see addPatch()
   */
  void
  addPolygon(const Polygon &polygon);
  /**
   * add a polygon to the PolyhedralSurface
   * @param polygon The polygon to add
   * @deprecated see addPatch()
   */
  void
  addPolygon(Polygon *polygon);
  /**
   * @brief add polygons from an other PolyhedralSurface
   * @param polyhedralSurface The polyhedral surface to add polygons from
   * @deprecated see addPatchs()
   */
  void
  addPolygons(const PolyhedralSurface &polyhedralSurface);

  /**
   * @brief Sets the n-th Geometry, starting at zero
   * @param geometry Geometry to set (must be a polygon)
   * @param idx Index of the patch to set
   */
  void
  setPatchN(const Geometry &geometry, size_t const &idx);

  /**
   * @brief Sets the n-th Geometry, starting at zero
   * @param geometry Pointer to geometry (must be a polygon)
   * @param idx Index of the patch to set
   * @note The ownership of the polygon is taken. The caller is not responsible
   * anymore of its deallocation.
   * @deprecated The unique_ptr version should be used instead
   */
  void
  setPatchN(Geometry *geometry, size_t const &idx);

  /**
   * @brief Sets the n-th patch in the PolyhedralSurface.
   *
   * Replaces the patch at the specified index with the provided one.
   * Only polygon geometries are allowed. Indexing starts at zero.
   *
   * @param geometry A unique pointer to the Geometry object to set. Ownership
   * is transferred to this class.
   * @param idx The zero-based index of the patch to replace.
   *
   */
  void
  setPatchN(std::unique_ptr<Geometry> geometry, size_t const &idx);

  /**
   * @brief Sets the n-th Patch, starting at zero
   * @param patch Polygon patch to set
   * @param idx Index of the patch to set
   */
  void
  setPatchN(const Polygon &patch, size_t const &idx);

  /**
   * @brief Sets the n-th Patch, starting at zero
   * @param patch Pointer to polygon patch
   * @param idx Index of the patch to set
   * @note The ownership of the polygon is taken. The caller is not responsible
   * anymore of its deallocation.
   * @deprecated The unique_ptr version should be used instead
   */
  void
  setPatchN(Polygon *patch, size_t const &idx);

  /**
   * @brief Sets the n-th patch in the PolyhedralSurface.
   *
   * Replaces the patch at the specified index with the provided one.
   * Indexing starts at zero.
   *
   * @param patch A unique pointer to the Polygon object to set. Ownership is
   * transferred to this class.
   * @param idx The zero-based index of the patch to replace.
   *
   */
  void
  setPatchN(std::unique_ptr<Polygon> patch, size_t const &idx);

  /**
   * @brief Convert to CGAL::Polyhedron_3
   * @return Unique pointer to CGAL polyhedron representation
   */
  template <typename Polyhedron>
  auto
  toPolyhedron_3() const -> std::unique_ptr<Polyhedron>
  {
    TriangulatedSurface tri;
    triangulate::triangulatePolygon3D(*this, tri);
    return tri.toPolyhedron_3<Polyhedron>();
  }

  //-- iterators

  /**
   * @brief Get iterator to beginning of patches
   * @return Iterator to first patch
   */
  auto
  begin() -> iterator
  {
    return dereference_iterator(_polygons.begin());
  }
  /**
   * @brief Get const iterator to beginning of patches
   * @return Const iterator to first patch
   */
  [[nodiscard]] auto
  begin() const -> const_iterator
  {
    return dereference_iterator(_polygons.begin());
  }

  /**
   * @brief Get iterator to end of patches
   * @return Iterator to past-the-end
   */
  auto
  end() -> iterator
  {
    return dereference_iterator(_polygons.end());
  }
  /**
   * @brief Get const iterator to end of patches
   * @return Const iterator to past-the-end
   */
  [[nodiscard]] auto
  end() const -> const_iterator
  {
    return dereference_iterator(_polygons.end());
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

  /**
   * @brief Serializer
   * @param ar Archive for serialization
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int /*version*/)
  {
    ar &boost::serialization::base_object<Geometry>(*this);
    ar & _polygons;
  }

  //-- iterators
  /// @brief Value type for container compatibility
  using value_type = Polygon;
  /// @brief Reference type for container compatibility
  using reference = Polygon &;
  /// @brief Const reference type for container compatibility
  using const_reference = const Polygon &;

  /**
   * @brief Add a polygon patch (container compatibility)
   * @param polygon The polygon to add as a patch
   */
  void
  push_back(const Polygon &polygon)
  {
    addPatch(polygon);
  }

private:
  std::vector<std::unique_ptr<Polygon>> _polygons;

  void
  swap(PolyhedralSurface &other) noexcept
  {
    std::swap(_polygons, other._polygons);
  }
};
} // namespace SFCGAL

#endif
