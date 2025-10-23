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
class SFCGAL_API PolyhedralSurface : public Surface {
public:
  using iterator =
      DereferenceIterator<std::vector<std::unique_ptr<Polygon>>::iterator>;
  using const_iterator = DereferenceIterator<
      std::vector<std::unique_ptr<Polygon>>::const_iterator>;

  /**
   * Empty PolyhedralSurface constructor
   */
  PolyhedralSurface();
  /**
   * Constructor with a vector of polygons
   */
  PolyhedralSurface(const std::vector<Polygon> &polygons);

  /**
   * Constructor with a Geometry
   */
  PolyhedralSurface(const std::unique_ptr<Geometry> &geometry);

  /**
   * Constructor from a Polyhedron (detail::MarkedPolyhedron or
   * CGAL::Polyhedron_3)
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
   */
  PolyhedralSurface(const Mesh &sm);

  PolyhedralSurface(const InexactMesh &inexactMesh);
  /**
   * Copy constructor
   */
  PolyhedralSurface(const PolyhedralSurface &other);
  /**
   * assign operator
   */
  PolyhedralSurface &
  operator=(PolyhedralSurface other);
  /**
   * destructor
   */
  ~PolyhedralSurface();

  //-- SFCGAL::Geometry
  PolyhedralSurface *
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
   * Convert PolyhedralSurface to TriangulatedSurface
   */
  TriangulatedSurface
  toTriangulatedSurface() const;

  /**
   * [SFA/OGC]Returns the number of patches
   * @warning PolyhedralSurface is treated as one geometry. numGeometries
   * returns 1 or 0 for empty PolyhedralSurface
   */
  inline auto
  numPatches() const -> size_t
  {
    return _polygons.size();
  }

  /**
   * [SFA/OGC]Returns the number of polygons
   *
   * @warning PolyhedralSurface is treated as one geometry. numGeometries
   * returns 1 or 0 for empty PolyhedralSurface
   * @deprecated see numPatches
   * @see numGeometries()
   */
  inline size_t
  numPolygons() const
  {
    return numPatches();
  }

  /**
   * [SFA/OGC]Returns the n-th patch
   */
  inline const Polygon &
  patchN(size_t const &n) const
  {
    BOOST_ASSERT(n < _polygons.size());
    return *_polygons[n];
  }
  /**
   * [SFA/OGC]Returns the n-th patch
   */
  inline Polygon &
  patchN(size_t const &n)
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
   */
  void
  addPatchs(const PolyhedralSurface &polyhedralSurface);

  /**
   * [SFA/OGC]Returns the n-th polygon
   * @deprecated see patchN()
   */
  inline const Polygon &
  polygonN(size_t const &n) const
  {
    return patchN(n);
  }
  /**
   * [SFA/OGC]Returns the n-th polygon
   * @deprecated see patchN()
   */
  inline Polygon &
  polygonN(size_t const &n)
  {
    return patchN(n);
  }
  /**
   * add a polygon to the PolyhedralSurface
   * @deprecated see addPatch()
   */
  void
  addPolygon(const Polygon &polygon);
  /**
   * add a polygon to the PolyhedralSurface
   * @deprecated see addPatch()
   */
  void
  addPolygon(Polygon *polygon);
  /**
   * add polygons from an other PolyhedralSurface
   * @deprecated see addPatchs()
   */
  void
  addPolygons(const PolyhedralSurface &polyhedralSurface);

  /**
   * Sets the n-th Geometry, starting at zero
   * It needs to be a polygon.
   */
  void
  setPatchN(const Geometry &geometry, size_t const &n);

  /**
   * Sets the n-th Geometry, starting at zero
   * It needs to be a polygon.
   * The ownership of the polygon is taken. The caller is not responsible
   * anymore of its deallocation.
   */
  void
  setPatchN(Geometry *geometry, size_t const &n);

  /**
   * Sets the n-th Patch, starting at zero
   */
  void
  setPatchN(const Polygon &patch, size_t const &n);

  /**
   * Sets the n-th Patch, starting at zero
   * The ownership of the polygon is taken. The caller is not responsible
   * anymore of its deallocation.
   */
  void
  setPatchN(Polygon *patch, size_t const &n);

  /**
   * Convert to CGAL::Polyhedron_3
   */
  template <typename K, typename Polyhedron>
  std::unique_ptr<Polyhedron>
  toPolyhedron_3() const
  {
    TriangulatedSurface tri;
    triangulate::triangulatePolygon3D(*this, tri);
    return tri.toPolyhedron_3<K, Polyhedron>();
  }

  //-- iterators

  inline iterator
  begin()
  {
    return dereference_iterator(_polygons.begin());
  }
  inline const_iterator
  begin() const
  {
    return dereference_iterator(_polygons.begin());
  }

  inline iterator
  end()
  {
    return dereference_iterator(_polygons.end());
  }
  inline const_iterator
  end() const
  {
    return dereference_iterator(_polygons.end());
  }

  //-- visitors

  //-- SFCGAL::Geometry
  void
  accept(GeometryVisitor &visitor) override;
  //-- SFCGAL::Geometry
  void
  accept(ConstGeometryVisitor &visitor) const override;

  /**
   * Serializer
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int /*version*/)
  {
    ar &boost::serialization::base_object<Geometry>(*this);
    ar & _polygons;
  }

  //-- iterators
  using value_type      = Polygon;
  using reference       = Polygon &;
  using const_reference = const Polygon &;

  void
  push_back(const Polygon &polygon)
  {
    addPatch(polygon);
  }

private:
  std::vector<std::unique_ptr<Polygon>> _polygons;

  void
  swap(PolyhedralSurface &other)
  {
    std::swap(_polygons, other._polygons);
  }
};
} // namespace SFCGAL

#endif
