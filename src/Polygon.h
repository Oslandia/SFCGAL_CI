// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_POLYGON_H_
#define SFCGAL_POLYGON_H_

#include <boost/assert.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/unique_ptr.hpp>
#include <memory>
#include <vector>

#include <CGAL/Polygon_with_holes_2.h>

#include "SFCGAL/DereferenceIterator.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Surface.h"

namespace SFCGAL {

/**
 * A Polygon in SFA with holes
 */
class SFCGAL_API Polygon : public GeometryImpl<Polygon, Surface> {
public:
  /// @brief Iterator type for polygon rings
  using iterator =
      DereferenceIterator<std::vector<std::unique_ptr<LineString>>::iterator>;
  /// @brief Const iterator type for polygon rings
  using const_iterator = DereferenceIterator<
      std::vector<std::unique_ptr<LineString>>::const_iterator>;

  /**
   * Empty Polygon constructor
   */
  Polygon();
  /**
   * Constructor with an exterior ring
   * @param rings Vector of rings (first is exterior, rest are holes)
   */
  Polygon(const std::vector<LineString> &rings);
  /**
   * Constructor with an exterior ring
   * @param exteriorRing The exterior ring of the polygon
   */
  Polygon(const LineString &exteriorRing);
  /**
   * Constructor with an exterior ring (takes ownership)
   * @param exteriorRing The exterior ring of the polygon
   */
  Polygon(LineString *exteriorRing);
  /**
   * Constructor with a Triangle
   * @param triangle The triangle to convert to a polygon
   */
  Polygon(const Triangle &triangle);
  /**
   * Copy constructor
   * @param other The polygon to copy from
   */
  Polygon(const Polygon &other);

  /**
   * Constructor from CGAL::Polygon_2<K>
   * @param other The CGAL polygon to convert
   */
  Polygon(const CGAL::Polygon_2<Kernel> &other);
  /**
   * Constructor from CGAL::Polygon_with_holes_2<K>
   * @param poly The CGAL polygon with holes to convert
   */
  Polygon(const CGAL::Polygon_with_holes_2<Kernel> &poly);

  /**
   * assign operator
   * @param other The polygon to assign from
   * @return Reference to this polygon
   */
  auto
  operator=(Polygon other) -> Polygon &;

  /**
   * destructor
   */
  ~Polygon() override;

  //-- SFCGAL::Geometry
  /// @brief Get the geometry type as string
  /// @return "Polygon"
  [[nodiscard]] auto
  geometryType() const -> std::string override;
  //-- SFCGAL::Geometry
  /// @brief Get the geometry type identifier
  /// @return TYPE_POLYGON
  [[nodiscard]] auto
  geometryTypeId() const -> GeometryType override;
  //-- SFCGAL::Geometry
  /// @brief Get the coordinate dimension
  /// @return Number of coordinates per point
  [[nodiscard]] auto
  coordinateDimension() const -> int override;
  //-- SFCGAL::Geometry
  /// @brief Check if the polygon is empty
  /// @return true if empty, false otherwise
  [[nodiscard]] auto
  isEmpty() const -> bool override;
  //-- SFCGAL::Geometry
  /// @brief Check if the polygon has 3D coordinates
  /// @return true if 3D, false otherwise
  [[nodiscard]] auto
  is3D() const -> bool override;
  //-- SFCGAL::Geometry
  /// @brief Check if the polygon has measured coordinates
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
   * Check whether the 2D polygon is pointing up
   * @return true if counter-clockwise oriented, false otherwise
   */
  [[nodiscard]] auto
  isCounterClockWiseOriented() const -> bool;

  /**
   * reverse Polygon orientation
   */
  void
  reverse();

  /**
   * [OGC/SFA]returns the exterior ring
   * @return Const reference to the exterior ring
   */
  [[nodiscard]] auto
  exteriorRing() const -> const LineString &
  {
    return *_rings.front();
  }
  /**
   * [OGC/SFA]returns the exterior ring
   * @return Reference to the exterior ring
   */
  auto
  exteriorRing() -> LineString &
  {
    return *_rings.front();
  }
  /**
   * Sets the exterior ring
   * @param ring The new exterior ring
   */
  void
  setExteriorRing(const LineString &ring)
  {
    setExteriorRing(std::unique_ptr<LineString>(ring.clone()));
  }
  /**
   * Sets the exterior ring (takes ownership)
   * @param ring The new exterior ring
   * @deprecated The unique_ptr version should be used instead
   */
  void
  setExteriorRing(LineString *ring)
  {
    setExteriorRing(std::unique_ptr<LineString>(ring));
  }
  /**
   * @brief Sets the exterior ring.
   *
   * Replaces the exterior ring with the provided one.
   *
   * @param ring A unique pointer to the ring object to set. Ownership is
   * transferred to this class.
   */
  void
  setExteriorRing(std::unique_ptr<LineString> ring)
  {
    _rings[0] = std::move(ring);
  }

  /**
   * Test if the polygon has interior rings
   * @return true if polygon has holes, false otherwise
   */
  [[nodiscard]] auto
  hasInteriorRings() const -> bool
  {
    return _rings.size() > 1;
  }

  /**
   * [OGC/SFA]returns the number of interior rings
   * @return Number of interior rings (holes)
   */
  [[nodiscard]] auto
  numInteriorRings() const -> size_t
  {
    return _rings.size() - 1;
  }
  /**
   * [OGC/SFA]returns the nth interior ring
   * @param n The index of the interior ring to get
   * @return Const reference to the nth interior ring
   */
  [[nodiscard]] auto
  interiorRingN(const size_t &n) const -> const LineString &
  {
    return *_rings[n + 1];
  }
  /**
   * [OGC/SFA]returns the nth interior ring
   * @param n The index of the interior ring to get
   * @return Reference to the nth interior ring
   */
  auto
  interiorRingN(const size_t &n) -> LineString &
  {
    return *_rings[n + 1];
  }

  /**
   * Returns the number of rings
   * @return Total number of rings (exterior + interior)
   */
  [[nodiscard]] auto
  numRings() const -> size_t
  {
    return _rings.size();
  }
  /**
   * Returns the n-th ring, 0 is exteriorRing
   * @param n The ring index (0 for exterior ring)
   * @return Const reference to the nth ring
   * @warning not standard, avoid conditionnal to access rings
   */
  [[nodiscard]] auto
  ringN(const size_t &n) const -> const LineString &
  {
    BOOST_ASSERT(n < _rings.size());
    return *_rings[n];
  }
  /**
   * Returns the n-th ring, 0 is exteriorRing
   * @param n The ring index (0 for exterior ring)
   * @return Reference to the nth ring
   * @warning not standard, avoid conditionnal to access rings
   */
  auto
  ringN(const size_t &n) -> LineString &
  {
    BOOST_ASSERT(n < _rings.size());
    return *_rings[n];
  }

  /**
   * append a ring to the Polygon
   * @param ls The line string to add as interior ring
   */
  void
  addInteriorRing(const LineString &ls)
  {
    addInteriorRing(ls.clone());
  }
  /**
   * append a ring to the Polygon (take ownership)
   * @param ls The line string to add as interior ring
   */
  void
  addInteriorRing(LineString *ls)
  {
    addInteriorRing(std::unique_ptr<LineString>(ls));
  }
  /**
   * @brief Adds a ring to the Polygon.
   *
   * @param ring A unique pointer to the LineString object representing the new
   * ring to add. Ownership of this ring is moved into the Polygon.
   */
  void
  addInteriorRing(std::unique_ptr<LineString> ring)
  {
    BOOST_ASSERT(ring != nullptr);
    _rings.push_back(std::move(ring));
  }

  /**
   * @brief append a ring to the Polygon
   * @param ls LineString to add as ring
   * @deprecated addInteriorRing
   */
  void
  addRing(const LineString &ls)
  {
    _rings.push_back(ls.clone());
  }
  /**
   * @brief append a ring to the Polygon (take ownership)
   * @param ls Pointer to LineString to add as ring
   * @deprecated addInteriorRing
   */
  void
  addRing(LineString *ls)
  {
    BOOST_ASSERT(ls != NULL);
    _rings.push_back(std::unique_ptr<LineString>(ls));
  }

  /**
   * @brief Get iterator to beginning of rings
   * @return Iterator to first ring
   */
  auto
  begin() -> iterator
  {
    return dereference_iterator(_rings.begin());
  }
  /**
   * @brief Get const iterator to beginning of rings
   * @return Const iterator to first ring
   */
  [[nodiscard]] auto
  begin() const -> const_iterator
  {
    return dereference_iterator(_rings.begin());
  }

  /**
   * @brief Get iterator to end of rings
   * @return Iterator to past-the-end
   */
  auto
  end() -> iterator
  {
    return dereference_iterator(_rings.end());
  }
  /**
   * @brief Get const iterator to end of rings
   * @return Const iterator to past-the-end
   */
  [[nodiscard]] auto
  end() const -> const_iterator
  {
    return dereference_iterator(_rings.end());
  }

  /**
   * @brief Convert to CGAL::Polygon_2. Does not consider holes, if any
   * @param fixOrientation force exterior ring orientation to counter clockwise
   * @return CGAL 2D polygon representation
   */
  [[nodiscard]] auto
  toPolygon_2(bool fixOrientation = true) const -> CGAL::Polygon_2<Kernel>;

  /**
   * @brief Convert to CGAL::Polygon_with_holes_2.
   * @param fixOrientation force exterior ring orientation to counter clockwise
   * and interior ring to clockwise.
   * @return CGAL 2D polygon with holes representation
   */
  [[nodiscard]] auto
  toPolygon_with_holes_2(bool fixOrientation = true) const
      -> CGAL::Polygon_with_holes_2<Kernel>;

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
    ar & _rings;
  }

private:
  /**
   * rings forming the polygon (size() >= 1)
   *
   * _ring[0] is the interior ring
   *
   * @warning never empty, empty LineString as exteriorRing for empty Polygon
   */
  std::vector<std::unique_ptr<LineString>> _rings;

  void
  swap(Polygon &other) noexcept
  {
    std::swap(_rings, other._rings);
  }
};

} // namespace SFCGAL

#endif
