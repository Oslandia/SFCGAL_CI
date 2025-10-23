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
  using iterator =
      DereferenceIterator<std::vector<std::unique_ptr<LineString>>::iterator>;
  using const_iterator = DereferenceIterator<
      std::vector<std::unique_ptr<LineString>>::const_iterator>;

  /**
   * Empty Polygon constructor
   */
  Polygon();
  /**
   * Constructor with an exterior ring
   */
  Polygon(const std::vector<LineString> &rings);
  /**
   * Constructor with an exterior ring
   */
  Polygon(const LineString &exteriorRing);
  /**
   * Constructor with an exterior ring (takes ownership)
   */
  Polygon(LineString *exteriorRing);
  /**
   * Constructor with a Triangle
   */
  Polygon(const Triangle &triangle);
  /**
   * Copy constructor
   */
  Polygon(const Polygon &other);

  /**
   * Constructor from CGAL::Polygon_with_holes_2<K>
   */
  Polygon(const CGAL::Polygon_2<Kernel> &other);
  /**
   * Constructor from CGAL::Polygon_with_holes_2<K>
   */
  Polygon(const CGAL::Polygon_with_holes_2<Kernel> &other);

  /**
   * assign operator
   */
  Polygon &
  operator=(Polygon other);

  /**
   * destructor
   */
  ~Polygon();

  //-- SFCGAL::Geometry
  std::string
  geometryType() const override;
  //-- SFCGAL::Geometry
  GeometryType
  geometryTypeId() const override;
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
   * Check whether the 2D polygon is pointing up
   */
  bool
  isCounterClockWiseOriented() const;

  /**
   * reverse Polygon orientation
   */
  void
  reverse();

  /**
   * [OGC/SFA]returns the exterior ring
   */
  inline const LineString &
  exteriorRing() const
  {
    return *_rings.front();
  }
  /**
   * [OGC/SFA]returns the exterior ring
   */
  inline LineString &
  exteriorRing()
  {
    return *_rings.front();
  }
  /**
   * Sets the exterior ring
   */
  void
  setExteriorRing(const LineString &ring)
  {
    setExteriorRing(std::unique_ptr<LineString>(ring.clone()));
  }
  /**
   * Sets the exterior ring (takes ownership)
   *
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
   */
  inline bool
  hasInteriorRings() const
  {
    return _rings.size() > 1;
  }

  /**
   * [OGC/SFA]returns the exterior ring
   */
  inline size_t
  numInteriorRings() const
  {
    return _rings.size() - 1;
  }
  /**
   * [OGC/SFA]returns the exterior ring
   */
  inline const LineString &
  interiorRingN(const size_t &n) const
  {
    return *_rings[n + 1];
  }
  /**
   * [OGC/SFA]returns the exterior ring
   */
  inline LineString &
  interiorRingN(const size_t &n)
  {
    return *_rings[n + 1];
  }

  /**
   * Returns the number of rings
   */
  inline size_t
  numRings() const
  {
    return _rings.size();
  }
  /**
   * Returns the n-th ring, 0 is exteriorRing
   * @warning not standard, avoid conditionnal to access rings
   */
  inline const LineString &
  ringN(const size_t &n) const
  {
    BOOST_ASSERT(n < _rings.size());
    return *_rings[n];
  }
  /**
   * Returns the n-th ring, 0 is exteriorRing
   * @warning not standard, avoid conditionnal to access rings
   */
  inline LineString &
  ringN(const size_t &n)
  {
    BOOST_ASSERT(n < _rings.size());
    return *_rings[n];
  }

  /**
   * append a ring to the Polygon
   */
  void
  addInteriorRing(const LineString &ls)
  {
    addInteriorRing(ls.clone());
  }
  /**
   * append a ring to the Polygon (take ownership)
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
   * append a ring to the Polygon
   * @deprecated addInteriorRing
   */
  inline void
  addRing(const LineString &ls)
  {
    _rings.push_back(ls.clone());
  }
  /**
   * append a ring to the Polygon (take ownership)
   * @deprecated addInteriorRing
   */
  inline void
  addRing(LineString *ls)
  {
    BOOST_ASSERT(ls != NULL);
    _rings.push_back(std::unique_ptr<LineString>(ls));
  }

  inline iterator
  begin()
  {
    return dereference_iterator(_rings.begin());
  }
  inline const_iterator
  begin() const
  {
    return dereference_iterator(_rings.begin());
  }

  inline iterator
  end()
  {
    return dereference_iterator(_rings.end());
  }
  inline const_iterator
  end() const
  {
    return dereference_iterator(_rings.end());
  }

  /*
   * @brief Convert to CGAL::Polygon_2. Does not consider holes, if any
   * @param forceCounterClocksize force exterior ring orientation to counter
   * clocksize
   */
  CGAL::Polygon_2<Kernel>
  toPolygon_2(bool fixOrientation = true) const;

  /*
   * @brief Convert to CGAL::Polygon_with_holes_2.
   * @param forceCounterClocksize force exterior ring orientation to counter
   * clocksize and interior ring to clocksize.
   */
  CGAL::Polygon_with_holes_2<Kernel>
  toPolygon_with_holes_2(bool fixOrientation = true) const;

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
  swap(Polygon &other)
  {
    std::swap(_rings, other._rings);
  }
};

} // namespace SFCGAL

#endif
