// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef _SFCGAL_POLYGON_H_
#define _SFCGAL_POLYGON_H_

#include <boost/assert.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/serialize_ptr_vector.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/serialization/base_object.hpp>
#include <vector>

#include <CGAL/Polygon_with_holes_2.h>

#include "SFCGAL/Kernel.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Surface.h"

namespace SFCGAL {

/**
 * A Polygon in SFA with holes
 * @ingroup public_api
 */
class SFCGAL_API Polygon : public Surface {
public:
  typedef boost::ptr_vector<LineString>::iterator       iterator;
  typedef boost::ptr_vector<LineString>::const_iterator const_iterator;

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
  virtual Polygon *
  clone() const;

  //-- SFCGAL::Geometry
  virtual std::string
  geometryType() const;
  //-- SFCGAL::Geometry
  virtual GeometryType
  geometryTypeId() const;
  //-- SFCGAL::Geometry
  virtual int
  coordinateDimension() const;
  //-- SFCGAL::Geometry
  virtual bool
  isEmpty() const;
  //-- SFCGAL::Geometry
  virtual bool
  is3D() const;
  //-- SFCGAL::Geometry
  virtual bool
  isMeasured() const;

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
    return _rings.front();
  }
  /**
   * [OGC/SFA]returns the exterior ring
   */
  inline LineString &
  exteriorRing()
  {
    return _rings.front();
  }
  /**
   * Sets the exterior ring
   */
  inline void
  setExteriorRing(const LineString &ring)
  {
    _rings.front() = ring;
  }
  /**
   * Sets the exterior ring (takes ownership)
   */
  inline void
  setExteriorRing(LineString *ring)
  {
    _rings.replace(0, ring);
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
    return _rings[n + 1];
  }
  /**
   * [OGC/SFA]returns the exterior ring
   */
  inline LineString &
  interiorRingN(const size_t &n)
  {
    return _rings[n + 1];
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
    return _rings[n];
  }
  /**
   * Returns the n-th ring, 0 is exteriorRing
   * @warning not standard, avoid conditionnal to access rings
   */
  inline LineString &
  ringN(const size_t &n)
  {
    BOOST_ASSERT(n < _rings.size());
    return _rings[n];
  }

  /**
   * append a ring to the Polygon
   */
  inline void
  addInteriorRing(const LineString &ls)
  {
    _rings.push_back(ls.clone());
  }
  /**
   * append a ring to the Polygon (take ownership)
   */
  inline void
  addInteriorRing(LineString *ls)
  {
    BOOST_ASSERT(ls != NULL);
    _rings.push_back(ls);
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
    _rings.push_back(ls);
  }
  /**
   * Remove all interior rings
   */
  inline void
  removeInteriorRings()
  {
    if (_rings.size() > 1) {
      _rings.erase(_rings.begin() + 1, _rings.end());
    }
  }

  /**
   * Remove the n-th interior ring
   * @param n index of the interior ring to remove (0-based)
   * @throws std::out_of_range if n is greater than or equal to the number of
   * interior rings
   */
  inline void
  removeInteriorRingN(size_t n)
  {
    if (n >= numInteriorRings()) {
      throw std::out_of_range("Interior ring index out of range");
    }
    _rings.erase(_rings.begin() + n + 1);
  }

  inline iterator
  begin()
  {
    return _rings.begin();
  }
  inline const_iterator
  begin() const
  {
    return _rings.begin();
  }

  inline iterator
  end()
  {
    return _rings.end();
  }
  inline const_iterator
  end() const
  {
    return _rings.end();
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
  virtual void
  accept(GeometryVisitor &visitor);
  //-- SFCGAL::Geometry
  virtual void
  accept(ConstGeometryVisitor &visitor) const;

  /**
   * Serializer
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int /*version*/)
  {
    ar &boost::serialization::base_object<Geometry>(*this);
    ar &_rings;
  }

  /**
   * Returns a const iterator to the first interior ring
   */
  inline const_iterator
  interiorRings_begin() const
  {
    return _rings.size() > 1 ? _rings.begin() + 1 : _rings.end();
  }

  /**
   * Returns a const iterator to the end of interior rings
   */
  inline const_iterator
  interiorRings_end() const
  {
    return _rings.end();
  }

  /**
   * Returns a const range of interior rings
   */
  inline boost::iterator_range<const_iterator>
  interiorRings() const
  {
    return boost::make_iterator_range(interiorRings_begin(),
                                      interiorRings_end());
  }

private:
  /**
   * rings forming the polygon (size() >= 1)
   *
   * _ring[0] is the interior ring
   *
   * @warning never empty, empty LineString as exteriorRing for empty Polygon
   */
  boost::ptr_vector<LineString> _rings;

  void
  swap(Polygon &other)
  {
    std::swap(_rings, other._rings);
  }
};

} // namespace SFCGAL

#endif
