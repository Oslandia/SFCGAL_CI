// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_GEOMETRYCOLLECTION_H_
#define SFCGAL_GEOMETRYCOLLECTION_H_

#include <boost/assert.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/serialize_ptr_vector.hpp>
#include <boost/serialization/base_object.hpp>
#include <vector>

#include "SFCGAL/Geometry.h"

namespace SFCGAL {

/**
 * A GeometryCollection in SFA.
 */
class SFCGAL_API GeometryCollection : public Geometry {
public:
  typedef boost::ptr_vector<Geometry>::iterator       iterator;
  typedef boost::ptr_vector<Geometry>::const_iterator const_iterator;

  /**
   * Empty GeometryCollection constructor
   */
  GeometryCollection();
  /**
   * Copy constructor
   */
  GeometryCollection(const GeometryCollection &other);
  /**
   * assign operator
   */
  GeometryCollection &
  operator=(GeometryCollection other);
  /**
   * destructor
   */
  virtual ~GeometryCollection();

  //-- SFCGAL::Geometry
  GeometryCollection *
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

  //-- SFCGAL::Geometry
  size_t
  numGeometries() const override;
  //-- SFCGAL::Geometry
  const Geometry &
  geometryN(size_t const &n) const override;
  //-- SFCGAL::Geometry
  Geometry &
  geometryN(size_t const &n) override;

  //-- SFCGAL::Geometry
  virtual void
  setGeometryN(const Geometry &geometry, size_t const &n) override;
  //-- SFCGAL::Geometry
  virtual void
  setGeometryN(Geometry *geometry, size_t const &n) override;

  /**
   * [SFA/OGC]add a geometry to the collection (takes ownership)
   */
  void
  addGeometry(Geometry *geometry);
  /**
   * [SFA/OGC]add a geometry to the collection (clone instance)
   */
  void
  addGeometry(Geometry const &geometry);

  //-- iterators

  inline iterator
  begin()
  {
    return _geometries.begin();
  }
  inline const_iterator
  begin() const
  {
    return _geometries.begin();
  }

  inline iterator
  end()
  {
    return _geometries.end();
  }
  inline const_iterator
  end() const
  {
    return _geometries.end();
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
    ar &_geometries;
  }

private:
  boost::ptr_vector<Geometry> _geometries;

protected:
  /**
   * Test if a geometry in the collection
   */
  virtual bool
  isAllowed(Geometry const &g);

  /**
   * Swap
   */
  void
  swap(GeometryCollection &other)
  {
    _geometries.swap(other._geometries);
  }
};

} // namespace SFCGAL

#endif
