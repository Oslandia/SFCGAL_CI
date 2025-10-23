// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_GEOMETRYCOLLECTION_H_
#define SFCGAL_GEOMETRYCOLLECTION_H_

#include <boost/assert.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/unique_ptr.hpp>
#include <memory>
#include <vector>

#include "SFCGAL/DereferenceIterator.h"
#include "SFCGAL/Geometry.h"

namespace SFCGAL {

/**
 * A GeometryCollection in SFA.
 */
class SFCGAL_API GeometryCollection : public Geometry {
public:
  using iterator =
      DereferenceIterator<std::vector<std::unique_ptr<Geometry>>::iterator>;
  using const_iterator = DereferenceIterator<
      std::vector<std::unique_ptr<Geometry>>::const_iterator>;

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

  auto
  swapXY() -> void override;

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
    return dereference_iterator(_geometries.begin());
  }
  inline const_iterator
  begin() const
  {
    return dereference_iterator(_geometries.begin());
  }

  inline iterator
  end()
  {
    return dereference_iterator(_geometries.end());
  }
  inline const_iterator
  end() const
  {
    return dereference_iterator(_geometries.end());
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
    ar & _geometries;
  }

private:
  std::vector<std::unique_ptr<Geometry>> _geometries;

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
