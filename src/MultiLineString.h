// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_MULTILINESTRING_H_
#define SFCGAL_MULTILINESTRING_H_

#include <boost/assert.hpp>
#include <vector>

#include <boost/serialization/base_object.hpp>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"

namespace SFCGAL {

/**
 * A MultiLineString in SFA.
 */
class SFCGAL_API MultiLineString : public GeometryCollection {
public:
  /**
   * Empty MultiLineString constructor
   */
  MultiLineString();
  /**
   * Copy constructor
   */
  MultiLineString(const MultiLineString &other);
  /**
   * assign operator
   */
  MultiLineString &
  operator=(MultiLineString other);
  /**
   * destructor
   */
  virtual ~MultiLineString();

  //-- SFCGAL::Geometry
  virtual MultiLineString *
  clone() const override;

  //-- SFCGAL::Geometry
  virtual std::string
  geometryType() const override;
  //-- SFCGAL::Geometry
  virtual GeometryType
  geometryTypeId() const override;

  /**
   * returns the n-th Geometry as a Polygon
   */
  inline LineString &
  lineStringN(const size_t &n)
  {
    return geometryN(n).as<LineString>();
  }
  /**
   * returns the n-th Geometry as a Polygon
   */
  inline const LineString &
  lineStringN(const size_t &n) const
  {
    return geometryN(n).as<LineString>();
  }

  //-- visitors

  //-- SFCGAL::Geometry
  virtual void
  accept(GeometryVisitor &visitor) override;
  //-- SFCGAL::Geometry
  virtual void
  accept(ConstGeometryVisitor &visitor) const override;

  /**
   * Serializer
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int /*version*/)
  {
    ar &boost::serialization::base_object<GeometryCollection>(*this);
  }

protected:
  //-- SFCGAL::GeometryCollection
  virtual bool
  isAllowed(Geometry const &g) override;
};

} // namespace SFCGAL

#endif
