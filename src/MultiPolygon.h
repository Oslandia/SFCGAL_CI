// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_MULTIPOLYGON_H_
#define SFCGAL_MULTIPOLYGON_H_

#include <boost/assert.hpp>
#include <vector>

#include <boost/serialization/base_object.hpp>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Polygon.h"

namespace SFCGAL {

/**
 * A MultiPolygon in SFA.
 * @ŧodo add polygon() etc.
 */
class SFCGAL_API MultiPolygon : public GeometryCollection {
public:
  /**
   * Empty MultiPolygon constructor
   */
  MultiPolygon();
  /**
   * Copy constructor
   */
  MultiPolygon(MultiPolygon const &other);
  /**
   * assign operator
   */
  MultiPolygon &
  operator=(MultiPolygon other);
  /**
   * destructor
   */
  virtual ~MultiPolygon();

  //-- SFCGAL::Geometry
  MultiPolygon *
  clone() const override;

  //-- SFCGAL::Geometry
  std::string
  geometryType() const override;
  //-- SFCGAL::Geometry
  GeometryType
  geometryTypeId() const override;

  /**
   * returns the n-th Geometry as a Polygon
   */
  inline Polygon &
  polygonN(const size_t &n)
  {
    return geometryN(n).as<Polygon>();
  }
  /**
   * returns the n-th Geometry as a Polygon
   */
  inline const Polygon &
  polygonN(const size_t &n) const
  {
    return geometryN(n).as<Polygon>();
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
    ar &boost::serialization::base_object<GeometryCollection>(*this);
  }

protected:
  //-- SFCGAL::GeometryCollection
  bool
  isAllowed(Geometry const &g) override;
};

} // namespace SFCGAL

#endif
