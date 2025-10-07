// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
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

  /**
   * @copydoc SFCGAL::Geometry::clone()
   */
  MultiPolygon *
  clone() const override;

  /**
   * @copydoc SFCGAL::Geometry::geometryType()
   */
  std::string
  geometryType() const override;
  /**
   * @copydoc SFCGAL::Geometry::geometryTypeId()
   */
  GeometryType
  geometryTypeId() const override;

  /**
   * returns the n-th Geometry as a Polygon
   * @param n index of the polygon
   */
  inline Polygon &
  polygonN(const size_t &n)
  {
    return geometryN(n).as<Polygon>();
  }
  /**
   * returns the n-th Geometry as a Polygon
   * @param n index of the polygon
   */
  inline const Polygon &
  polygonN(const size_t &n) const
  {
    return geometryN(n).as<Polygon>();
  }

  //-- visitors

  /**
   * @copydoc SFCGAL::Geometry::accept()
   */
  void
  accept(GeometryVisitor &visitor) override;
  /**
   * @copydoc SFCGAL::Geometry::accept()
   */
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
  /**
   * @copydoc SFCGAL::GeometryCollection::isAllowed()
   */
  bool
  isAllowed(Geometry const &g) override;
};

} // namespace SFCGAL

#endif
