// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_MULTISOLID_H_
#define SFCGAL_MULTISOLID_H_

#include <boost/assert.hpp>
#include <vector>

#include <boost/serialization/base_object.hpp>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Solid.h"

namespace SFCGAL {

/**
 * A MultiSolid
 */
class SFCGAL_API MultiSolid
    : public GeometryImpl<MultiSolid, GeometryCollection> {
public:
  /**
   * Empty MultiSolid constructor
   */
  MultiSolid();
  /**
   * Copy constructor
   */
  MultiSolid(const MultiSolid &other);
  /**
   * assign operator
   */
  MultiSolid &
  operator=(MultiSolid other);
  /**
   * destructor
   */
  virtual ~MultiSolid();

  //-- SFCGAL::Geometry
  std::string
  geometryType() const override;
  //-- SFCGAL::Geometry
  GeometryType
  geometryTypeId() const override;

  /**
   * returns the n-th Geometry as a Solid
   */
  inline Solid &
  solidN(const size_t &n)
  {
    return geometryN(n).as<Solid>();
  }
  /**
   * returns the n-th Geometry as a Solid
   */
  inline const Solid &
  solidN(const size_t &n) const
  {
    return geometryN(n).as<Solid>();
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
