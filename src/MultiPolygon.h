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
class SFCGAL_API MultiPolygon
    : public GeometryImpl<MultiPolygon, GeometryCollection> {
public:
  /**
   * Empty MultiPolygon constructor
   */
  MultiPolygon();
  /**
   * Copy constructor
   * @param other The multi-polygon to copy from
   */
  MultiPolygon(MultiPolygon const &other);
  /**
   * assign operator
   * @param other The multi-polygon to assign from
   * @return Reference to this multi-polygon
   */
  auto
  operator=(MultiPolygon other) -> MultiPolygon &;
  /**
   * destructor
   */
  ~MultiPolygon() override;

  //-- SFCGAL::Geometry
  /// @brief Get the geometry type as string
  /// @return "MultiPolygon"
  [[nodiscard]] auto
  geometryType() const -> std::string override;
  //-- SFCGAL::Geometry
  /// @brief Get the geometry type identifier
  /// @return TYPE_MULTIPOLYGON
  [[nodiscard]] auto
  geometryTypeId() const -> GeometryType override;

  /**
   * returns the n-th Geometry as a Polygon
   * @param n The index of the polygon to get
   * @return Reference to the nth polygon
   */
  inline Polygon &
  polygonN(const size_t &n)
  {
    return geometryN(n).as<Polygon>();
  }
  /**
   * returns the n-th Geometry as a Polygon
   * @param n The index of the polygon to get
   * @return Const reference to the nth polygon
   */
  inline const Polygon &
  polygonN(const size_t &n) const
  {
    return geometryN(n).as<Polygon>();
  }

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
    ar &boost::serialization::base_object<GeometryCollection>(*this);
  }

protected:
  //-- SFCGAL::GeometryCollection
  /// @brief Check if geometry is allowed in this collection
  /// @param geometry Geometry to test
  /// @return true if geometry is a Polygon
  bool
  isAllowed(Geometry const &geometry) override;
};

} // namespace SFCGAL

#endif
