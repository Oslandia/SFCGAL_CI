// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_MULTIPOINT_H_
#define SFCGAL_MULTIPOINT_H_

#include <boost/assert.hpp>
#include <vector>

#include <boost/serialization/base_object.hpp>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Point.h"

namespace SFCGAL {

/**
 * A MultiPoint in SFA.
 */
class SFCGAL_API MultiPoint
    : public GeometryImpl<MultiPoint, GeometryCollection> {
public:
  /**
   * Empty MultiPoint constructor
   */
  MultiPoint();
  /**
   * Copy constructor
   * @param other The multi-point to copy from
   */
  MultiPoint(const MultiPoint &other);
  /**
   * assign operator
   * @param other The multi-point to assign from
   * @return Reference to this multi-point
   */
  auto
  operator=(MultiPoint other) -> MultiPoint &;
  /**
   * destructor
   */
  ~MultiPoint() override;

  //-- SFCGAL::Geometry
  /// @brief Get the geometry type as string
  /// @return "MultiPoint"
  [[nodiscard]] auto
  geometryType() const -> std::string override;
  //-- SFCGAL::Geometry
  /// @brief Get the geometry type identifier
  /// @return TYPE_MULTIPOINT
  [[nodiscard]] auto
  geometryTypeId() const -> GeometryType override;

  /**
   * returns the n-th Geometry as a Point
   * @param n The index of the point to get
   * @return Reference to the nth point
   */
  inline auto
  pointN(const size_t &n) -> Point &
  {
    return geometryN(n).as<Point>();
  }
  /**
   * returns the n-th Geometry as a Point
   * @param n The index of the point to get
   * @return Const reference to the nth point
   */
  inline auto
  pointN(const size_t &n) const -> const Point &
  {
    return geometryN(n).as<Point>();
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
  /// @return true if geometry is a Point
  auto
  isAllowed(Geometry const &geometry) -> bool override;
};

} // namespace SFCGAL

#endif
