// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
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
class SFCGAL_API MultiLineString
    : public GeometryImpl<MultiLineString, GeometryCollection> {
public:
  /**
   * Empty MultiLineString constructor
   */
  MultiLineString();
  /**
   * Copy constructor
   * @param other The multi-linestring to copy from
   */
  MultiLineString(const MultiLineString &other);
  /**
   * assign operator
   * @param other The multi-linestring to assign from
   * @return Reference to this multi-linestring
   */
  MultiLineString &
  operator=(MultiLineString other);
  /**
   * destructor
   */
  ~MultiLineString() override;

  //-- SFCGAL::Geometry
  /// @brief Get the geometry type as string
  /// @return "MultiLineString"
  [[nodiscard]] auto
  geometryType() const -> std::string override;
  //-- SFCGAL::Geometry
  /// @brief Get the geometry type identifier
  /// @return TYPE_MULTILINESTRING
  [[nodiscard]] auto
  geometryTypeId() const -> GeometryType override;

  /**
   * returns the n-th Geometry as a LineString
   * @param n The index of the linestring to get
   * @return Reference to the nth linestring
   */
  inline auto
  lineStringN(const size_t &n) -> LineString &
  {
    return geometryN(n).as<LineString>();
  }
  /**
   * returns the n-th Geometry as a LineString
   * @param n The index of the linestring to get
   * @return Const reference to the nth linestring
   */
  inline auto
  lineStringN(const size_t &n) const -> const LineString &
  {
    return geometryN(n).as<LineString>();
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
  /// @return true if geometry is a LineString
  auto
  isAllowed(Geometry const &geometry) -> bool override;
};

} // namespace SFCGAL

#endif
