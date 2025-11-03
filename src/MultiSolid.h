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
   * @param other The multi-solid to copy from
   */
  MultiSolid(const MultiSolid &other);
  /**
   * assign operator
   * @param other The multi-solid to assign from
   * @return Reference to this multi-solid
   */
  auto
  operator=(MultiSolid other) -> MultiSolid &;
  /**
   * destructor
   */
  ~MultiSolid() override;

  //-- SFCGAL::Geometry
  /// @brief Get the geometry type as string
  /// @return "MultiSolid"
  [[nodiscard]] auto
  geometryType() const -> std::string override;
  //-- SFCGAL::Geometry
  /// @brief Get the geometry type identifier
  /// @return TYPE_MULTISOLID
  [[nodiscard]] auto
  geometryTypeId() const -> GeometryType override;

  /**
   * returns the n-th Geometry as a Solid
   * @param n The index of the solid to get
   * @return Reference to the nth solid
   */
  inline auto
  solidN(const size_t &n) -> Solid &
  {
    return geometryN(n).as<Solid>();
  }
  /**
   * returns the n-th Geometry as a Solid
   * @param n The index of the solid to get
   * @return Const reference to the nth solid
   */
  inline auto
  solidN(const size_t &n) const -> const Solid &
  {
    return geometryN(n).as<Solid>();
  }

  //-- visitors

  //-- SFCGAL::Geometry
  /// @brief Accept a geometry visitor
  /// @param visitor Visitor to accept
  auto
  accept(GeometryVisitor &visitor) -> void override;
  //-- SFCGAL::Geometry
  /// @brief Accept a const geometry visitor
  /// @param visitor Const visitor to accept
  auto
  accept(ConstGeometryVisitor &visitor) const -> void override;

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
  /// @return true if geometry is a Solid
  auto
  isAllowed(Geometry const &geometry) -> bool override;
};

} // namespace SFCGAL

#endif
