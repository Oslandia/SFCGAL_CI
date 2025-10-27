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
class SFCGAL_API GeometryCollection
    : public GeometryImpl<GeometryCollection, Geometry> {
public:
  /// @brief Iterator type for geometry collection
  using iterator =
      DereferenceIterator<std::vector<std::unique_ptr<Geometry>>::iterator>;
  /// @brief Const iterator type for geometry collection
  using const_iterator = DereferenceIterator<
      std::vector<std::unique_ptr<Geometry>>::const_iterator>;

  /**
   * Empty GeometryCollection constructor
   */
  GeometryCollection();
  /**
   * Copy constructor
   * @param other The geometry collection to copy from
   */
  GeometryCollection(const GeometryCollection &other);
  /**
   * assign operator
   * @param other The geometry collection to assign from
   * @return Reference to this geometry collection
   */
  auto
  operator=(GeometryCollection other) -> GeometryCollection &;
  /**
   * destructor
   */
  ~GeometryCollection() override;

  //-- SFCGAL::Geometry
  /// @brief Get the geometry type as string
  /// @return "GeometryCollection"
  [[nodiscard]] auto
  geometryType() const -> std::string override;
  //-- SFCGAL::Geometry
  /// @brief Get the geometry type identifier
  /// @return TYPE_GEOMETRYCOLLECTION
  [[nodiscard]] auto
  geometryTypeId() const -> GeometryType override;
  //-- SFCGAL::Geometry
  /// @brief Get the dimension of the collection
  /// @return Maximum dimension of contained geometries
  [[nodiscard]] auto
  dimension() const -> int override;
  //-- SFCGAL::Geometry
  /// @brief Get the coordinate dimension
  /// @return Number of coordinates per point
  [[nodiscard]] auto
  coordinateDimension() const -> int override;
  //-- SFCGAL::Geometry
  /// @brief Check if the collection is empty
  /// @return true if empty, false otherwise
  [[nodiscard]] auto
  isEmpty() const -> bool override;
  //-- SFCGAL::Geometry
  /// @brief Check if the collection has 3D coordinates
  /// @return true if any geometry is 3D, false otherwise
  [[nodiscard]] auto
  is3D() const -> bool override;
  //-- SFCGAL::Geometry
  /// @brief Check if the collection has measured coordinates
  /// @return true if any geometry is measured, false otherwise
  [[nodiscard]] auto
  isMeasured() const -> bool override;

  /// @brief Drop Z coordinate from all geometries
  /// @return true if Z was dropped, false otherwise
  auto
  dropZ() -> bool override;

  /// @brief Drop M coordinate from all geometries
  /// @return true if M was dropped, false otherwise
  auto
  dropM() -> bool override;

  /// @brief Swap X and Y coordinates of all geometries in collection
  auto
  swapXY() -> void override;

  //-- SFCGAL::Geometry
  /// @brief Get the number of geometries in the collection
  /// @return Number of geometries
  [[nodiscard]] auto
  numGeometries() const -> size_t override;
  //-- SFCGAL::Geometry
  /// @brief Get the n-th geometry in the collection (const version)
  /// @param n Index of the geometry to retrieve
  /// @return Const reference to the geometry
  auto
  geometryN(size_t const &n) const -> const Geometry & override;
  //-- SFCGAL::Geometry
  /// @brief Get the n-th geometry in the collection (non-const version)
  /// @param n Index of the geometry to retrieve
  /// @return Reference to the geometry
  auto
  geometryN(size_t const &n) -> Geometry & override;

  //-- SFCGAL::Geometry
  /// @brief Set the geometry at the given index (copy version)
  /// @param geometry Geometry to copy and set
  /// @param idx Index where to set the geometry
  void
  setGeometryN(const Geometry &geometry, size_t const &idx) override;
  //-- SFCGAL::Geometry
  /// @brief Set the geometry at the given index (pointer version)
  /// @param geometry Pointer to geometry to set
  /// @param idx Index where to set the geometry
  void
  setGeometryN(Geometry *geometry, size_t const &idx) override;
  /**
   * @brief [SFA/OGC] Sets the n-th geometry using unique_ptr
   * @param geometry Unique pointer to geometry to set
   * @param idx Index where to set the geometry
   */
  void
  setGeometryN(std::unique_ptr<Geometry> geometry, size_t const &idx) override;

  /**
   * @brief [SFA/OGC] Add a geometry to the collection.
   *
   * @param geometry A unique pointer to the Geometry object to add. Ownership
   * of this geometry is moved into the GeometryCollection.
   */
  void
  addGeometry(std::unique_ptr<Geometry> geometry);
  /**
   * @brief [SFA/OGC]add a geometry to the collection (takes ownership)
   * @param geometry Pointer to geometry to add
   * @deprecated The unique_ptr version should be used instead
   */
  void
  addGeometry(Geometry *geometry);
  /**
   * @brief [SFA/OGC]add a geometry to the collection (clone instance)
   * @param geometry Const reference to geometry to clone and add
   */
  void
  addGeometry(Geometry const &geometry);

  //-- iterators

  /**
   * @brief Get iterator to beginning of geometries
   * @return Iterator to first geometry
   */
  inline auto
  begin() -> iterator
  {
    return dereference_iterator(_geometries.begin());
  }
  /**
   * @brief Get const iterator to beginning of geometries
   * @return Const iterator to first geometry
   */
  inline auto
  begin() const -> const_iterator
  {
    return dereference_iterator(_geometries.begin());
  }

  /**
   * @brief Get iterator to end of geometries
   * @return Iterator to past-the-end
   */
  inline auto
  end() -> iterator
  {
    return dereference_iterator(_geometries.end());
  }
  /**
   * @brief Get const iterator to end of geometries
   * @return Const iterator to past-the-end
   */
  inline auto
  end() const -> const_iterator
  {
    return dereference_iterator(_geometries.end());
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
    ar &boost::serialization::base_object<Geometry>(*this);
    ar & _geometries;
  }

private:
  std::vector<std::unique_ptr<Geometry>> _geometries;

protected:
  /**
   * @brief Test if a geometry is allowed in the collection
   * @param geometry Geometry to test
   * @return true if geometry is allowed in this collection
   */
  virtual auto
  isAllowed(Geometry const &geometry) -> bool;

  /**
   * @brief Swap the contents of two geometry collections
   * @param other GeometryCollection to swap with
   */
  void
  swap(GeometryCollection &other)
  {
    _geometries.swap(other._geometries);
  }
};

} // namespace SFCGAL

#endif
