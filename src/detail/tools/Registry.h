// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_REGISTRY_H_
#define SFCGAL_REGISTRY_H_

#include "SFCGAL/DereferenceIterator.h"
#include "SFCGAL/config.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace SFCGAL {
class Geometry;
}

namespace SFCGAL::tools {

/**
 * Registry for dynamic information about SFCGAL library
 */
class SFCGAL_API Registry {
public:
  using prototype_iterator =
      DereferenceIterator<std::vector<std::unique_ptr<Geometry>>>::
          iterator; ///< Iterator type for prototypes
  using const_prototype_iterator =
      DereferenceIterator<std::vector<std::unique_ptr<Geometry>>>::
          const_iterator; ///< Const iterator type for prototypes

  /**
   * destructor
   */
  ~Registry() = default;

  // Delete copy operations
  Registry(const Registry &) = delete;
  auto
  operator=(const Registry &) -> Registry & = delete;

  /**
   * @brief Register a new Geometry type
   * @param g The geometry prototype to register
   */
  void
  addPrototype(const Geometry &g);

  /**
   * @brief Returns the list of the geometry types
   * @return Vector of geometry type names
   */
  std::vector<std::string>
  getGeometryTypes() const;

  /**
   * @brief Returns a new instance of the given geometryTypeName
   * @param geometryTypeName The name of the geometry type
   * @return std::unique_ptr to new geometry instance
   */
  [[nodiscard]] auto
  newGeometryByTypeName(const std::string &geometryTypeName) const
      -> std::unique_ptr<Geometry>;

  /**
   * @brief Returns a new instance of the given geometryType
   * @param typeId The geometry type ID
   * @return Pointer to new geometry instance
   */
  [[nodiscard]] auto
  newGeometryByTypeId(int typeId) const -> std::unique_ptr<Geometry>;

  /**
   * @brief Returns the instance of the registry
   * @return Reference to the singleton registry instance
   */
  static Registry &
  instance();

  //-- iterators

  /**
   * @brief Returns a mutable iterator to the beginning of the registry
   * collection.
   *
   * @return prototype_iterator Iterator to the first element.
   */
  auto
  begin() -> prototype_iterator
  {
    return dereference_iterator(_prototypes.begin());
  }

  /**
   * @brief Returns a const iterator to the beginning of the registry
   * collection.
   *
   * @return const_prototype_iterator Const iterator to the first element.
   */
  [[nodiscard]] auto
  begin() const -> const_prototype_iterator
  {
    return dereference_iterator(_prototypes.begin());
  }

  /**
   * @brief Returns a mutable iterator to the end of the registry collection.
   *
   * @return prototype_iterator Iterator pointing past the last element.
   */
  auto
  end() -> prototype_iterator
  {
    return dereference_iterator(_prototypes.end());
  }

  /**
   * @brief Returns a const iterator to the end of the registry collection.
   *
   * @return const_prototype_iterator Const iterator pointing past the last
   * element.
   */
  [[nodiscard]] auto
  end() const -> const_prototype_iterator
  {
    return dereference_iterator(_prototypes.end());
  }

private:
  /**
   * static instance of the Singleton
   */
  static Registry *_instance;
  /**
   * prototypes of the geometries
   */
  std::vector<std::unique_ptr<Geometry>> _prototypes;

  /**
   * init registry
   */
  Registry();
};

} // namespace SFCGAL::tools

#endif
