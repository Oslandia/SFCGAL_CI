// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_REGISTRY_H_
#define SFCGAL_REGISTRY_H_

#include "SFCGAL/config.h"

#include <map>
#include <string>
#include <vector>

namespace SFCGAL {
class Geometry;
}

namespace SFCGAL {
namespace tools {

/**
 * Registry for dynamic information about SFCGAL library
 */
class SFCGAL_API Registry {
public:
  typedef std::vector<Geometry *>::iterator       prototype_iterator; ///< Iterator type for prototypes
  typedef std::vector<Geometry *>::const_iterator const_prototype_iterator; ///< Const iterator type for prototypes

  /**
   * destructor
   */
  ~Registry();

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
   * @return Pointer to new geometry instance
   */
  Geometry *
  newGeometryByTypeName(const std::string &geometryTypeName) const;

  /**
   * @brief Returns a new instance of the given geometryType
   * @param typeId The geometry type ID
   * @return Pointer to new geometry instance
   */
  Geometry *
  newGeometryByTypeId(int typeId) const;

  /**
   * @brief Returns the instance of the registry
   * @return Reference to the singleton registry instance
   */
  static Registry &
  instance();

private:
  /**
   * static instance of the Singleton
   */
  static Registry *_instance;
  /**
   * prototypes of the geometries
   */
  std::vector<Geometry *> _prototypes;

  /**
   * init registry
   */
  Registry();
};

} // namespace tools
} // namespace SFCGAL

#endif
