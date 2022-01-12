/**
 *   SFCGAL
 *
 *   Copyright (C) 2012-2013 Oslandia <infos@oslandia.com>
 *   Copyright (C) 2012-2013 IGN (http://www.ign.fr)
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Library General Public
 *   License as published by the Free Software Foundation; either
 *   version 2 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Library General Public License for more details.

 *   You should have received a copy of the GNU Library General Public
 *   License along with this library; if not, see
 <http://www.gnu.org/licenses/>.
 */

#ifndef _SFCGAL_REGISTRY_H_
#define _SFCGAL_REGISTRY_H_

#include <SFCGAL/config.h>

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
  typedef std::vector<Geometry *>::iterator       prototype_iterator;
  typedef std::vector<Geometry *>::const_iterator const_prototype_iterator;

  /**
   * destructor
   */
  ~Registry();

  /**
   * Register a new Geometry type
   */
  void
  addPrototype(const Geometry &g);

  /**
   * returns the list of the geometry types
   */
  std::vector<std::string>
  getGeometryTypes() const;

  /**
   * returns a new instance of the given geometryTypeName
   */
  Geometry *
  newGeometryByTypeName(const std::string &geometryTypeName) const;

  /**
   * returns a new instance of the given geometryType
   */
  Geometry *
  newGeometryByTypeId(int typeId) const;

  /**
   * returns the instance of the registry
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
