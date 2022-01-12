// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/detail/tools/Registry.h>

#include <SFCGAL/Geometry.h>
#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/MultiSolid.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>

#include <SFCGAL/Exception.h>
#include <SFCGAL/detail/tools/Log.h>

namespace SFCGAL {
namespace tools {

Registry *Registry::_instance = 0;

///
///
///
Registry::~Registry()
{
  for (prototype_iterator it = _prototypes.begin(); it != _prototypes.end();
       ++it) {
    delete *it;
  }
}

///
///
///
void
Registry::addPrototype(const Geometry &g)
{
  // find prototype by name
  const_prototype_iterator it = _prototypes.begin();

  for (; it != _prototypes.end(); ++it) {
    if ((*it)->geometryTypeId() == g.geometryTypeId()) {
      break;
    }
  }

  if (it != _prototypes.end()) {
    return;
  }

  _prototypes.push_back(g.clone());
}

///
///
///
std::vector<std::string>
Registry::getGeometryTypes() const
{
  std::vector<std::string> names;

  for (const_prototype_iterator it = _prototypes.begin();
       it != _prototypes.end(); ++it) {
    names.push_back((*it)->geometryType());
  }

  return names;
}

///
///
///
Geometry *
Registry::newGeometryByTypeName(const std::string &geometryTypeName) const
{
  for (const_prototype_iterator it = _prototypes.begin();
       it != _prototypes.end(); ++it) {
    if (geometryTypeName == (*it)->geometryType()) {
      return (*it)->clone();
    }
  }

  SFCGAL_WARNING(boost::format("Registry can't create a new Geometry for the "
                               "type '%s' (returning null pointer)") %
                 geometryTypeName);
  return NULL;
}

///
///
///
Geometry *
Registry::newGeometryByTypeId(int typeId) const
{
  for (const_prototype_iterator it = _prototypes.begin();
       it != _prototypes.end(); ++it) {
    if (typeId == (*it)->geometryTypeId()) {
      return (*it)->clone();
    }
  }

  SFCGAL_WARNING(boost::format("Registry can't create a new Geometry for the "
                               "type '%d' (returning null pointer)") %
                 typeId);
  return NULL;
}

///
///
///
Registry &
Registry::instance()
{
  if (!Registry::_instance) {
    Registry::_instance = new Registry();
  }

  return *_instance;
}

///
///
///
Registry::Registry() : _prototypes()
{
  addPrototype(Point());
  addPrototype(LineString());
  addPrototype(Polygon());
  addPrototype(Triangle());
  addPrototype(Solid());

  addPrototype(GeometryCollection());

  addPrototype(MultiPoint());
  addPrototype(MultiLineString());
  addPrototype(MultiPolygon());
  addPrototype(MultiSolid());

  addPrototype(TriangulatedSurface());
  addPrototype(PolyhedralSurface());
}

} // namespace tools
} // namespace SFCGAL
