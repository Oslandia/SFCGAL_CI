// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/tools/Registry.h"

#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/detail/tools/Log.h"

namespace SFCGAL::tools {

Registry *Registry::_instance = nullptr;

Registry::~Registry()
{
  for (auto &_prototype : _prototypes) {
    delete _prototype;
  }
}

void
Registry::addPrototype(const Geometry &g)
{
  // find prototype by name
  auto it = _prototypes.begin();

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

auto
Registry::getGeometryTypes() const -> std::vector<std::string>
{
  std::vector<std::string> names;

  names.reserve(_prototypes.size());
  for (auto *_prototype : _prototypes) {
    names.push_back(_prototype->geometryType());
  }

  return names;
}

auto
Registry::newGeometryByTypeName(const std::string &geometryTypeName) const
    -> Geometry *
{
  for (auto *_prototype : _prototypes) {
    if (geometryTypeName == _prototype->geometryType()) {
      return _prototype->clone();
    }
  }

  SFCGAL_WARNING(boost::format("Registry can't create a new Geometry for the "
                               "type '%s' (returning null pointer)") %
                 geometryTypeName);
  return nullptr;
}

auto
Registry::newGeometryByTypeId(int typeId) const -> Geometry *
{
  for (auto *_prototype : _prototypes) {
    if (typeId == _prototype->geometryTypeId()) {
      return _prototype->clone();
    }
  }

  SFCGAL_WARNING(boost::format("Registry can't create a new Geometry for the "
                               "type '%d' (returning null pointer)") %
                 typeId);
  return nullptr;
}

auto
Registry::instance() -> Registry &
{
  if (Registry::_instance == nullptr) {
    Registry::_instance = new Registry();
  }

  return *_instance;
}

Registry::Registry()
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

} // namespace SFCGAL::tools
