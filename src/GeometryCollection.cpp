// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/GeometryVisitor.h>

#include <SFCGAL/Exception.h>

namespace SFCGAL {

///
///
///
GeometryCollection::GeometryCollection()  {}

///
///
///
GeometryCollection::GeometryCollection(GeometryCollection const &other)
    : Geometry(other)
{
  for (size_t i = 0; i < other.numGeometries(); i++) {
    addGeometry(other.geometryN(i).clone());
  }
}

///
///
///
auto
GeometryCollection::operator=(GeometryCollection other) -> GeometryCollection &
{
  swap(other);
  return *this;
}

///
///
///
GeometryCollection::~GeometryCollection() = default;

///
///
///
auto
GeometryCollection::clone() const -> GeometryCollection *
{
  return new GeometryCollection(*this);
}

///
///
///
auto
GeometryCollection::geometryType() const -> std::string
{
  return "GeometryCollection";
}

///
///
///
auto
GeometryCollection::geometryTypeId() const -> GeometryType
{
  return TYPE_GEOMETRYCOLLECTION;
}

///
///
///
auto
GeometryCollection::dimension() const -> int
{
  int maxDimension = 0;

  for (const auto &_geometrie : _geometries) {
    maxDimension = std::max(maxDimension, _geometrie.dimension());
  }

  return maxDimension;
}

///
///
///
auto
GeometryCollection::coordinateDimension() const -> int
{
  if (isEmpty()) {
    return 0;
  }     return _geometries.front().coordinateDimension();
 
}

///
///
///
auto
GeometryCollection::isEmpty() const -> bool
{
  return _geometries.empty();
}

///
///
///
auto
GeometryCollection::is3D() const -> bool
{
  return !isEmpty() && _geometries.front().is3D();
}

///
///
///
auto
GeometryCollection::isMeasured() const -> bool
{
  return !isEmpty() && _geometries.front().isMeasured();
}

///
///
///
auto
GeometryCollection::numGeometries() const -> size_t
{
  return _geometries.size();
}

///
///
///
auto
GeometryCollection::geometryN(size_t const &n) const -> const Geometry &
{
  return _geometries[n];
}

///
///
///
auto
GeometryCollection::geometryN(size_t const &n) -> Geometry &
{
  return _geometries[n];
}

///
///
///
void
GeometryCollection::addGeometry(Geometry *geometry)
{
  BOOST_ASSERT(geometry != NULL);

  if (!isAllowed(*geometry)) {
    std::ostringstream oss;
    oss << "try a add a '" << geometry->geometryType() << "' in a '"
        << geometryType() << "'";
    delete geometry; // we are responsible for the resource here
    BOOST_THROW_EXCEPTION(std::runtime_error(oss.str()));
  }

  _geometries.push_back(geometry);
}

///
///
///
void
GeometryCollection::addGeometry(Geometry const &geometry)
{
  addGeometry(geometry.clone());
}

///
///
///
auto
GeometryCollection::isAllowed(Geometry const & /*unused*/) -> bool
{
  // GeometryCollection accepts all subtypes
  return true;
}

///
///
///
void
GeometryCollection::accept(GeometryVisitor &visitor)
{
  return visitor.visit(*this);
}

///
///
///
void
GeometryCollection::accept(ConstGeometryVisitor &visitor) const
{
  return visitor.visit(*this);
}

} // namespace SFCGAL
