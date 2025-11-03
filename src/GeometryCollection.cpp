// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/format.hpp>
#include <memory>

#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/GeometryVisitor.h"

#include "SFCGAL/Exception.h"

namespace SFCGAL {

GeometryCollection::GeometryCollection() = default;

GeometryCollection::GeometryCollection(GeometryCollection const &other)
    : GeometryImpl(other)
{
  _geometries.reserve(other._geometries.size());
  for (const auto &geometry : other._geometries) {
    _geometries.emplace_back(geometry->clone());
  }
}

auto
GeometryCollection::operator=(GeometryCollection other) -> GeometryCollection &
{
  swap(other);
  return *this;
}

GeometryCollection::~GeometryCollection() = default;

auto
GeometryCollection::geometryType() const -> std::string
{
  return "GeometryCollection";
}

auto
GeometryCollection::geometryTypeId() const -> GeometryType
{
  return TYPE_GEOMETRYCOLLECTION;
}

auto
GeometryCollection::dimension() const -> int
{
  int maxDimension = 0;

  for (const auto &geometry : _geometries) {
    maxDimension = std::max(maxDimension, geometry->dimension());
  }

  return maxDimension;
}

auto
GeometryCollection::coordinateDimension() const -> int
{
  if (isEmpty()) {
    return 0;
  }
  return _geometries.front()->coordinateDimension();
}

auto
GeometryCollection::isEmpty() const -> bool
{
  return _geometries.empty();
}

auto
GeometryCollection::is3D() const -> bool
{
  return !isEmpty() && _geometries.front()->is3D();
}

auto
GeometryCollection::isMeasured() const -> bool
{
  return !isEmpty() && _geometries.front()->isMeasured();
}

auto
GeometryCollection::dropZ() -> bool
{
  if (!is3D()) {
    return false;
  }

  for (auto &geometry : _geometries) {
    geometry->dropZ();
  }

  return true;
}

auto
GeometryCollection::dropM() -> bool
{
  if (!isMeasured()) {
    return false;
  }

  for (auto &geometry : _geometries) {
    geometry->dropM();
  }

  return true;
}

auto
GeometryCollection::swapXY() -> void
{
  for (auto &geometry : _geometries) {
    geometry->swapXY();
  }
}

auto
GeometryCollection::numGeometries() const -> size_t
{
  return _geometries.size();
}

auto
GeometryCollection::geometryN(size_t const &n) const -> const Geometry &
{
  if (n >= numGeometries()) {
    BOOST_THROW_EXCEPTION(
        Exception((boost::format("Cannot access geometry at position %s. "
                                 "GeometryCollection has only %d geometries.") %
                   n % numGeometries())
                      .str()));
  }

  return *_geometries[n];
}

auto
GeometryCollection::geometryN(size_t const &n) -> Geometry &
{
  if (n >= numGeometries()) {
    BOOST_THROW_EXCEPTION(
        Exception((boost::format("Cannot access geometry at position %s. "
                                 "GeometryCollection has only %d geometries.") %
                   n % numGeometries())
                      .str()));
  }

  return *_geometries[n];
}

void
GeometryCollection::setGeometryN(std::unique_ptr<Geometry> geometry,
                                 size_t const             &idx)
{
  BOOST_ASSERT(geometry != nullptr);

  if (idx >= numGeometries()) {
    BOOST_THROW_EXCEPTION(
        Exception((boost::format("Cannot set geometry at position %s. "
                                 "GeometryCollection has only %d geometries.") %
                   idx % numGeometries())
                      .str()));
  }

  if (!isAllowed(*geometry)) {
    std::ostringstream oss;
    oss << "try to add a '" << geometry->geometryType() << "' in a '"
        << geometryType() << "'\n";
    BOOST_THROW_EXCEPTION(InappropriateGeometryException(oss.str()));
  }

  _geometries[idx] = std::move(geometry);
}

void
GeometryCollection::setGeometryN(Geometry *geometry, size_t const &idx)
{
  setGeometryN(std::unique_ptr<Geometry>(geometry), idx);
}

void
GeometryCollection::setGeometryN(const Geometry &geometry, size_t const &idx)
{
  setGeometryN(geometry.clone(), idx);
}

void
GeometryCollection::addGeometry(std::unique_ptr<Geometry> geometry)
{
  BOOST_ASSERT(geometry != nullptr);

  if (!isAllowed(*geometry)) {
    std::ostringstream oss;
    oss << "try to add a '" << geometry->geometryType() << "' in a '"
        << geometryType() << "'\n";
    BOOST_THROW_EXCEPTION(InappropriateGeometryException(oss.str()));
  }

  _geometries.emplace_back(std::move(geometry));
}

void
GeometryCollection::addGeometry(Geometry *geometry)
{
  addGeometry(std::unique_ptr<Geometry>(geometry));
}

void
GeometryCollection::addGeometry(Geometry const &geometry)
{
  addGeometry(geometry.clone());
}

auto
GeometryCollection::isAllowed(Geometry const & /*unused*/) -> bool
{
  // GeometryCollection accepts all subtypes
  return true;
}

void
GeometryCollection::accept(GeometryVisitor &visitor)
{
  return visitor.visit(*this);
}

void
GeometryCollection::accept(ConstGeometryVisitor &visitor) const
{
  return visitor.visit(*this);
}

} // namespace SFCGAL
