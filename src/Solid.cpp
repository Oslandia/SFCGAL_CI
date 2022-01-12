// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/GeometryVisitor.h>
#include <SFCGAL/Solid.h>

namespace SFCGAL {

///
///
///
Solid::Solid() { _shells.push_back(new PolyhedralSurface()); }

///
///
///
Solid::Solid(const PolyhedralSurface &exteriorShell)
{
  _shells.push_back(exteriorShell.clone());
}

///
///
///
Solid::Solid(PolyhedralSurface *exteriorShell)
{
  _shells.push_back(exteriorShell);
}

///
///
///
Solid::Solid(const std::vector<PolyhedralSurface> &shells)
{
  if (shells.empty()) {
    _shells.resize(1, new PolyhedralSurface());
  } else {
    for (size_t i = 0; i < shells.size(); i++) {
      _shells.push_back(shells[i].clone());
    }
  }
}

///
///
///
Solid::Solid(const Solid &other) : Geometry(other)
{
  for (size_t i = 0; i < other.numShells(); i++) {
    _shells.push_back(other.shellN(i).clone());
  }
}

///
///
///
Solid &
Solid::operator=(Solid other)
{
  swap(other);
  return *this;
}

///
///
///
Solid::~Solid() {}

///
///
///
Solid *
Solid::clone() const
{
  return new Solid(*this);
}

///
///
///
std::string
Solid::geometryType() const
{
  return "Solid";
}

///
///
///
GeometryType
Solid::geometryTypeId() const
{
  return TYPE_SOLID;
}

///
///
///
int
Solid::dimension() const
{
  return 3;
}

///
///
///
int
Solid::coordinateDimension() const
{
  return exteriorShell().coordinateDimension();
}

///
///
///
bool
Solid::isEmpty() const
{
  return exteriorShell().isEmpty();
}

///
///
///
bool
Solid::is3D() const
{
  return exteriorShell().is3D();
}

///
///
///
bool
Solid::isMeasured() const
{
  return exteriorShell().isMeasured();
}

///
///
///
void
Solid::accept(GeometryVisitor &visitor)
{
  return visitor.visit(*this);
}

///
///
///
void
Solid::accept(ConstGeometryVisitor &visitor) const
{
  return visitor.visit(*this);
}

} // namespace SFCGAL
