// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/primitive3d/Cube.h"

namespace SFCGAL {

Cube::Cube(const Kernel::FT &size) : m_box(Box(size, size, size))
{
  // m_parameters["size"] is not necessary because the information is stored in
  // m_box However, it is needed to be compatible with Primitive loginc
  m_parameters["size"] = m_box.xExtent();
}

auto
Cube::primitiveType() const -> std::string
{
  return "Cube";
}

auto
Cube::primitiveTypeId() const -> PrimitiveType
{
  return PrimitiveType::TYPE_CUBE;
}

void
Cube::setSize(const Kernel::FT &size)
{
  validateAndSetParameter("size", size);
}

void
Cube::validateParameters(
    std::unordered_map<std::string, PrimitiveParameter> const &tempParameters)
    const
{
  const double size =
      CGAL::to_double(std::get<Kernel::FT>(tempParameters.at("size")));

  if (size < 0.) {
    BOOST_THROW_EXCEPTION(Exception("Cube size cannot be negative."));
  }
}

void
Cube::onValidatedAndSetParameter(const std::string        &name,
                                 const PrimitiveParameter &parameter)
{
  (void)parameter; // unused

  if (name == "size") {
    const Kernel::FT &value = size();

    m_box.setXExtent(value);
    m_box.setYExtent(value);
    m_box.setZExtent(value);
  }
}

auto
Cube::generatePolyhedralSurface() const -> PolyhedralSurface
{
  return m_box.generatePolyhedralSurface();
}

auto
Cube::area3D(bool withDiscretization) const -> double
{
  return m_box.area3D(withDiscretization);
}

auto
Cube::volume(bool withDiscretization) const -> double
{
  return m_box.volume(withDiscretization);
}

auto
Cube::size() const -> const Kernel::FT &
{
  return std::get<Kernel::FT>(m_parameters.at("size"));
}

auto
Cube::toString() const -> std::string
{
  std::ostringstream stringStream;
  stringStream << "[Primitive type: " << primitiveType() << ", size: " << size()
               << ", box: ";
  stringStream << m_box;
  stringStream << "]";

  return stringStream.str();
}

} // namespace SFCGAL
