// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/Curve.h"
#include "SFCGAL/Exception.h" // AJOUT : nécessaire pour BOOST_THROW_EXCEPTION

namespace SFCGAL {

///
/// Curve
///
Curve::Curve() = default;

///
/// Curve
///
Curve::Curve(const Curve &other) = default;

///
/// operator=
///
auto
Curve::operator=(const Curve &other) -> Curve & = default;

///
/// ~Curve
///
Curve::~Curve() = default;

///
/// startPoint
///
auto
Curve::startPoint() const -> Point
{
  auto bounds = parameterBounds();
  return evaluate(bounds.first);
}

///
/// endPoint
///
auto
Curve::endPoint() const -> Point
{
  auto bounds = parameterBounds();
  return evaluate(bounds.second);
}

///
/// midPoint
///
auto
Curve::midPoint() const -> Point
{
  auto bounds   = parameterBounds();
  auto midParam = (bounds.first + bounds.second) / FT(2);
  return evaluate(midParam);
}

///
/// isValidParameter
///
auto
Curve::isValidParameter(Parameter parameter) const -> bool
{
  auto bounds = parameterBounds();
  return parameter >= bounds.first && parameter <= bounds.second;
}

///
/// normalizeParameter
///
auto
Curve::normalizeParameter(Parameter parameter) const -> FT
{
  auto bounds = parameterBounds();
  FT   range  = bounds.second - bounds.first;

  if (CGAL::abs(range) < FT(1e-10)) {
    return FT(0);
  }

  return (parameter - bounds.first) / range;
}

///
/// denormalizeParameter
///
auto
Curve::denormalizeParameter(FT normalizedParameter) const -> Parameter
{
  auto bounds = parameterBounds();
  return bounds.first + normalizedParameter * (bounds.second - bounds.first);
}

} // namespace SFCGAL
