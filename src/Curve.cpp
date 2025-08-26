// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/Curve.h"
#include "SFCGAL/Exception.h"

namespace SFCGAL {

Curve::Curve() = default;

Curve::Curve(const Curve &other) = default;

auto
Curve::operator=(const Curve &other) -> Curve & = default;

Curve::~Curve() = default;

auto
Curve::startPoint() const -> Point
{
  auto bounds = parameterBounds();
  return evaluate(bounds.first);
}

auto
Curve::endPoint() const -> Point
{
  auto bounds = parameterBounds();
  return evaluate(bounds.second);
}

auto
Curve::midPoint() const -> Point
{
  auto bounds   = parameterBounds();
  auto midParam = (bounds.first + bounds.second) / FT(2);
  return evaluate(midParam);
}

auto
Curve::isValidParameter(const Parameter &parameter) const -> bool
{
  const auto      bounds = parameterBounds();
  const Parameter startParam =
      (bounds.first < bounds.second) ? bounds.first : bounds.second;
  const Parameter endParam =
      (bounds.first < bounds.second) ? bounds.second : bounds.first;
  return parameter >= startParam && parameter <= endParam;
}

auto
Curve::normalizeParameter(const Parameter &parameter) const -> FT
{
  const auto      bounds = parameterBounds();
  const Parameter startParam =
      (bounds.first < bounds.second) ? bounds.first : bounds.second;
  const Parameter endParam =
      (bounds.first < bounds.second) ? bounds.second : bounds.first;
  FT range = endParam - startParam;
  if (range == FT(0)) {
    return {0};
  }
  // Optional: keep epsilon guard for inexact kernels
  if (CGAL::abs(range) < FT(1e-10)) {
    return {0};
  }
  return (parameter - startParam) / range;
}

auto
Curve::denormalizeParameter(const FT &normalizedParameter) const -> Parameter
{
  const auto      bounds = parameterBounds();
  const Parameter startParam =
      (bounds.first < bounds.second) ? bounds.first : bounds.second;
  const Parameter endParam =
      (bounds.first < bounds.second) ? bounds.second : bounds.first;
  return startParam + normalizedParameter * (endParam - startParam);
}

} // namespace SFCGAL
