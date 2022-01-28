// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/detail/Interval.h>
#include <SFCGAL/numeric.h>

#include <algorithm>

namespace SFCGAL {
namespace detail {

///
///
///
Interval::Interval() : _lower(NaN()), _upper(NaN()) {}

///
///
///
Interval::Interval(const double &value) : _lower(value), _upper(value) {}

///
///
///
Interval::Interval(const double &v1, const double &v2)
    : _lower(std::min(v1, v2)), _upper(std::max(v1, v2))
{
}

///
///
///
Interval::Interval(const Interval &other)

    = default;

///
///
///
auto
Interval::operator=(const Interval &other) -> Interval & = default;

///
///
///
Interval::~Interval() = default;

///
///
///
auto
Interval::isEmpty() const -> bool
{
  return std::isnan(_lower) || std::isnan(_upper);
}

///
///
///
void
Interval::expandBy(const double &d)
{
  if (isEmpty()) {
    return;
  }

  _lower = _lower - d;
  _upper = _upper + d;
}

///
///
///
void
Interval::expandToInclude(const Interval &other)
{
  // ignore empty interval
  if (other.isEmpty()) {
    return;
  }

  if (isEmpty()) {
    (*this) = other;
  } else {
    _lower = std::min(_lower, other._lower);
    _upper = std::max(_upper, other._upper);
  }
}

///
///
///
void
Interval::expandToInclude(const double &value)
{
  if (std::isnan(value)) {
    return;
  }

  if (isEmpty()) {
    _lower = value;
    _upper = value;
  } else {
    _lower = std::min(_lower, value);
    _upper = std::max(_upper, value);
  }
}

///
///
///
auto
Interval::intersects(const Interval &other) const -> bool
{
  // empty intervals never intersects
  if (isEmpty() || other.isEmpty()) {
    return false;
  }

  return !(_lower > other._upper || _upper < other._lower);
}

///
///
///
auto
Interval::operator==(const Interval &other) const -> bool
{
  if (isEmpty() && other.isEmpty()) {
    return true;
  }

  return _lower == other._lower && _upper == other._upper;
}

///
///
///
auto
Interval::operator!=(const Interval &other) const -> bool
{
  return !((*this) == other);
}

} // namespace detail
} // namespace SFCGAL
