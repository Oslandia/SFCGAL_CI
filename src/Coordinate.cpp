// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/Coordinate.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/numeric.h"

namespace SFCGAL {

Coordinate::Coordinate() : _storage(Coordinate::Empty()) {}

Coordinate::Coordinate(const Kernel::FT &x, const Kernel::FT &y)
    : _storage(Kernel::Point_2(x, y))
{
}

Coordinate::Coordinate(const Kernel::FT &x, const Kernel::FT &y,
                       const Kernel::FT &z)
    : _storage(Kernel::Point_3(x, y, z))
{
}

Coordinate::Coordinate(const double &x, const double &y)
{
  if (!std::isfinite(x) || !std::isfinite(y)) {
    BOOST_THROW_EXCEPTION(NonFiniteValueException(
        "cannot create coordinate with non finite value"));
  }

  _storage = Kernel::Point_2(x, y);
}

Coordinate::Coordinate(const double &x, const double &y, const double &z)
{
  if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
    BOOST_THROW_EXCEPTION(NonFiniteValueException(
        "cannot create coordinate with non finite value"));
  }

  _storage = Kernel::Point_3(x, y, z);
}

Coordinate::Coordinate(const Kernel::Point_2 &other) : _storage(other) {}

Coordinate::Coordinate(const Kernel::Point_3 &other) : _storage(other) {}

Coordinate::Coordinate(const Coordinate &other) = default;

auto
Coordinate::operator=(const Coordinate &other) -> Coordinate & = default;

Coordinate::~Coordinate() = default;

/**
 * @brief Visitor to get coordinate dimension
 */
class CoordinateDimensionVisitor : public boost::static_visitor<int> {
public:
  /// @brief Handle empty coordinates
  /// @return 0 for empty coordinates
  auto
  operator()(const Coordinate::Empty & /*unused*/) const -> int
  {
    return 0;
  }
  /// @brief Handle 2D points
  /// @return 2 for 2D coordinates
  auto
  operator()(const Kernel::Point_2 & /*unused*/) const -> int
  {
    return 2;
  }
  /// @brief Handle 3D points
  /// @return 3 for 3D coordinates
  auto
  operator()(const Kernel::Point_3 & /*unused*/) const -> int
  {
    return 3;
  }
};

auto
Coordinate::coordinateDimension() const -> int
{
  CoordinateDimensionVisitor visitor;
  return boost::apply_visitor(visitor, _storage);
}

auto
Coordinate::isEmpty() const -> bool
{
  return _storage.which() == 0;
}

auto
Coordinate::is3D() const -> bool
{
  return _storage.which() == 2;
}

/**
 * @brief Visitor to get X coordinate value
 */
class GetXVisitor : public boost::static_visitor<Kernel::FT> {
public:
  /// @brief Handle empty coordinates
  /// @throws Exception when trying to get X from empty coordinate
  /// @return Never returns (throws exception)
  auto
  operator()(const Coordinate::Empty & /*unused*/) const -> Kernel::FT
  {
    BOOST_THROW_EXCEPTION(
        Exception("trying to get an empty coordinate x value"));
    return 0;
  }
  /// @brief Get X coordinate from 2D point
  /// @param storage The 2D point
  /// @return X coordinate value
  auto
  operator()(const Kernel::Point_2 &storage) const -> Kernel::FT
  {
    return storage.x();
  }
  /// @brief Get X coordinate from 3D point
  /// @param storage The 3D point
  /// @return X coordinate value
  auto
  operator()(const Kernel::Point_3 &storage) const -> Kernel::FT
  {
    return storage.x();
  }
};

auto
Coordinate::x() const -> Kernel::FT
{
  GetXVisitor visitor;
  return boost::apply_visitor(visitor, _storage);
}

/**
 * @brief Visitor to get Y coordinate value
 */
class GetYVisitor : public boost::static_visitor<Kernel::FT> {
public:
  /// @brief Handle empty coordinates
  /// @throws Exception when trying to get Y from empty coordinate
  /// @return Never returns (throws exception)
  auto
  operator()(const Coordinate::Empty & /*unused*/) const -> Kernel::FT
  {
    BOOST_THROW_EXCEPTION(
        Exception("trying to get an empty coordinate y value"));
    return 0;
  }
  /// @brief Get Y coordinate from 2D point
  /// @param storage The 2D point
  /// @return Y coordinate value
  auto
  operator()(const Kernel::Point_2 &storage) const -> Kernel::FT
  {
    return storage.y();
  }
  /// @brief Get Y coordinate from 3D point
  /// @param storage The 3D point
  /// @return Y coordinate value
  auto
  operator()(const Kernel::Point_3 &storage) const -> Kernel::FT
  {
    return storage.y();
  }
};

auto
Coordinate::y() const -> Kernel::FT
{
  GetYVisitor visitor;
  return boost::apply_visitor(visitor, _storage);
}

/**
 * @brief Visitor to get Z coordinate value
 */
class GetZVisitor : public boost::static_visitor<Kernel::FT> {
public:
  /// @brief Handle empty coordinates
  /// @throws Exception when trying to get Z from empty coordinate
  /// @return Never returns (throws exception)
  auto
  operator()(const Coordinate::Empty & /*unused*/) const -> Kernel::FT
  {
    BOOST_THROW_EXCEPTION(
        Exception("trying to get an empty coordinate z value"));
    return 0;
  }
  /// @brief Handle 2D points (no Z coordinate)
  /// @return 0 for 2D points (no Z coordinate)
  auto
  operator()(const Kernel::Point_2 & /*unused*/) const -> Kernel::FT
  {
    return 0;
  }
  /// @brief Get Z coordinate from 3D point
  /// @param storage The 3D point
  /// @return Z coordinate value
  auto
  operator()(const Kernel::Point_3 &storage) const -> Kernel::FT
  {
    return storage.z();
  }
};

auto
Coordinate::z() const -> Kernel::FT
{
  GetZVisitor visitor;
  return boost::apply_visitor(visitor, _storage);
}

//----------------------

/**
 * @brief Visitor to round coordinate values
 */
class RoundVisitor : public boost::static_visitor<> {
public:
  /// @brief Constructor with scale factor
  /// @param scaleFactor The scaling factor for rounding
  RoundVisitor(const long &scaleFactor) : _scaleFactor(scaleFactor) {}

  /// @brief Handle empty coordinates (no operation)
  void
  operator()(Coordinate::Empty & /*unused*/) const
  {
  }
  /// @brief Round 2D point coordinates
  /// @param storage The 2D point to round
  void
  operator()(Kernel::Point_2 &storage) const
  {
    storage = Kernel::Point_2(_roundFT(storage.x()), _roundFT(storage.y()));
  }
  /// @brief Round 3D point coordinates
  /// @param storage The 3D point to round
  void
  operator()(Kernel::Point_3 &storage) const
  {
    storage = Kernel::Point_3(_roundFT(storage.x()), _roundFT(storage.y()),
                              _roundFT(storage.z()));
  }

private:
  long _scaleFactor;

  [[nodiscard]] auto
  _roundFT(const Kernel::FT &v) const -> Kernel::FT
  {

#ifdef CGAL_USE_GMPXX
    ::mpq_class q(SFCGAL::round(v.exact() * _scaleFactor), _scaleFactor);
    q.canonicalize();
    return Kernel::FT(q);
#else
    return Kernel::FT(
        CGAL::Gmpq(SFCGAL::round(v.exact() * _scaleFactor), _scaleFactor));
#endif
  }
};

auto
Coordinate::round(const long &scaleFactor) -> Coordinate &
{
  RoundVisitor roundVisitor(scaleFactor);
  boost::apply_visitor(roundVisitor, _storage);
  return *this;
}

//----------------------

/**
 * @brief Visitor to convert coordinate to CGAL 2D point
 */
class ToPoint2Visitor : public boost::static_visitor<Kernel::Point_2> {
public:
  /// @brief Handle empty coordinates
  /// @return Origin point for empty coordinates
  auto
  operator()(const Coordinate::Empty & /*unused*/) const -> Kernel::Point_2
  {
    return Kernel::Point_2(CGAL::ORIGIN);
  }
  /// @brief Convert 2D point to 2D point (identity)
  /// @param storage The 2D point
  /// @return The same 2D point
  auto
  operator()(const Kernel::Point_2 &storage) const -> Kernel::Point_2
  {
    return storage;
  }
  /// @brief Convert 3D point to 2D point (drop Z)
  /// @param storage The 3D point
  /// @return 2D point with X,Y coordinates
  auto
  operator()(const Kernel::Point_3 &storage) const -> Kernel::Point_2
  {
    return Kernel::Point_2(storage.x(), storage.y());
  }
};

auto
Coordinate::toPoint_2() const -> Kernel::Point_2
{
  ToPoint2Visitor visitor;
  return boost::apply_visitor(visitor, _storage);
}

/**
 * @brief Visitor to convert coordinate to CGAL 3D point
 */
class ToPoint3Visitor : public boost::static_visitor<Kernel::Point_3> {
public:
  /// @brief Handle empty coordinates
  /// @return Origin point for empty coordinates
  auto
  operator()(const Coordinate::Empty & /*storage*/) const -> Kernel::Point_3
  {
    return Kernel::Point_3(CGAL::ORIGIN);
  }
  /// @brief Convert 2D point to 3D point (Z=0)
  /// @param storage The 2D point
  /// @return 3D point with X,Y coordinates and Z=0
  auto
  operator()(const Kernel::Point_2 &storage) const -> Kernel::Point_3
  {
    return Kernel::Point_3(storage.x(), storage.y(), 0.0);
  }
  /// @brief Convert 3D point to 3D point (identity)
  /// @param storage The 3D point
  /// @return The same 3D point
  auto
  operator()(const Kernel::Point_3 &storage) const -> Kernel::Point_3
  {
    return storage;
  }
};

auto
Coordinate::toPoint_3() const -> Kernel::Point_3
{
  ToPoint3Visitor visitor;
  return boost::apply_visitor(visitor, _storage);
}

auto
Coordinate::operator<(const Coordinate &other) const -> bool
{
  // no empty comparison
  if (isEmpty() || other.isEmpty()) {
    BOOST_THROW_EXCEPTION(
        Exception("try to compare empty points using a < b "));
  }

  // no mixed dimension comparison
  if ((is3D() && !other.is3D()) || (!is3D() && other.is3D())) {
    BOOST_THROW_EXCEPTION(
        Exception("try to compare empty points with different coordinate "
                  "dimension using a < b"));
  }

  // comparison along x
  if (x() < other.x()) {
    return true;
  }
  if (other.x() < x()) {
    return false;
  }

  // comparison along y
  if (y() < other.y()) {
    return true;
  }
  if (other.y() < y()) {
    return false;
  }

  // comparison along z if possible
  if (is3D()) {
    if (z() < other.z()) {
      return true;
    }
    if (other.z() < z()) {
      return false;
    }
  }

  // points are equals
  return false;
}

auto
Coordinate::operator==(const Coordinate &other) const -> bool
{
  if (isEmpty()) {
    return other.isEmpty();
  }

  if (is3D() || other.is3D()) {
    return x() == other.x() && y() == other.y() && z() == other.z();
  }
  return x() == other.x() && y() == other.y();
}

auto
Coordinate::operator!=(const Coordinate &other) const -> bool
{
  return !(*this == other);
}

auto
Coordinate::almostEqual(const Coordinate &other, const double tolerance) const
    -> bool
{
  bool result{true};
  if (isEmpty()) {
    return result;
  }

  // no mixed dimension comparison
  if ((is3D() && !other.is3D()) || (!is3D() && other.is3D())) {
    BOOST_THROW_EXCEPTION(
        Exception("try to compare points with different coordinate "
                  "dimension using a.almostEqual(b)"));
  }

  result = SFCGAL::almostEqual(x(), other.x(), Kernel::FT(tolerance)) &&
           SFCGAL::almostEqual(y(), other.y(), Kernel::FT(tolerance));
  if (is3D()) {
    result &= SFCGAL::almostEqual(z(), other.z(), Kernel::FT(tolerance));
  }
  return result;
}

auto
Coordinate::dropZ() -> bool
{
  if (!is3D()) {
    return false;
  }

  _storage = Kernel::Point_2(x(), y());
  return true;
}

auto
Coordinate::swapXY() -> void
{
  if (_storage.which() == 1) {
    _storage = Kernel::Point_2(y(), x());
  } else if (_storage.which() == 2) {
    _storage = Kernel::Point_3(y(), x(), z());
  }
}

} // namespace SFCGAL
