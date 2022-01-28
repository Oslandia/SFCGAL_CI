// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/Coordinate.h>

#include <SFCGAL/Exception.h>
#include <SFCGAL/Kernel.h>
#include <SFCGAL/numeric.h>

namespace SFCGAL {

///
///
///
Coordinate::Coordinate() : _storage(Coordinate::Empty()) {}

///
///
///
Coordinate::Coordinate(const Kernel::FT &x, const Kernel::FT &y)
    : _storage(Kernel::Point_2(x, y))
{
}

///
///
///
Coordinate::Coordinate(const Kernel::FT &x, const Kernel::FT &y,
                       const Kernel::FT &z)
    : _storage(Kernel::Point_3(x, y, z))
{
}

///
///
///
Coordinate::Coordinate(const double &x, const double &y)
{
  if (!std::isfinite(x) || !std::isfinite(y)) {
    BOOST_THROW_EXCEPTION(NonFiniteValueException(
        "cannot create coordinate with non finite value"));
  }

  _storage = Kernel::Point_2(x, y);
}

///
///
///
Coordinate::Coordinate(const double &x, const double &y, const double &z)
{
  if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
    BOOST_THROW_EXCEPTION(NonFiniteValueException(
        "cannot create coordinate with non finite value"));
  }

  _storage = Kernel::Point_3(x, y, z);
}

///
///
///
Coordinate::Coordinate(const Kernel::Point_2 &other) : _storage(other) {}

///
///
///
Coordinate::Coordinate(const Kernel::Point_3 &other) : _storage(other) {}

///
///
///
Coordinate::Coordinate(const Coordinate &other) = default;

///
///
///
auto
Coordinate::operator=(const Coordinate &other) -> Coordinate & = default;

///
///
///
Coordinate::~Coordinate() = default;

class CoordinateDimensionVisitor : public boost::static_visitor<int> {
public:
  auto
  operator()(const Coordinate::Empty &) const -> int
  {
    return 0;
  }
  auto
  operator()(const Kernel::Point_2 &) const -> int
  {
    return 2;
  }
  auto
  operator()(const Kernel::Point_3 &) const -> int
  {
    return 3;
  }
};

///
///
///
auto
Coordinate::coordinateDimension() const -> int
{
  CoordinateDimensionVisitor visitor;
  return boost::apply_visitor(visitor, _storage);
}

///
///
///
auto
Coordinate::isEmpty() const -> bool
{
  return _storage.which() == 0;
}

///
///
///
auto
Coordinate::is3D() const -> bool
{
  return _storage.which() == 2;
}

class GetXVisitor : public boost::static_visitor<Kernel::FT> {
public:
  auto
  operator()(const Coordinate::Empty &) const -> Kernel::FT
  {
    BOOST_THROW_EXCEPTION(
        Exception("trying to get an empty coordinate x value"));
    return 0;
  }
  auto
  operator()(const Kernel::Point_2 &storage) const -> Kernel::FT
  {
    return storage.x();
  }
  auto
  operator()(const Kernel::Point_3 &storage) const -> Kernel::FT
  {
    return storage.x();
  }
};

///
///
///
auto
Coordinate::x() const -> Kernel::FT
{
  GetXVisitor visitor;
  return boost::apply_visitor(visitor, _storage);
}

class GetYVisitor : public boost::static_visitor<Kernel::FT> {
public:
  auto
  operator()(const Coordinate::Empty &) const -> Kernel::FT
  {
    BOOST_THROW_EXCEPTION(
        Exception("trying to get an empty coordinate y value"));
    return 0;
  }
  auto
  operator()(const Kernel::Point_2 &storage) const -> Kernel::FT
  {
    return storage.y();
  }
  auto
  operator()(const Kernel::Point_3 &storage) const -> Kernel::FT
  {
    return storage.y();
  }
};

///
///
///
auto
Coordinate::y() const -> Kernel::FT
{
  GetYVisitor visitor;
  return boost::apply_visitor(visitor, _storage);
}

class GetZVisitor : public boost::static_visitor<Kernel::FT> {
public:
  auto
  operator()(const Coordinate::Empty &) const -> Kernel::FT
  {
    BOOST_THROW_EXCEPTION(
        Exception("trying to get an empty coordinate z value"));
    return 0;
  }
  auto
  operator()(const Kernel::Point_2 &) const -> Kernel::FT
  {
    return 0;
  }
  auto
  operator()(const Kernel::Point_3 &storage) const -> Kernel::FT
  {
    return storage.z();
  }
};

///
///
///
auto
Coordinate::z() const -> Kernel::FT
{
  GetZVisitor visitor;
  return boost::apply_visitor(visitor, _storage);
}

//----------------------

class RoundVisitor : public boost::static_visitor<> {
public:
  RoundVisitor(const long &scaleFactor) : _scaleFactor(scaleFactor) {}

  void
  operator()(Coordinate::Empty &) const
  {
  }
  void
  operator()(Kernel::Point_2 &storage) const
  {
    storage = Kernel::Point_2(_roundFT(storage.x()), _roundFT(storage.y()));
  }
  void
  operator()(Kernel::Point_3 &storage) const
  {
    storage = Kernel::Point_3(_roundFT(storage.x()), _roundFT(storage.y()),
                              _roundFT(storage.z()));
  }

private:
  long _scaleFactor;

  auto
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

class ToPoint2Visitor : public boost::static_visitor<Kernel::Point_2> {
public:
  auto
  operator()(const Coordinate::Empty &) const -> Kernel::Point_2
  {
    return Kernel::Point_2(CGAL::ORIGIN);
  }
  auto
  operator()(const Kernel::Point_2 &storage) const -> Kernel::Point_2
  {
    return storage;
  }
  auto
  operator()(const Kernel::Point_3 &storage) const -> Kernel::Point_2
  {
    return Kernel::Point_2(storage.x(), storage.y());
  }
};

///
///
///
auto
Coordinate::toPoint_2() const -> Kernel::Point_2
{
  ToPoint2Visitor visitor;
  return boost::apply_visitor(visitor, _storage);
}

class ToPoint3Visitor : public boost::static_visitor<Kernel::Point_3> {
public:
  auto
  operator()(const Coordinate::Empty & /*storage*/) const -> Kernel::Point_3
  {
    return Kernel::Point_3(CGAL::ORIGIN);
  }
  auto
  operator()(const Kernel::Point_2 &storage) const -> Kernel::Point_3
  {
    return Kernel::Point_3(storage.x(), storage.y(), 0.0);
  }
  auto
  operator()(const Kernel::Point_3 &storage) const -> Kernel::Point_3
  {
    return storage;
  }
};

///
///
///
auto
Coordinate::toPoint_3() const -> Kernel::Point_3
{
  ToPoint3Visitor visitor;
  return boost::apply_visitor(visitor, _storage);
}

///
///
///
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
  } else if (other.x() < x()) {
    return false;
  }

  // comparison along y
  if (y() < other.y()) {
    return true;
  } else if (other.y() < y()) {
    return false;
  }

  // comparison along z if possible
  if (is3D()) {
    if (z() < other.z()) {
      return true;
    } else if (other.z() < z()) {
      return false;
    }
  }

  // points are equals
  return false;
}

///
///
///
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

///
///
///
auto
Coordinate::operator!=(const Coordinate &other) const -> bool
{
  return !(*this == other);
}

auto
Coordinate::almostEqual(const Coordinate &other, const double tolerance) const
    -> bool
{
  if (isEmpty()) {
    return other.isEmpty();
  }

  const bool xEquality =
      SFCGAL::almostEqual(x(), other.x(), Kernel::FT(tolerance));
  const bool yEquality =
      SFCGAL::almostEqual(y(), other.y(), Kernel::FT(tolerance));
  if (is3D() || other.is3D()) {
    const bool zEquality =
        SFCGAL::almostEqual(z(), other.z(), Kernel::FT(tolerance));
    return xEquality && yEquality && zEquality;
  }
  return xEquality && yEquality;
}

} // namespace SFCGAL
