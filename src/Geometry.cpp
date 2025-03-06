// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/Geometry.h"

#include "SFCGAL/GeometryVisitor.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/detail/GetPointsVisitor.h"
#include "SFCGAL/detail/io/WkbWriter.h"
#include "SFCGAL/detail/io/WktWriter.h"

#include "SFCGAL/algorithm/BoundaryVisitor.h"
#include "SFCGAL/algorithm/centroid.h"
#include "SFCGAL/algorithm/distance.h"
#include "SFCGAL/algorithm/distance3d.h"

#include "SFCGAL/detail/EnvelopeVisitor.h"

#include "SFCGAL/detail/transform/RoundTransform.h"

#include "SFCGAL/Kernel.h"

namespace SFCGAL {

///
///
///
auto
Geometry::asText(const int &numDecimals) const -> std::string
{
  std::ostringstream oss;

  if (numDecimals >= 0) {
    oss << std::fixed;
    oss.precision(numDecimals);
  }

  detail::io::WktWriter writer(oss);
  bool                  exact = false;

  if (numDecimals == -1) {
    exact = true;
  }

  writer.write(*this, exact);
  return oss.str();
}

auto
Geometry::asWkb(boost::endian::order wkbOrder, bool asHex) const -> std::string
{
  std::ostringstream    oss;
  detail::io::WkbWriter writer(oss, asHex);
  writer.write(*this, wkbOrder);
  return oss.str();
}
///
///
///
auto
Geometry::envelope() const -> Envelope
{
  Envelope                box;
  detail::EnvelopeVisitor envelopeVisitor(box);
  accept(envelopeVisitor);
  return box;
}

///
///
///
auto
Geometry::boundary() const -> std::unique_ptr<Geometry>
{
  algorithm::BoundaryVisitor visitor;
  accept(visitor);
  return std::unique_ptr<Geometry>(visitor.releaseBoundary());
}

///
///
///
auto
Geometry::distance(const Geometry &other) const -> double
{
  return algorithm::distance(*this, other);
}

///
///
///
auto
Geometry::distance3D(const Geometry &other) const -> double
{
  return algorithm::distance3D(*this, other);
}

///
///
///
void
Geometry::round(const long &scale)
{
  transform::RoundTransform roundTransform(scale);
  accept(roundTransform);
}

///
///
///
auto
Geometry::numGeometries() const -> size_t
{
  return isEmpty() ? 0 : 1;
}

///
///
///
auto
Geometry::geometryN(size_t const &n) const -> const Geometry &
{
  BOOST_ASSERT(n == 0);
  (void)n;
  return *this;
}

///
///
///
auto
Geometry::geometryN(size_t const &n) -> Geometry &
{
  BOOST_ASSERT(n == 0);
  (void)n;
  return *this;
}

///
///
///
Geometry::Geometry() = default;

auto
Geometry::hasValidityFlag() const -> bool
{
  return validityFlag_;
}

void
Geometry::forceValidityFlag(bool valid)
{
  validityFlag_ = valid;
}

auto
Geometry::almostEqual(const Geometry &other, const double tolerance) const
    -> bool
{
  if (geometryTypeId() != other.geometryTypeId()) {
    return false;
  }

  detail::GetPointsVisitor get_points_a;
  detail::GetPointsVisitor get_points_b;
  accept(get_points_a);
  other.accept(get_points_b);

  if (get_points_a.points.size() != get_points_b.points.size()) {
    return false;
  }

  for (auto &point : get_points_a.points) {
    bool found = false;

    for (auto &j : get_points_b.points) {
      const Point &pta = *point;
      const Point &ptb = *j;

      if (pta.almostEqual(ptb, tolerance)) {
        found = true;
        break;
      }
    }

    if (!found) {
      return false;
    }
  }

  return true;
}

///
///
///
auto
Geometry::centroid() const -> Point
{
  std::unique_ptr<Point> out = algorithm::centroid(*this);

  return *(out.get());
}

///
///
///
auto
Geometry::centroid3D() const -> Point
{
  std::unique_ptr<Point> out = algorithm::centroid3D(*this);

  return *(out.get());
}

///
/// Function used to compare geometries
/// FIXME
/// Since we do not have (yet) a real "equals" operator, we only compare points
/// coordinates
auto
operator==(const Geometry &ga, const Geometry &gb) -> bool
{
  if (ga.geometryTypeId() != gb.geometryTypeId()) {
    return false;
  }

  detail::GetPointsVisitor get_points_a;
  detail::GetPointsVisitor get_points_b;
  ga.accept(get_points_a);
  gb.accept(get_points_b);

  if (get_points_a.points.size() != get_points_b.points.size()) {
    return false;
  }

  for (auto &point : get_points_a.points) {
    bool found = false;

    for (auto &j : get_points_b.points) {
      const Point &pta = *point;
      const Point &ptb = *j;

      if (pta == ptb) {
        found = true;
        break;
      }
    }

    if (!found) {
      return false;
    }
  }

  return true;
}

} // namespace SFCGAL
