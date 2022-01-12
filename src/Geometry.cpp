// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/Geometry.h>

#include <SFCGAL/GeometryVisitor.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/detail/GetPointsVisitor.h>
#include <SFCGAL/detail/io/WktWriter.h>

#include <SFCGAL/algorithm/BoundaryVisitor.h>
#include <SFCGAL/algorithm/distance.h>
#include <SFCGAL/algorithm/distance3d.h>

#include <SFCGAL/detail/EnvelopeVisitor.h>

#include <SFCGAL/detail/transform/RoundTransform.h>

#include <SFCGAL/Kernel.h>

namespace SFCGAL {

///
///
///
std::string
Geometry::asText(const int &numDecimals) const
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

///
///
///
Envelope
Geometry::envelope() const
{
  Envelope                box;
  detail::EnvelopeVisitor envelopeVisitor(box);
  accept(envelopeVisitor);
  return box;
}

///
///
///
std::unique_ptr<Geometry>
Geometry::boundary() const
{
  algorithm::BoundaryVisitor visitor;
  accept(visitor);
  return std::unique_ptr<Geometry>(visitor.releaseBoundary());
}

///
///
///
double
Geometry::distance(const Geometry &other) const
{
  return algorithm::distance(*this, other);
}

///
///
///
double
Geometry::distance3D(const Geometry &other) const
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
size_t
Geometry::numGeometries() const
{
  return 1;
}

///
///
///
const Geometry &
Geometry::geometryN(size_t const &n) const
{
  BOOST_ASSERT(n == 0);
  (void)n;
  return *this;
}

///
///
///
Geometry &
Geometry::geometryN(size_t const &n)
{
  BOOST_ASSERT(n == 0);
  (void)n;
  return *this;
}

///
///
///
Geometry::Geometry() : validityFlag_(false) {}

bool
Geometry::hasValidityFlag() const
{
  return validityFlag_;
}

void
Geometry::forceValidityFlag(bool valid)
{
  validityFlag_ = valid;
}

///
/// Function used to compare geometries
/// FIXME
/// Since we do not have (yet) a real "equals" operator, we only compare points
/// coordinates
bool
operator==(const Geometry &ga, const Geometry &gb)
{
  if (ga.geometryTypeId() != gb.geometryTypeId()) {
    return false;
  }

  detail::GetPointsVisitor get_points_a, get_points_b;
  ga.accept(get_points_a);
  gb.accept(get_points_b);

  if (get_points_a.points.size() != get_points_b.points.size()) {
    return false;
  }

  for (size_t i = 0; i < get_points_a.points.size(); ++i) {
    bool found = false;

    for (size_t j = 0; j < get_points_b.points.size(); ++j) {
      const Point &pta = *(get_points_a.points[i]);
      const Point &ptb = *(get_points_b.points[j]);

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
