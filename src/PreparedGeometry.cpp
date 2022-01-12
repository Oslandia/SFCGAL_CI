// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/PreparedGeometry.h>

#include <SFCGAL/detail/io/WktWriter.h>

namespace SFCGAL {
PreparedGeometry::PreparedGeometry() : _srid(0) {}

PreparedGeometry::PreparedGeometry(std::unique_ptr<Geometry> &&geometry,
                                   srid_t                      srid)
    : _geometry(geometry.release()), _srid(srid)
{
}

PreparedGeometry::PreparedGeometry(Geometry *geometry, srid_t srid)
    : _geometry(geometry), _srid(srid)
{
}

PreparedGeometry::~PreparedGeometry() = default;

auto
PreparedGeometry::geometry() const -> const Geometry &
{
  BOOST_ASSERT(_geometry.get());
  return *_geometry;
}

auto
PreparedGeometry::geometry() -> Geometry &
{
  BOOST_ASSERT(_geometry.get());
  return *_geometry;
}

void
PreparedGeometry::resetGeometry(Geometry *geom)
{
  _geometry.reset(geom);
  invalidateCache();
}

auto
PreparedGeometry::envelope() const -> const Envelope &
{
  if (!_envelope) {
    _envelope.reset(_geometry->envelope());
  }

  return *_envelope;
}

void
PreparedGeometry::invalidateCache()
{
  _envelope.reset();
}

auto
PreparedGeometry::asEWKT(const int &numDecimals) const -> std::string
{
  std::ostringstream oss;

  if (numDecimals >= 0) {
    oss << std::fixed;
    oss.precision(numDecimals);
  }

  if (_srid != 0) {
    oss << "SRID=" << _srid << ";";
  }

  detail::io::WktWriter writer(oss);
  bool                  exactWrite = false;

  if (numDecimals == -1) {
    exactWrite = true;
  }

  writer.write(*_geometry, exactWrite);
  return oss.str();
}
} // namespace SFCGAL
