// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_PREPARED_GEOMETRY_H_
#define SFCGAL_PREPARED_GEOMETRY_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/config.h"

#include "SFCGAL/Envelope.h"

#include <boost/endian/conversion.hpp>
#include <boost/noncopyable.hpp>
#include <boost/optional.hpp>
#include <boost/serialization/split_member.hpp>

#include <stdint.h> // uint32_t

namespace SFCGAL {

class Geometry;

/**
 * A PreparedGeometry is a shell around a SFCGAL::Geometry.
 * It is used to store annex data, like SRID or cached computations
 *
 * It is noncopyable since it stores a std::unique_ptr<SFCGAL::Geometry>
 *
 */
class SFCGAL_API PreparedGeometry : public boost::noncopyable {
public:
  /**
   * Default constructor
   */
  PreparedGeometry();

  /**
   * Constructor
   * @param geometry pointer to the underlying SFCGAL::Geometry. Takes ownership
   * @param srid The spatial reference system identifier
   */
  PreparedGeometry(std::unique_ptr<Geometry> &&geometry, srid_t srid = 0);

  /**
   * Constructor
   * @param geometry pointer to the underlying SFCGAL::Geometry. Takes ownership
   * @param srid The spatial reference system identifier
   */
  PreparedGeometry(Geometry *geometry, srid_t srid = 0);

  ~PreparedGeometry();

  /**
   * Geometry accessors
   * @return Const reference to the geometry
   */
  auto
  geometry() const -> const Geometry &;
  /**
   * Geometry accessor
   * @return Reference to the geometry
   */
  auto
  geometry() -> Geometry &;

  /**
   * Geometry setter
   * @param geom The new geometry to set
   */
  void
  resetGeometry(Geometry *geom);

  /**
   * SRID read only accessor
   * @return The spatial reference system identifier
   */
  [[nodiscard]] auto
  SRID() const -> const srid_t &
  {
    return _srid;
  }

  /**
   * SRID accessor
   * @return Reference to the spatial reference system identifier
   */
  auto
  SRID() -> srid_t &
  {
    return _srid;
  }

  /**
   * Envelope accessor (using cache)
   * @return The envelope of the geometry
   */
  [[nodiscard]] auto
  envelope() const -> const Envelope &;

  /**
   * Resets the cache
   */
  void
  invalidateCache();

  /**
   * @brief Convert to an extended WKT (with SRID)
   * @param numDecimals number of decimals, -1 for keeping the exact rational
   * representation, if possible
   * @return Extended WKT representation with SRID
   */
  [[nodiscard]] auto
  asEWKT(const int &numDecimals = -1) const -> std::string;

  /**
   * @brief Convert to Extended Well-Known Binary
   * @param wkbOrder Byte order for WKB output
   * @param asHex Return as hexadecimal string if true
   * @return Extended WKB representation
   */
  auto
  asEWKB(boost::endian::order wkbOrder = boost::endian::order::native,
         bool                 asHex    = false) const -> std::string;

  /**
   * @brief Save prepared geometry to archive
   * @param ar Archive to save to
   */
  template <class Archive>
  void
  save(Archive &ar, const unsigned int /*version*/) const
  {
    ar & _srid;
    const Geometry *pgeom = _geometry.get();
    ar & pgeom;
  }

  /**
   * @brief Load prepared geometry from archive
   * @param ar Archive to load from
   */
  template <class Archive>
  void
  load(Archive &ar, const unsigned int /*version*/)
  {
    ar & _srid;
    Geometry *pgeom;
    ar & pgeom;
    _geometry.reset(pgeom);
  }

  /**
   * @brief Serialize prepared geometry to/from archive
   * @param ar Archive for serialization
   * @param version Serialization version
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version)
  {
    boost::serialization::split_member(ar, *this, version);
  }

protected:
  /// @brief Pointer to underlying Geometry
  std::unique_ptr<Geometry> _geometry;

  /// @brief SRID of the geometry
  srid_t _srid;

  /// @brief Cached bounding box of the geometry
  mutable boost::optional<Envelope> _envelope;
};

} // namespace SFCGAL

#endif
