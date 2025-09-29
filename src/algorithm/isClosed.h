// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later
#ifndef SFCGAL_ALGORITHM_ISCLOSED_H_
#define SFCGAL_ALGORITHM_ISCLOSED_H_
#include "SFCGAL/Geometry.h"
#include "SFCGAL/config.h"
namespace SFCGAL::algorithm {
/**
 * @brief Result type for closure check
 */
class SFCGAL_API Closure {
public:
  /**
   * @brief Default constructor (closed)
   */
  Closure() : _closed(true) {}
  /**
   * @brief Constructor with status and reason
   * @param closed true if the geometry is closed. false otherwise.
   * @param reason A string explaining why the geometry is open. Required only
   * if the Geometry is open.
   */
  Closure(bool closed, std::string reason = "")
      : _closed(closed), _reason(std::move(reason))
  {
  }
  /**
   * @brief Check if geometry is closed
   */
  operator bool() const { return _closed; }
  /**
   * @brief Get the reason why geometry is not closed
   * @return A string describing why the Geometry is open.
   */
  [[nodiscard]] auto
  reason() const -> const std::string &
  {
    return _reason;
  }
  /**
   * @brief Create a closed result
   * @return A closed result
   */
  static auto
  closed() -> Closure
  {
    return {true};
  }
  /**
   * @brief Create an open result with reason
   * @param reason A string describing why the Geometry is open.
   * @return An open result
   */
  static auto
  open(const std::string &reason) -> Closure
  {
    return {false, reason};
  }

private:
  bool        _closed;
  std::string _reason;
};
/**
 * @brief Check if a geometry is closed
 * @param g The geometry to check
 * @return Closure object with status and optional reason
 *
 * @note Definition of "closed" varies by geometry type:
 * - Point: Always closed
 * - LineString: Closed if first and last points are identical
 * - Polygon: Always closed (rings are closed by definition)
 * - Triangle: Always closed
 * - PolyhedralSurface: Closed if it forms a closed volume (no boundary edges)
 * - TriangulatedSurface: Closed if it forms a closed volume
 * - Solid: Always closed (by definition, but we test if the shells are closed)
 * - MultiPoint: Always closed
 * - MultiLineString: Closed if all LineStrings are closed
 * - MultiPolygon: Always closed
 * - MultiSolid: Always closed (by definition, cf Solid)
 * - GeometryCollection: Closed if all contained geometries are closed
 */
SFCGAL_API auto
isClosed(const Geometry &g) -> Closure;
} // namespace SFCGAL::algorithm
#endif // SFCGAL_ALGORITHM_ISCLOSED_H_
