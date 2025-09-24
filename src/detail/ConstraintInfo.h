// Copyright (c) 2025-2025, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_DETAIL_CONSTRAINTINFO_H_
#define SFCGAL_DETAIL_CONSTRAINTINFO_H_

#include <algorithm>
#include <cstddef>
#include <vector>

namespace SFCGAL {
namespace detail {

/**
 * @brief Structure to identify constraint sources for topology-preserving mode
 */
struct ConstraintInfo {
  /**
   * @brief Constraint type enumeration
   */
  enum class Type {
    LINESTRING,
    POLYGON_EXTERIOR,
    POLYGON_INTERIOR,
    POLYHEDRALSURFACE_EXTERIOR,
    POLYHEDRALSURFACE_INTERIOR,
    MULTIPOLYGON_EXTERIOR,
    MULTIPOLYGON_INTERIOR
  };

  Type   type;      ///< The constraint type
  size_t geomIndex; ///< Index of the source geometry
  size_t ringIndex; ///< Index of the ring within the geometry

  /**
   * @brief Constructor for constraint info
   * @param t The constraint type
   * @param gIdx Index of the source geometry
   * @param rIdx Index of the ring within the geometry
   */
  ConstraintInfo(Type t, size_t gIdx, size_t rIdx = 0)
      : type(t), geomIndex(gIdx), ringIndex(rIdx)
  {
  }
};

/**
 * @brief Structure to store constraint order information
 * Allows recovering the original insertion order after simplification
 */
template <typename Constraint_id>
struct ConstraintOrderInfo {
  Constraint_id cid;       ///< Constraint ID returned by insertion
  size_t        geomIndex; ///< Index of the geometry in the collection
  size_t polyIndex; ///< Index of the polygon in MultiPolygon (if applicable)
  size_t ringIndex; ///< Ring index (0=exterior, >0=hole)
  size_t order; ///< Local order (for multiple constraints of the same polygon)
  ConstraintInfo::Type
      type; ///< Geometry type (linestring, polygon exterior, polygon interior)

  /**
   * @brief Constructor for constraint order info
   * @param cid_ Constraint ID returned by insertion
   * @param geomIndex_ Index of the geometry in the collection
   * @param polyIndex_ Index of the polygon in MultiPolygon (if applicable)
   * @param ringIndex_ Ring index (0=exterior, >0=hole)
   * @param order_ Local order (for multiple constraints of the same polygon)
   * @param type_ Geometry type
   */
  ConstraintOrderInfo(
      Constraint_id cid_, size_t geomIndex_, size_t polyIndex_ = 0,
      size_t ringIndex_ = 0, size_t order_ = 0,
      ConstraintInfo::Type type_ = ConstraintInfo::Type::LINESTRING)
      : cid(cid_), geomIndex(geomIndex_), polyIndex(polyIndex_),
        ringIndex(ringIndex_), order(order_), type(type_)
  {
  }
};

/**
 * @brief Comparison function for sorting constraints by original order
 */
template <typename Constraint_id>
struct ConstraintInfoCompare {
  /**
   * @brief Compare two constraint order infos
   * @param a First constraint order info
   * @param b Second constraint order info
   * @return True if a should be ordered before b
   */
  bool
  operator()(const ConstraintOrderInfo<Constraint_id> &a,
             const ConstraintOrderInfo<Constraint_id> &b) const
  {
    if (a.geomIndex != b.geomIndex)
      return a.geomIndex < b.geomIndex;
    if (a.polyIndex != b.polyIndex)
      return a.polyIndex < b.polyIndex;
    if (a.ringIndex != b.ringIndex)
      return a.ringIndex < b.ringIndex;
    return a.order < b.order;
  }
};

} // namespace detail
} // namespace SFCGAL

#endif // SFCGAL_DETAIL_CONSTRAINTINFO_H_
