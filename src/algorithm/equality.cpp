// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/equality.h"

#include "SFCGAL/algorithm/covers.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/detail/GetPointsVisitor.h"

#include <boost/format.hpp>

#include <vector>

namespace SFCGAL::algorithm {

auto
EqualityStrictness::operator|(Flag flag) -> EqualityStrictness &
{
  // can not have CheckCoverOrPoint and any of InternalPoint* checks
  if ((((_flags & CheckCoverOrPoint) != 0) &&
       flag >= InternalPointOrdered) || //
      (((flag & CheckCoverOrPoint) != 0) && _flags >= InternalPointOrdered)) {
    BOOST_THROW_EXCEPTION(
        Exception((boost::format("Conflict in EqualityStrictness flags: can "
                                 "not have CheckCoverOrPoint and any of "
                                 "InternalPoint* checks ('%s' vs '%s')") %
                   toString() % EqualityStrictness(flag).toString())
                      .str()));
  }

  // can not have multiple InternalPoint* checks
  if (_flags >= InternalPointOrdered && flag >= InternalPointOrdered) {
    BOOST_THROW_EXCEPTION(Exception(
        (boost::format("Conflict in EqualityStrictness flags: can not have "
                       "multiple InternalPoint* checks ('%s' vs '%s')") %
         toString() % EqualityStrictness(flag).toString())
            .str()));
  }

  _flags |= flag;
  return *this;
}

auto
EqualityStrictness::toString() const -> std::string
{
  std::string out;
  if ((_flags & CheckCoverOrPoint) != 0) {
    out = "CheckCover";
  } else {
    out = "CheckPoint";
  }
  if ((_flags & SubGeomOrdered) != 0) {
    out += " | SubGeomOrdered";
  }
  if ((_flags & SubPartOrdered) != 0) {
    out += " | SubPartOrdered";
  }
  if ((_flags & InternalPointOrdered) != 0) {
    out += " | InternalPointOrdered";
  }
  if ((_flags & InternalPointShifted) != 0) {
    out += " | InternalPointShifted";
  }
  if ((_flags & InternalPointInverted) != 0) {
    out += " | InternalPointInverted";
  }

  return out;
}

/// @{
/// @privatesection
#ifndef DOXYGEN_SHOULD_SKIP_THIS

/**
 * @brief check if the sub parts respect strictness, tolerance and the sub parts
 * are in the same order.
 * @param geomA geometry to compare with
 * @param geomB geometry to compare with
 * @param tolerance allowed distance between same points.
 * @param strictness will refer to EqualityStrictness flags
 * @pre geomA and geomB have the same type and this type has sub parts (not a
 * collection of geometries)
 * @return true if geomA and geomB are almost equal regard to strictness and
 * tolerance
 */
template <class G>
auto
compareAnySubPartOrdered(const G &geomA, const G &geomB, const double tolerance,
                         EqualityStrictness strictness) -> bool
{
  auto iteB = geomB.begin();
  for (auto iteA = geomA.begin(); iteA != geomA.end(); ++iteA) {
    bool found = almostEqual((*iteA), (*iteB), tolerance, strictness);
    if (!found) {
      return false;
    }
    ++iteB;
  }

  return true;
}

/**
 * @brief check if the sub parts respect strictness, tolerance and the sub parts
 * are in any order.
 * @param geomA geometry to compare with
 * @param geomB geometry to compare with
 * @param tolerance allowed distance between same points.
 * @param strictness will refer to EqualityStrictness flags
 * @pre geomA and geomB have the same type and this type has sub parts (not a
 * collection of geometries)
 * @return true if geomA and geomB are almost equal regard to strictness and
 * tolerance
 */
template <class G>
auto
compareAnySubPartNonOrdered(size_t numPart, const G &geomA, const G &geomB,
                            const double       tolerance,
                            EqualityStrictness strictness) -> bool
{
  std::vector<bool> hasGoodMatch(numPart);
  for (auto iteA = geomA.begin(); iteA != geomA.end(); ++iteA) {
    bool found = false;
    int  bIdx  = 0;
    for (auto iteB = geomB.begin(); !found && iteB != geomB.end(); ++iteB) {
      if (hasGoodMatch[bIdx]) { // already watched this geom
        ++bIdx;
        continue;
      }
      found              = almostEqual((*iteA), (*iteB), tolerance, strictness);
      hasGoodMatch[bIdx] = found;
      ++bIdx;
    }
    if (!found) {
      return false;
    }
  }

  return true;
}

/**
 * @brief check if the sub parts respect strictness, tolerance and the sub parts
 * are in the same order.
 * @param geomA geometry to compare with
 * @param geomB geometry to compare with
 * @param tolerance allowed distance between same points.
 * @param strictness will refer to EqualityStrictness flags
 * @return true if geomA and geomB are geometries with sub parts (not a
 * collection of geometries) and are almost equal regard to strictness and
 * tolerance
 * @see hasSubPart()
 */
auto
compareSubPartOrdered(const Geometry &geomA, const Geometry &geomB,
                      const double tolerance, EqualityStrictness strictness)
    -> bool
{
  switch (geomA.geometryTypeId()) {
  case TYPE_POINT:
  case TYPE_LINESTRING:
  case TYPE_TRIANGLE:
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
    // geom type is not valid for this kind of check
    return false;

  case TYPE_POLYGON:
    return compareAnySubPartOrdered<Polygon>(
        geomA.as<Polygon>(), geomB.as<Polygon>(), tolerance, strictness);

  case TYPE_TRIANGULATEDSURFACE:
    return compareAnySubPartOrdered<TriangulatedSurface>(
        geomA.as<TriangulatedSurface>(), geomB.as<TriangulatedSurface>(),
        tolerance, strictness);

  case TYPE_POLYHEDRALSURFACE:
    return compareAnySubPartOrdered<PolyhedralSurface>(
        geomA.as<PolyhedralSurface>(), geomB.as<PolyhedralSurface>(), tolerance,
        strictness);

  case TYPE_SOLID:
    return compareAnySubPartOrdered<Solid>(geomA.as<Solid>(), geomB.as<Solid>(),
                                           tolerance, strictness);

  default:
    BOOST_THROW_EXCEPTION(NotImplementedException(
        (boost::format("Unmanaged geometry type: %d.") % geomA.geometryType())
            .str()));
  }
}

/**
 * @brief check if the sub parts respect strictness, tolerance and the sub parts
 * are in any order.
 * @param geomA geometry to compare with
 * @param geomB geometry to compare with
 * @param tolerance allowed distance between same points.
 * @param strictness will refer to EqualityStrictness flags
 * @return true if geomA and geomB are geometries with sub parts (not a
 * collection of geometries) and are almost equal regard to strictness and
 * tolerance
 * @see hasSubPart()
 */
auto
compareSubPartNonOrdered(const Geometry &geomA, const Geometry &geomB,
                         const double tolerance, EqualityStrictness strictness)
    -> bool
{
  switch (geomA.geometryTypeId()) {
  case TYPE_POINT:
  case TYPE_LINESTRING:
  case TYPE_TRIANGLE:
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
    // geom type is not valid for this kind of check
    return false;

  case TYPE_POLYGON: {
    const auto &tempGA = geomA.as<Polygon>();
    return compareAnySubPartNonOrdered<Polygon>(
        tempGA.numRings(), tempGA, geomB.as<Polygon>(), tolerance, strictness);
  }

  case TYPE_TRIANGULATEDSURFACE: {
    const auto &tempGA = geomA.as<TriangulatedSurface>();
    return compareAnySubPartNonOrdered<TriangulatedSurface>(
        tempGA.numPatches(), tempGA, geomB.as<TriangulatedSurface>(), tolerance,
        strictness);
  }

  case TYPE_POLYHEDRALSURFACE: {
    const auto &tempGA = geomA.as<PolyhedralSurface>();
    return compareAnySubPartNonOrdered<PolyhedralSurface>(
        tempGA.numPatches(), tempGA, geomB.as<PolyhedralSurface>(), tolerance,
        strictness);
  }

  case TYPE_SOLID: {
    const auto &tempGA = geomA.as<Solid>();
    return compareAnySubPartNonOrdered<Solid>(
        tempGA.numShells(), tempGA, geomB.as<Solid>(), tolerance, strictness);
  }

  default:
    BOOST_THROW_EXCEPTION(NotImplementedException(
        (boost::format("Unmanaged geometry type: %d.") % geomA.geometryType())
            .str()));
  }
}

/**
 * @brief check if the sub geometries respect strictness, tolerance and the sub
 * geometries are in the same order.
 * @param geomA geometry to compare with
 * @param geomB geometry to compare with
 * @param tolerance allowed distance between same points.
 * @param strictness will refer to EqualityStrictness flags
 * @pre geomA and geomB have the same type and this type has sub geometries (not
 * sub parts)
 * @return true if geomA and geomB are almost equal regard to strictness and
 * tolerance
 */
auto
compareSubGeometryOrdered(const Geometry &geomA, const Geometry &geomB,
                          const double tolerance, EqualityStrictness strictness)
    -> bool
{
  for (size_t i = 0; i < geomA.numGeometries(); ++i) {
    bool found = almostEqual(geomA.geometryN(i), geomB.geometryN(i), tolerance,
                             strictness);
    if (!found) {
      return false;
    }
  }

  return true;
}

/**
 * @brief check if the sub geometries respect strictness, tolerance and the sub
 * geometries are in any order.
 * @param geomA geometry to compare with
 * @param geomB geometry to compare with
 * @param tolerance allowed distance between same points.
 * @param strictness will refer to EqualityStrictness flags
 * @pre geomA and geomB have the same type and this type has sub geometries (not
 * sub parts)
 * @return true if geomA and geomB are almost equal regard to strictness and
 * tolerance
 */
auto
compareSubGeometryNonOrdered(const Geometry &geomA, const Geometry &geomB,
                             const double       tolerance,
                             EqualityStrictness strictness) -> bool
{
  std::vector<bool> hasGoodMatch(geomA.numGeometries());
  for (size_t i = 0; i < geomA.numGeometries(); ++i) {
    bool found = false;
    for (size_t j = 0; !found && j < geomA.numGeometries(); ++j) {
      if (hasGoodMatch[j]) { // already watched this geom
        continue;
      }
      found = almostEqual(geomA.geometryN(i), geomB.geometryN(j), tolerance,
                          strictness);
      hasGoodMatch[j] = found;
    }
    if (!found) {
      return false;
    }
  }

  return true;
}

/**
 * @brief check if the geometry points respect strictness, tolerance and are in
 * the same order.
 * @param geomA geometry to compare with
 * @param geomB geometry to compare with
 * @param tolerance allowed distance between same points.
 * @param strictness will refer to EqualityStrictness flags
 * @pre geomA and geomB have the same type, the same number of point
 * @return true if geomA and geomB are almost equal regard to strictness and
 * tolerance
 */
auto
comparePointsOrdered(detail::GetPointsVisitor &getPointsA,
                     detail::GetPointsVisitor &getPointsB,
                     const double              tolerance) -> bool
{
  for (size_t i = 0; i < getPointsA.points.size(); ++i) {
    const Point &pta = *getPointsA.points[i];
    const Point &ptb = *getPointsB.points[i];
    bool         found =
        (tolerance < 0.0 && pta == ptb) || pta.almostEqual(ptb, tolerance);
    if (!found) {
      return false;
    }
  }

  return true;
}

/**
 * @brief check if the geometry points respect strictness, tolerance and are in
 * any order.
 * @param geomA geometry to compare with
 * @param geomB geometry to compare with
 * @param tolerance allowed distance between same points.
 * @param strictness will refer to EqualityStrictness flags
 * @pre geomA and geomB have the same type, the same number of point
 * @return true if geomA and geomB are almost equal regard to strictness and
 * tolerance
 */
auto
comparePointsNonOrdered(detail::GetPointsVisitor &getPointsA,
                        detail::GetPointsVisitor &getPointsB,
                        const double              tolerance) -> bool
{
  bool isClosed = *(getPointsA.points[0]) ==
                  *(getPointsA.points[getPointsA.points.size() - 1]);
  std::vector<bool> hasGoodMatch(getPointsA.points.size());
  for (size_t i = 0;                                      //
       i < getPointsA.points.size() - (isClosed ? 1 : 0); //
       ++i) {
    bool         found = false;
    const Point &pta   = *getPointsA.points[i];

    for (size_t j = 0; j < getPointsB.points.size(); ++j) {
      const Point &ptb = *getPointsB.points[j];
      if (hasGoodMatch[j]) { // already watched this point
        continue;
      }

      if ((tolerance < 0.0 && pta == ptb) || pta.almostEqual(ptb, tolerance)) {
        found           = true;
        hasGoodMatch[j] = true;
        break;
      }
    }

    if (!found) {
      return false;
    }
  }

  return true;
}

/**
 * @brief check if the geometry points respect strictness, tolerance and can be
 * shifted.
 *
 * Shifted means for opened linestring:
 * geom(1, 2, 3, 4) <==> geom(2, 3, 4, 1) <==> geom(3, 4, 1, 2) <==> geom(4, 1,
 * 2, 3)
 *
 * or for closed linestring:
 * geom(1, 2, 3, 4, 1) <==> geom(2, 3, 4, 1, 2) <==> geom(3, 4, 1, 2, 3) <==>
 * geom(4, 1, 2, 3, 4)
 *
 * @param geomA geometry to compare with
 * @param geomB geometry to compare with
 * @param tolerance allowed distance between same points.
 * @param strictness will refer to EqualityStrictness flags
 * @pre geomA and geomB have the same type, the same number of point
 * @return true if geomA and geomB are almost equal regard to strictness and
 * tolerance
 */
auto
comparePointsShifted(detail::GetPointsVisitor &getPointsA,
                     detail::GetPointsVisitor &getPointsB,
                     const double              tolerance) -> bool
{
  long startPos = -1;
  bool isClosed = *(getPointsA.points[0]) ==
                  *(getPointsA.points[getPointsA.points.size() - 1]);

  for (size_t i = 0;                                      //
       i < getPointsA.points.size() - (isClosed ? 1 : 0); //
       ++i) {
    // search first matching point
    if (startPos == -1) {
      const Point &pta = *getPointsA.points[i];
      for (size_t j = 0; j < getPointsB.points.size(); ++j) {

        const Point &ptb = *getPointsB.points[j];

        bool found = ((tolerance < 0.0 && pta == ptb) //
                      || pta.almostEqual(ptb, tolerance));
        if (found) {
          startPos = j; // found first point
          ++i;          // go for second point
          break;
        }
      }
    }
    // no first point found
    if (startPos == -1) {
      return false;
    }

    // handle next pos and array limits
    ++startPos;
    if (static_cast<size_t>(startPos) == getPointsB.points.size()) {
      startPos = (isClosed ? 1 : 0);
    }

    const Point &pta = *getPointsA.points[i];
    const Point &ptb = *getPointsB.points[startPos];

    bool found = ((tolerance < 0.0 && pta == ptb) //
                  || pta.almostEqual(ptb, tolerance));
    if (!found) {
      return false;
    }
  }

  return true;
}

/**
 * @brief check if the geometry points respect strictness, tolerance and are
 * inverted.
 *
 * Inverted means:
 * geom(1, 2, 3, 4) <==> geom(4, 3, 2, 1)
 *
 * @param geomA geometry to compare with
 * @param geomB geometry to compare with
 * @param tolerance allowed distance between same points.
 * @param strictness will refer to EqualityStrictness flags
 * @pre geomA and geomB have the same type, the same number of point
 * @return true if geomA and geomB are almost equal regard to strictness and
 * tolerance
 */
auto
comparePointsInverted(detail::GetPointsVisitor &getPointsA,
                      detail::GetPointsVisitor &getPointsB,
                      const double              tolerance) -> bool
{
  // first try with same order
  if (!comparePointsOrdered(getPointsA, getPointsB, tolerance)) {
    // second try with inverted order
    for (size_t i = 0; i < getPointsA.points.size(); ++i) {
      const Point &pta = *getPointsA.points[getPointsA.points.size() - 1 - i];
      const Point &ptb = *getPointsB.points[i];
      bool         found =
          (tolerance < 0.0 && pta == ptb) || pta.almostEqual(ptb, tolerance);
      if (!found) {
        return false;
      }
    }
  }
  return true;
}

/**
 * @param geom a geometry to check
 * @return true if geom has sub parts (not a collection of geometries)
 */
auto
hasSubPart(const Geometry &geom) -> bool
{
  switch (geom.geometryTypeId()) {
  case TYPE_POINT:
  case TYPE_LINESTRING:
  case TYPE_TRIANGLE:
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
    return false;

  case TYPE_POLYGON:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
  case TYPE_SOLID:
    return true;

  default:
    BOOST_THROW_EXCEPTION(NotImplementedException(
        (boost::format("Unmanaged geometry type: %d.") % geom.geometryType())
            .str()));
  }
}

#endif // ifndef DOXYGEN_SHOULD_SKIP_THIS
/// @} end of private section
/// @publicsection

// ============================================================
// === main algorithm
// NOLINTBEGIN(readability-function-cognitive-complexity)
auto
almostEqual(const Geometry &geomA, const Geometry &geomB, double tolerance,
            EqualityStrictness strictness) -> bool
{
  if (geomA.geometryTypeId() != geomB.geometryTypeId() &&
      geomA.numGeometries() == geomB.numGeometries()) {
    return false;
  }

  if (geomA.isEmpty()) {
    return true;
  }

  bool out;
  if (geomA.numGeometries() > 1) {
    if (strictness & EqualityStrictness::SubGeomOrdered) {
      out = compareSubGeometryOrdered(geomA, geomB, tolerance, strictness);
    } else {
      out = compareSubGeometryNonOrdered(geomA, geomB, tolerance, strictness);
    }

  } else { // no sub geometries

    if (strictness & EqualityStrictness::CheckCoverOrPoint) {
      out = algorithm::covers3D(geomA, geomB);

    } else if (hasSubPart(geomB)) {
      if (strictness & EqualityStrictness::SubPartOrdered) {
        out = compareSubPartOrdered(geomA, geomB, tolerance, strictness);
      } else {
        out = compareSubPartNonOrdered(geomA, geomB, tolerance, strictness);
      }

    } else {
      detail::GetPointsVisitor getPointsA;
      detail::GetPointsVisitor getPointsB;
      geomA.accept(getPointsA);
      geomB.accept(getPointsB);

      if (getPointsA.points.size() != getPointsB.points.size()) {
        out = false;
      } else {
        if (strictness & EqualityStrictness::InternalPointOrdered) {
          out = comparePointsOrdered(getPointsA, getPointsB, tolerance);
        } else if (strictness & EqualityStrictness::InternalPointShifted) {
          out = comparePointsShifted(getPointsA, getPointsB, tolerance);
        } else if (strictness & EqualityStrictness::InternalPointInverted) {
          out = comparePointsInverted(getPointsA, getPointsB, tolerance);
        } else {
          out = comparePointsNonOrdered(getPointsA, getPointsB, tolerance);
        }
      }
    }
  }

  return out;
}
// NOLINTEND(readability-function-cognitive-complexity)
} // namespace SFCGAL::algorithm
