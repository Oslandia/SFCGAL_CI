// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <cmath>
#include <memory>

// SFCGAL
#include "SFCGAL/Exception.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/algorithm/length.h"
#include "SFCGAL/algorithm/lineSubstring.h"

namespace SFCGAL::algorithm {

namespace {
const double tol = 1.0e-9;

auto
find_position(const LineString &ls, const long N, const double target_length,
              const double offset, const double tol, std::size_t &idx,
              double &frac, bool &on_point, double &len_to_idx) -> Point
{
  BOOST_ASSERT(!(offset < 0.0));
  BOOST_ASSERT(!(target_length < 0.0));
  BOOST_ASSERT(offset < (target_length + tol));

  if (std::fabs(offset - target_length) < tol) {
    // Point(idx) is the target position.

    frac       = 0.0;
    on_point   = true;
    len_to_idx = offset;
    return ls.pointN(idx);
  }

  double cur_length = offset;
  double seg_length = 0.0;
  on_point          = false;

  for (; idx < static_cast<size_t>(N); ++idx) {
    const Point &p = ls.pointN(idx);
    const Point &q = ls.pointN(idx + 1);

    double seg_length_sq = std::pow(CGAL::to_double(p.x() - q.x()), 2.0) +
                           std::pow(CGAL::to_double(p.y() - q.y()), 2.0);
    if (ls.is3D()) {
      seg_length_sq += std::pow(CGAL::to_double(p.z() - q.z()), 2.0);
    }

    seg_length = std::sqrt(seg_length_sq);

    cur_length += seg_length;

    if (std::fabs(cur_length - target_length) < tol) {
      // Adjust idx to be that of the Point coincident
      // with the desired position.

      ++idx;
      on_point = true;

      break;
    }
    if (cur_length > target_length) {
      // We went too far. Subtract seg_length so
      // cur_length is the distance along ls
      // to the idx'th point.
      cur_length -= seg_length;
      break;
    }
  }

  // Return distance to point immediately before
  // desired position.
  len_to_idx = cur_length;

  // Calculate fraction between idx and idx + 1 where
  // the desired position resides.

  frac = 0.0;
  if (!on_point) {
    BOOST_ASSERT(seg_length > tol);
    frac = (target_length - cur_length) / seg_length;
  }

  // Calculate point.

  Point ret;
  if (on_point) {
    ret = ls.pointN(idx);
  } else {
    const Point &p = ls.pointN(idx);
    const Point &q = ls.pointN(idx + 1);

    const Kernel::RT x = p.x() + (frac * (q.x() - p.x()));
    const Kernel::RT y = p.y() + (frac * (q.y() - p.y()));

    Kernel::RT z;
    if (ls.is3D()) {
      z   = p.z() + (frac * (q.z() - p.z()));
      ret = Point(x, y, z);
    } else {
      ret = Point(x, y);
    }

    if (ls.isMeasured()) {
      ret.setM(p.m() + (frac * (q.m() - p.m())));
    }
  }

  return ret;
}

} // namespace

auto
lineSubstring(const LineString &ls, double start, double end)
    -> std::unique_ptr<LineString>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY(ls);

  if (ls.isEmpty()) {
    // Empty line, therefore start and end are
    // irrelevant.
    return std::make_unique<LineString>();
  }

  // Check for out of range start, end.

  if (std::fabs(start) > 1.0) {
    BOOST_THROW_EXCEPTION(Exception(
        "SFCGAL::algorithm::lineSubstring: start value out of range."));
  }

  if (std::fabs(end) > 1.0) {
    BOOST_THROW_EXCEPTION(
        Exception("SFCGAL::algorithm::lineSubstring: end value out of range."));
  }

  // Convert start and end into their equivalent positive values.

  if (start < 0.0) {
    start = 1.0 + start;
  }

  if (end < 0.0) {
    end = 1.0 + end;
  }

  // Check for equal start and end.
  if (std::fabs(start - end) < tol) {
    // start and end are equal, hence return an empty line substring.
    return std::make_unique<LineString>();
  }

  const unsigned long N = static_cast<unsigned long>(ls.numPoints());

  const bool closed = ls.isClosed();

  bool reverse = false;
  if (start > end) {
    if (closed && std::fabs(start - end - 1.0) < tol) {
      // We desire a complement line substring of a closed
      // line which has zero length, hence return empty
      // substring.
      return std::make_unique<LineString>();
    }

    // Swap the start and end positions so that they
    // define a positive range on ls, setting the
    // reverse flag to reverse to indicate that
    // we have done so. Note that reverse being set
    // has different consequences for open and
    // closed lines.

    std::swap(start, end);
    reverse = true;
  } else if (closed && std::fabs(end - start - 1.0) < tol) {
    // The desired line substring is the entire line.
    return std::unique_ptr<LineString>(ls.clone());
  }

  // Retrieve length of the line.
  double len = 0.0;
  if (ls.is3D()) {
    len = algorithm::length3D(ls);
  } else {
    len = algorithm::length(ls);
  }

  // Find Point immediately before/on start position.
  std::size_t start_idx        = 0; // Must initialise first.
  double      start_frac       = 0.0;
  bool        on_start         = false;
  double      len_to_start_idx = 0.0;
  Point       pstart = find_position(ls, N, len * start, 0.0, tol, start_idx,
                                     start_frac, on_start, len_to_start_idx);

  // Find Point immediately before/on end position.
  std::size_t end_idx        = start_idx; // Must initialise first.
  double      end_frac       = 0.0;
  bool        on_end         = false;
  double      len_to_end_idx = 0.0;
  Point       pend =
      find_position(ls, N, len * end, len_to_start_idx, tol, end_idx, end_frac,
                    on_end, len_to_end_idx // This result is not used.
      );

  if (reverse && closed) {
    // For closed lines we always want to follow the
    // direction of the original line. A set reversed
    // flag indicates that we are going to cross the
    // join, and hence we add N to the end_idx and
    // use modulus operation when dereferncing the points
    // to be added to the desired line substring.

    std::swap(start_idx, end_idx);
    std::swap(start_frac, end_frac);
#if defined(_MSC_VER) && (_MSC_VER < 1921)
    // https://developercommunity.visualstudio.com/content/problem/431904/vc-2017-stdswap-is-not-conforming.html
    Point temporary = std::move(pstart);
    pstart          = std::move(pend);
    pend            = std::move(temporary);
#else
    std::swap<Point>(pstart, pend);
#endif
    std::swap(on_start, on_end);
    end_idx += N;
  }

  // Construct the desired line substring.

  std::unique_ptr<LineString> substring = std::make_unique<LineString>();

  // Add start point.

  substring->addPoint(pstart);

  // Add intermediate points.

  bool skipped_duplicate = false;
  for (std::size_t i = start_idx + 1; i <= end_idx; ++i) {
    // Ensure that we don't add a duplicate point to match the
    // one at the join of a closed line if our desired line
    // substring crosses that join.
    if (closed &&
        // Check the required line substring is one that crosses
        // the join.
        reverse &&
        // Check we have not already skipped the duplicate point
        // at the join.
        (!skipped_duplicate) &&
        // check current point is one of the two duplicates at
        // the join, we encounter first.
        (((i % N) == 0) || ((i % N) == (N - 1)))) {
      // Skip the duplicate point. If this was
      // the last point and we ended up with
      // an invalid single point substring, the
      // fact that it is now made invalid is
      // not a concern since if added we would
      // have obtained a zero-length substring,
      // which is also invalid.
      skipped_duplicate = true;
      continue;
    }

    const Point &p = ls.pointN(i % N);

    substring->addPoint(p);
  }

  // Add end point if we have not already.

  if (!on_end) {
    substring->addPoint(pend);
  }

  if (reverse && (!closed)) {
    // Reverse the constructed substring.
    substring->reverse();
  }

  return substring;
}

} // namespace SFCGAL::algorithm
