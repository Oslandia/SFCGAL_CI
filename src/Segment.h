// Copyright (c) 2025-2025, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef _SFCGAL_SEGMENT_H_
#define _SFCGAL_SEGMENT_H_

#include "SFCGAL/Kernel.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/export.h"
#include "SFCGAL/numeric.h"
#include <type_traits>
#include <variant>

namespace SFCGAL {

/**
 * @class Segment
 * @brief Represents a line segment in space
 *
 * This class provides methods for working with line segments and
 * performing calculations like distance and interpolation. Unlike LineString
 * which inherits from Geometry, this is a lightweight primitive.
 */
class SFCGAL_API Segment {
public:
  /**
   * @brief Default constructor creating an empty segment
   */
  Segment();

  /**
   * @brief Constructor with two points
   * @param p1 The first endpoint of the segment
   * @param p2 The second endpoint of the segment
   * @throws Exception if points have inconsistent dimensions
   */
  Segment(const Point &p1, const Point &p2);

  /**
   * @brief Constructor with CGAL points (2D or 3D)
   * @param p1 The first endpoint of the segment
   * @param p2 The second endpoint of the segment
   * @note For CGAL points, dimensional consistency is guaranteed by the type
   * system
   */
  template <
      typename PointType,
      typename = std::enable_if_t<std::is_same_v<PointType, Kernel::Point_2> ||
                                  std::is_same_v<PointType, Kernel::Point_3>>>
  Segment(const PointType &p1, const PointType &p2) : _source(p1), _target(p2)
  {
  }

  /**
   * @brief Constructor with CGAL segment (2D or 3D)
   * @param segment A CGAL segment
   * @note For CGAL segments, dimensional consistency is guaranteed by the type
   * system
   */
  template <typename SegmentType,
            typename = std::enable_if_t<
                std::is_same_v<SegmentType, Kernel::Segment_2> ||
                std::is_same_v<SegmentType, Kernel::Segment_3>>>
  explicit Segment(const SegmentType &segment)
      : _source(segment.source()), _target(segment.target())
  {
  }

  /**
   * @brief Copy constructor
   */
  Segment(const Segment &other) = default;

  /**
   * @brief Assignment operator
   */
  auto
  operator=(const Segment &other) -> Segment & = default;

  /**
   * @brief Destructor
   */
  ~Segment() = default;

  /**
   * @brief Gets the first endpoint
   * @return The first endpoint
   */
  auto
  source() const -> const Point &
  {
    return _source;
  }

  /**
   * @brief Gets the second endpoint
   * @return The second endpoint
   */
  auto
  target() const -> const Point &
  {
    return _target;
  }

  /**
   * @brief Sets the first endpoint
   * @param p The new first endpoint
   * @throws Exception if the segment is empty
   * @throws Exception if the new point has inconsistent dimensions with target
   */
  auto
  setSource(const Point &p) -> void;

  /**
   * @brief Sets the second endpoint
   * @param p The new second endpoint
   * @throws Exception if the segment is empty
   * @throws Exception if the new point has inconsistent dimensions with source
   */
  auto
  setTarget(const Point &p) -> void;

  /**
   * @brief Sets both endpoints at once
   * @param source The new first endpoint
   * @param target The new second endpoint
   * @throws Exception if points have inconsistent dimensions
   * @note This is the only way to set points on an empty segment
   */
  auto
  setPoints(const Point &source, const Point &target) -> void;

  /**
   * @brief Sets both CGAL endpoints at once
   * @param source The new first endpoint
   * @param target The new second endpoint
   * @note Type safety ensures dimensional consistency
   */
  template <
      typename PointType,
      typename = std::enable_if_t<std::is_same_v<PointType, Kernel::Point_2> ||
                                  std::is_same_v<PointType, Kernel::Point_3>>>
  auto
  setPoints(const PointType &source, const PointType &target) -> void
  {
    _source = Point(source);
    _target = Point(target);
  }

  /**
   * @brief Reset the segment to empty
   */
  auto
  clear() -> void;

  /**
   * @brief Converts to a CGAL::Segment_2
   * @return A CGAL::Segment_2 representation of this segment
   */
  auto
  toSegment_2() const -> Kernel::Segment_2;

  /**
   * @brief Converts to a CGAL::Segment_3
   * @return A CGAL::Segment_3 representation of this segment
   */
  auto
  toSegment_3() const -> Kernel::Segment_3;

  /**
   * @brief Checks if the segment is empty (either endpoint is empty)
   * @return True if the segment is empty
   */
  auto
  isEmpty() const -> bool;

  /**
   * @brief Checks if the segment has consistent dimension
   * @return True if both endpoints have the same dimension
   * @note This method is redundant since consistency is enforced by
   * constructors and setters
   * @private
   */
  auto
  hasSameDimension() const -> bool;

  /**
   * @brief Checks if the segment is degenerate (zero length)
   * @return True if the segment has zero length
   */
  auto
  isDegenerate() const -> bool;

  /**
   * @brief Returns the squared length of the segment
   * @return The squared length
   */
  auto
  squaredLength() const -> Kernel::FT;

  /**
   * @brief Returns the length of the segment
   * @return The length
   */
  auto
  length() const -> Kernel::FT;

  /**
   * @brief Checks if the segment is in 3D space
   * @return True if either endpoint has Z coordinate
   */
  auto
  is3D() const -> bool;

  /**
   * @brief Checks if the segment is measured
   * @return True if either endpoint has M coordinate
   */
  auto
  isMeasured() const -> bool;

  /**
   * @brief Gets the squared distance from a point to the segment
   * @param p The point (SFCGAL::Point, CGAL::Point_2, or CGAL::Point_3)
   * @return The squared distance (exact calculation)
   * @throws Exception if the segment is empty
   */
  template <typename PointType>
  auto
  squaredDistanceToPoint(const PointType &p) const -> Kernel::FT;

  /**
   * @brief Distance from a point to the segment
   * @param p The point (can be SFCGAL::Point, coordinate pair, or CGAL points)
   * @return The shortest distance from the point to the segment
   * @throws Exception if the segment is empty
   */
  template <typename... Args>
  auto
  distanceToPoint(Args &&...args) const -> double;

  /**
   * @brief Calculates the exact parameter for a projected point on the segment
   * @param p The point to project (SFCGAL::Point, CGAL::Point_2, or
   * CGAL::Point_3)
   * @return Parameter value between 0.0 and 1.0 (clamped)
   * @throws Exception if the segment is empty
   */
  template <typename PointType>
  auto
  exactInterpolationParameter(const PointType &p) const -> Kernel::FT;

  /**
   * @brief Calculates the parameter for a projected point on the segment
   * @param args Point coordinates or point object
   * @return Parameter value between 0.0 (source) and 1.0 (target)
   * @throws Exception if the segment is empty
   */
  template <typename... Args>
  auto
  interpolationParameter(Args &&...args) const -> double;

  /**
   * @brief Gets a point on the segment at given parameter
   * @param t Parameter value between 0.0 (source) and 1.0 (target)
   * @return The interpolated point with M value if points are measured
   * @throws Exception if the segment is empty
   */
  auto
  interpolate(double t) const -> Point;

  /**
   * @brief Checks if a point is on the segment (works in 2D or 3D)
   * @param p The point to check
   * @param tolerance Optional tolerance value
   * @return True if the point lies on the segment
   * @throws Exception if the segment is empty
   */
  template <typename PointType>
  auto
  hasOn(const PointType &p, double tolerance = EPSILON) const -> bool;

  /**
   * @brief Returns the midpoint of the segment
   * @return A point at the middle of the segment
   * @throws Empty point if the segment is empty
   */
  auto
  midpoint() const -> Point;

  /**
   * @brief Reverses the segment direction
   * @throws Exception if the segment is empty
   */
  auto
  reverse() -> void;

private:
  Point _source; ///< First endpoint
  Point _target; ///< Second endpoint

  /**
   * @brief Helper to apply a function depending on the dimensionality
   * @param func2D Function to call for 2D computation
   * @param func3D Function to call for 3D computation
   */
  template <typename Func2D, typename Func3D>
  auto
  applyByDimension(Func2D func2D, Func3D func3D) const
      -> std::invoke_result_t<Func2D>;

  /// Private helper method to calculate parameter from a point
  template <typename SegmentType, typename PointType>
  auto
  calculateParameterFromPoint(const SegmentType &segment,
                              const PointType   &point) const -> Kernel::FT;
};

// Template implementation

template <typename Func2D, typename Func3D>
auto
Segment::applyByDimension(Func2D func2D, Func3D func3D) const
    -> std::invoke_result_t<Func2D>
{
  if (is3D()) {
    return func3D();
  }
  return func2D();
}

template <typename PointType>
auto
Segment::squaredDistanceToPoint(const PointType &p) const -> Kernel::FT
{
  if (isEmpty()) {
    return 0;
  }

  if constexpr (std::is_same_v<PointType, Point>) {
    return applyByDimension(
        [this, &p]() { return squaredDistanceToPoint(p.toPoint_2()); },
        [this, &p]() { return squaredDistanceToPoint(p.toPoint_3()); });
  } else if constexpr (std::is_same_v<PointType, Kernel::Point_2>) {
    const auto segment = toSegment_2();
    if (segment.is_degenerate()) {
      return CGAL::squared_distance(p, segment.source());
    }
    return CGAL::squared_distance(p, segment);
  } else if constexpr (std::is_same_v<PointType, Kernel::Point_3>) {
    const auto segment = toSegment_3();
    if (segment.is_degenerate()) {
      return CGAL::squared_distance(p, segment.source());
    }
    return CGAL::squared_distance(p, segment);
  }
}

template <typename... Args>
auto
Segment::distanceToPoint(Args &&...args) const -> double
{
  if constexpr (sizeof...(Args) == 1) {
    // Single argument: either Point or Point_2 or Point_3
    using FirstArg = std::tuple_element_t<0, std::tuple<Args...>>;
    if constexpr (std::is_same_v<std::decay_t<FirstArg>, Point> ||
                  std::is_same_v<std::decay_t<FirstArg>, Kernel::Point_2> ||
                  std::is_same_v<std::decay_t<FirstArg>, Kernel::Point_3>) {
      auto squared = squaredDistanceToPoint(std::forward<Args>(args)...);
      return std::sqrt(CGAL::to_double(squared));
    }
  } else if constexpr (sizeof...(Args) == 2) {
    // Handle case with x, y coordinates
    Kernel::Point_2 point(std::forward<Args>(args)...);
    return std::sqrt(CGAL::to_double(squaredDistanceToPoint(point)));
  }

  // Default case (should not happen)
  throw std::invalid_argument("Invalid arguments for distanceToPoint");
}

template <typename SegmentType, typename PointType>
auto
Segment::calculateParameterFromPoint(const SegmentType &segment,
                                     const PointType &point) const -> Kernel::FT
{
  const auto v = segment.to_vector();
  const auto w = point - segment.source();
  auto       t = (w * v) / v.squared_length();
  return std::clamp(t, Kernel::FT(0), Kernel::FT(1));
}

template <typename PointType>
auto
Segment::exactInterpolationParameter(const PointType &p) const -> Kernel::FT
{
  if (isEmpty()) {
    return 0;
  }

  if constexpr (std::is_same_v<PointType, Point>) {
    return applyByDimension(
        [this, &p]() { return exactInterpolationParameter(p.toPoint_2()); },
        [this, &p]() { return exactInterpolationParameter(p.toPoint_3()); });
  } else if constexpr (std::is_same_v<PointType, Kernel::Point_2>) {
    const auto segment = toSegment_2();
    if (segment.is_degenerate()) {
      return 0; // Return source point parameter
    }

    return calculateParameterFromPoint(segment, p);
  } else if constexpr (std::is_same_v<PointType, Kernel::Point_3>) {
    const auto segment = toSegment_3();
    if (segment.is_degenerate()) {
      return 0; // Return source point parameter
    }

    const auto line       = segment.supporting_line();
    const auto projection = line.projection(p);

    return calculateParameterFromPoint(segment, projection);
  }
}

template <typename... Args>
auto
Segment::interpolationParameter(Args &&...args) const -> double
{
  if constexpr (sizeof...(Args) == 1) {
    // Single argument: Point
    using FirstArg = std::tuple_element_t<0, std::tuple<Args...>>;
    if constexpr (std::is_same_v<std::decay_t<FirstArg>, Point>) {
      return CGAL::to_double(
          exactInterpolationParameter(std::forward<Args>(args)...));
    }
  } else if constexpr (sizeof...(Args) == 2) {
    // Handle case with x, y coordinates
    Kernel::Point_2 point(std::forward<Args>(args)...);
    return CGAL::to_double(exactInterpolationParameter(point));
  }

  // Default case (should not happen)
  throw std::invalid_argument("Invalid arguments for interpolationParameter");
}

template <typename PointType>
auto
Segment::hasOn(const PointType &p, double tolerance) const -> bool
{
  if (isEmpty()) {
    return false;
  }

  if constexpr (std::is_same_v<PointType, Point>) {
    return applyByDimension(
        [this, &p, tolerance]() { return hasOn(p.toPoint_2(), tolerance); },
        [this, &p, tolerance]() { return hasOn(p.toPoint_3(), tolerance); });
  } else if constexpr (std::is_same_v<PointType, Kernel::Point_2>) {
    const auto segment = toSegment_2();

    if (segment.is_degenerate()) {
      return CGAL::squared_distance(p, segment.source()) <=
             tolerance * tolerance;
    }

    if (segment.has_on(p)) {
      return true;
    }

    return CGAL::squared_distance(p, segment) <= tolerance * tolerance;
  } else if constexpr (std::is_same_v<PointType, Kernel::Point_3>) {
    const auto segment = toSegment_3();

    if (segment.is_degenerate()) {
      return CGAL::squared_distance(p, segment.source()) <=
             tolerance * tolerance;
    }

    if (segment.has_on(p)) {
      return true;
    }

    return CGAL::squared_distance(p, segment) <= tolerance * tolerance;
  }

  return false;
}

} // namespace SFCGAL

#endif // _SFCGAL_SEGMENT_H_
