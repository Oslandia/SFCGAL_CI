// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_LINESTRING_H_
#define SFCGAL_LINESTRING_H_

#include <memory>
#include <vector>

#include <boost/assert.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/unique_ptr.hpp>

#include "SFCGAL/DereferenceIterator.h"
#include "SFCGAL/Point.h"

#include <CGAL/Polygon_2.h>

namespace SFCGAL {

/**
 * A LineString in SFA
 */
class SFCGAL_API LineString : public GeometryImpl<LineString, Geometry> {
public:
  /// @brief Iterator type for linestring points
  using iterator =
      DereferenceIterator<std::vector<std::unique_ptr<Point>>::iterator>;
  /// @brief Const iterator type for linestring points
  using const_iterator =
      DereferenceIterator<std::vector<std::unique_ptr<Point>>::const_iterator>;

  /**
   * Empty LineString constructor
   */
  LineString();
  /**
   * Constructor with a point vector
   * @param points Vector of points to initialize the linestring
   */
  LineString(const std::vector<Point> &points);
  /**
   * LineString constructor
   * @param startPoint The starting point of the linestring
   * @param endPoint The ending point of the linestring
   */
  LineString(const Point &startPoint, const Point &endPoint);
  /**
   * Copy constructor
   * @param other The linestring to copy from
   */
  LineString(LineString const &other);

  /**
   * assign operator
   * @param other The linestring to assign from
   * @return Reference to this linestring
   */
  auto
  operator=(LineString other) -> LineString &;

  /**
   * destructor
   */
  ~LineString() override;

  //-- SFCGAL::Geometry
  /// @brief Get the geometry type as string
  /// @return "LineString"
  [[nodiscard]] auto
  geometryType() const -> std::string override;
  //-- SFCGAL::Geometry
  /// @brief Get the geometry type identifier
  /// @return TYPE_LINESTRING
  [[nodiscard]] auto
  geometryTypeId() const -> GeometryType override;
  //-- SFCGAL::Geometry
  /// @brief Get the dimension of the linestring
  /// @return 1 (linestrings are 1-dimensional)
  [[nodiscard]] auto
  dimension() const -> int override;
  //-- SFCGAL::Geometry
  /// @brief Get the coordinate dimension
  /// @return Number of coordinates per point
  [[nodiscard]] auto
  coordinateDimension() const -> int override;
  //-- SFCGAL::Geometry
  /// @brief Check if the linestring is empty
  /// @return true if empty, false otherwise
  [[nodiscard]] auto
  isEmpty() const -> bool override;
  //-- SFCGAL::Geometry
  /// @brief Check if the linestring has 3D coordinates
  /// @return true if 3D, false otherwise
  [[nodiscard]] auto
  is3D() const -> bool override;
  //-- SFCGAL::Geometry
  /// @brief Check if the linestring has measured coordinates
  /// @return true if measured, false otherwise
  [[nodiscard]] auto
  isMeasured() const -> bool override;

  /// @brief Drop Z coordinate from all points
  /// @return true if Z was dropped, false otherwise
  auto
  dropZ() -> bool override;

  /// @brief Drop M coordinate from all points
  /// @return true if M was dropped, false otherwise
  auto
  dropM() -> bool override;

  /// @brief Swap X and Y coordinates of all points
  auto
  swapXY() -> void override;

  /**
   * remove all points from the LineString
   */
  void
  clear();

  /**
   * reverse LineString orientation
   */
  void
  reverse();

  /**
   * [SFA/OGC]Returns the number of points
   * @return Number of points in the linestring
   */
  [[nodiscard]] auto
  numPoints() const -> size_t
  {
    return _points.size();
  }
  /**
   * Returns the number of segments
   * @return Number of segments in the linestring
   * @warning not standard, returns zero if LineString contains only one point
   */
  [[nodiscard]] auto
  numSegments() const -> size_t;

  /**
   * [SFA/OGC]Returns the n-th point
   * @param n The index of the point to get
   * @return Const reference to the nth point
   */
  [[nodiscard]] auto
  pointN(size_t const &n) const -> const Point &
  {
    BOOST_ASSERT(n < numPoints());
    return *_points[n];
  }
  /**
   * [SFA/OGC]Returns the n-th point
   * @param n The index of the point to get
   * @return Reference to the nth point
   */
  auto
  pointN(size_t const &n) -> Point &
  {
    BOOST_ASSERT(n < numPoints());
    return *_points[n];
  }

  /**
   * [SFA/OGC]Returns the first point
   * @return Const reference to the first point
   */
  [[nodiscard]] auto
  startPoint() const -> const Point &
  {
    return *_points.front();
  }
  /**
   * [SFA/OGC]Returns the first point
   * @return Reference to the first point
   */
  auto
  startPoint() -> Point &
  {
    return *_points.front();
  }

  /**
   * @brief [SFA/OGC]Returns the last point
   * @return Const reference to the last point
   */
  [[nodiscard]] auto
  endPoint() const -> const Point &
  {
    return *_points.back();
  }
  /**
   * @brief [SFA/OGC]Returns the last point
   * @return Reference to the last point
   */
  auto
  endPoint() -> Point &
  {
    return *_points.back();
  }

  /**
   * @brief Appends a Point to the LineString
   *
   * @param point The Point object to append.
   */
  void
  addPoint(const Point &point)
  {
    addPoint(point.clone());
  }
  /**
   * @brief Appends a Point to the LineString and takes ownership
   *
   * @param point A raw pointer to the Point object to append. Ownership of
   * this point is moved into the LineString.
   * @deprecated The unique_ptr version should be used instead
   */
  void
  addPoint(Point *point)
  {
    addPoint(std::unique_ptr<Point>(point));
  }
  /**
   * @brief Appends a Point to the LineString
   *
   * @param point A unique pointer to the Point object to append. Ownership of
   * this point is moved into the LineString.
   */
  void
  addPoint(std::unique_ptr<Point> point)
  {
    _points.push_back(std::move(point));
  }

  //-- methods

  /**
   * @brief test if the LineString is closed
   * @return true if linestring is closed, false otherwise
   */
  [[nodiscard]] auto
  isClosed() const -> bool;

  /**
   * closes the LineString
   */
  void
  closes()
  {
    _points.push_back(_points.front()->clone());
  }

  //-- iterators

  /**
   * @brief Get iterator to beginning of points
   * @return Iterator to first point
   */
  auto
  begin() -> iterator
  {
    return dereference_iterator(_points.begin());
  }
  /**
   * @brief Get const iterator to beginning of points
   * @return Const iterator to first point
   */
  [[nodiscard]] auto
  begin() const -> const_iterator
  {
    return dereference_iterator(_points.begin());
  }

  /**
   * @brief Get iterator to end of points
   * @return Iterator to past-the-end
   */
  auto
  end() -> iterator
  {
    return dereference_iterator(_points.end());
  }
  /**
   * @brief Get const iterator to end of points
   * @return Const iterator to past-the-end
   */
  [[nodiscard]] auto
  end() const -> const_iterator
  {
    return dereference_iterator(_points.end());
  }

  //-- optimization

  /**
   * @brief Reserve space for points
   * @param n Number of points to reserve space for
   */
  void
  reserve(const size_t &n);

  /**
   * Const iterator to 2D points
   * TODO: replace by boost::tranform_iterator ?
   */
  class Point_2_const_iterator
      : public boost::iterator_facade<Point_2_const_iterator,
                                      Kernel::Point_2 const,
                                      boost::bidirectional_traversal_tag> {
  public:
    /// @brief Default constructor
    Point_2_const_iterator() = default;
    /// @brief Constructor from const_iterator
    /// @param it Iterator to construct from
    explicit Point_2_const_iterator(const_iterator it) : it_(it) {}
    // Point_2_const_iterator( const Point_2_const_iterator<K>& other ) :
    // it_(other.it_) {}
  private:
    friend class boost::iterator_core_access;
    void
    increment()
    {
      it_++;
    }
    void
    decrement()
    {
      it_--;
    }
    auto
    equal(const Point_2_const_iterator &other) const -> bool
    {
      return this->it_ == other.it_;
    }
    auto
    dereference() const -> const Kernel::Point_2 &
    {
      p_ = it_->toPoint_2();
      return p_;
    }
    mutable Kernel::Point_2 p_;
    const_iterator          it_;
  };
  /**
   * @brief Get iterator to beginning of 2D points
   * @return Iterator to first 2D point
   */
  [[nodiscard]] auto
  points_2_begin() const -> Point_2_const_iterator
  {
    return Point_2_const_iterator(begin());
  }
  /**
   * @brief Get iterator to end of 2D points
   * @return Iterator to past-the-end of 2D points
   */
  [[nodiscard]] auto
  points_2_end() const -> Point_2_const_iterator
  {
    return Point_2_const_iterator(end());
  }
  /**
   * @brief Get pair of iterators for 2D points
   * @return Pair of begin and end iterators for 2D points
   */
  [[nodiscard]] auto
  points_2() const -> std::pair<Point_2_const_iterator, Point_2_const_iterator>
  {
    return std::make_pair(points_2_begin(), points_2_end());
  }

  /**
   * Const iterator to 3D points
   * TODO: replace by boost::tranform_iterator ?
   */
  class Point_3_const_iterator
      : public boost::iterator_facade<Point_3_const_iterator,
                                      Kernel::Point_3 const,
                                      boost::bidirectional_traversal_tag> {
  public:
    /// @brief Default constructor
    Point_3_const_iterator() = default;
    /// @brief Constructor from const_iterator
    /// @param it Iterator to construct from
    explicit Point_3_const_iterator(const_iterator it) : it_(it) {}
    /// @brief Copy constructor
    /// @param other Point_3_const_iterator to copy from
    Point_3_const_iterator(const Point_3_const_iterator &other) : it_(other.it_)
    {
    }

  private:
    friend class boost::iterator_core_access;
    void
    increment()
    {
      it_++;
    }
    void
    decrement()
    {
      it_--;
    }
    auto
    equal(const Point_3_const_iterator &other) const -> bool
    {
      return this->it_ == other.it_;
    }
    auto
    dereference() const -> const Kernel::Point_3 &
    {
      p_ = it_->toPoint_3();
      return p_;
    }
    mutable Kernel::Point_3 p_;
    const_iterator          it_;
  };

  /**
   * @brief Get iterator to beginning of 3D points
   * @return Iterator to first 3D point
   */
  [[nodiscard]] auto
  points_3_begin() const -> Point_3_const_iterator
  {
    return Point_3_const_iterator(begin());
  }
  /**
   * @brief Get iterator to end of 3D points
   * @return Iterator to past-the-end of 3D points
   */
  [[nodiscard]] auto
  points_3_end() const -> Point_3_const_iterator
  {
    return Point_3_const_iterator(end());
  }
  /**
   * @brief Get pair of iterators for 3D points
   * @return Pair of begin and end iterators for 3D points
   */
  [[nodiscard]] auto
  points_3() const -> std::pair<Point_3_const_iterator, Point_3_const_iterator>
  {
    return std::make_pair(points_3_begin(), points_3_end());
  }

  /**
   * @brief Convert to CGAL::Polygon_2
   * @param fixOrientation Force exterior ring orientation (counter clockwise)
   * @return CGAL 2D polygon representation
   */
  [[nodiscard]] auto
  toPolygon_2(bool fixOrientation = true) const -> CGAL::Polygon_2<Kernel>;

  //-- visitors

  //-- SFCGAL::Geometry
  /// @brief Accept a geometry visitor
  /// @param visitor Visitor to accept
  void
  accept(GeometryVisitor &visitor) override;
  //-- SFCGAL::Geometry
  /// @brief Accept a const geometry visitor
  /// @param visitor Const visitor to accept
  void
  accept(ConstGeometryVisitor &visitor) const override;

  /**
   * @brief Serializer
   * @param ar Archive for serialization
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int /*version*/)
  {
    ar &boost::serialization::base_object<Geometry>(*this);
    ar & _points;
  }

private:
  std::vector<std::unique_ptr<Point>> _points;

  void
  swap(LineString &other) noexcept
  {
    std::swap(_points, other._points);
  }
};

} // namespace SFCGAL

#endif
