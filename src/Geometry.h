// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_GEOMETRY_H_
#define SFCGAL_GEOMETRY_H_

#include "SFCGAL/config.h"

#include <boost/endian/conversion.hpp>
#include <boost/shared_ptr.hpp>

#include <memory>
#include <sstream>
#include <string>

#include <boost/assert.hpp>

namespace CGAL {
class Object;
}

/**
 * SRID type
 */
using srid_t = uint32_t;

namespace SFCGAL {

/// @{
/// @privatesection
class Geometry;
class Point;
class LineString;
class Polygon;
class GeometryCollection;
class MultiPoint;
class MultiLineString;
class MultiPolygon;

class Triangle;
class TriangulatedSurface;
class PolyhedralSurface;

// not SFA, appears in GML/CityGML
class Solid;
// not SFA, appears in GML/CityGML
class MultiSolid;

// Curves
class NURBSCurve;

class Envelope;

class GeometryVisitor;
class ConstGeometryVisitor;

const uint32_t wkbSRID = 0x20000000;
const uint32_t wkbM    = 0x40000000;
const uint32_t wkbZ    = 0x80000000;
/// @} end of private section

// NOLINTBEGIN(performance-enum-size)
/**
 * [OGC/SFA]8.2.3 "A common list of codes for geometric types"
 *
 * @todo solid and triangles as non OGC/SFA geometric types?
 * @warning codes for abstract classes and unimplemented classes are hidden
 * @warning code values have are important for WKB
 */
enum GeometryType {
  // TYPE_GEOMETRY = 0, //abstract
  TYPE_POINT              = 1,
  TYPE_LINESTRING         = 2,
  TYPE_POLYGON            = 3,
  TYPE_MULTIPOINT         = 4,
  TYPE_MULTILINESTRING    = 5,
  TYPE_MULTIPOLYGON       = 6,
  TYPE_GEOMETRYCOLLECTION = 7,
  // TYPE_CIRCULARSTRING = 8, // not yet supported
  // TYPE_COMPOUNDCURVE = 9, // not yet supported
  // TYPE_CURVEPOLYGON = 10, // not yet supported
  // TYPE_MULTICURVE = 11, //abstract
  // TYPE_MULTISURFACE = 12, //abstract
  // TYPE_CURVE = 13, // abstract
  // TYPE_SURFACE = 14, //abstract
  TYPE_POLYHEDRALSURFACE   = 15,
  TYPE_TRIANGULATEDSURFACE = 16,
  TYPE_TRIANGLE            = 17,
  // TYPE_CIRCLE = 18,
  // TYPE_GEODESICSTRING = 19,
  // TYPE_ELLIPTICALCURVE = 20,
  TYPE_NURBSCURVE = 21,
  // TYPE_CLOTHOID = 22,
  // TYPE_SPIRALCURVE = 23,

  //-- not official codes
  TYPE_SOLID      = 101,
  TYPE_MULTISOLID = 102,
  // AffinePlacement 102 1102
};

/**
 * @brief coordinate types (XY, XYZ, XYM, etc.)
 * @see SFA 2.8.3 LineStringZ = 1003 ( coordinateType + geometryType)
 */
enum CoordinateType {
  COORDINATE_XY   = 0,
  COORDINATE_XYZ  = 1000,
  COORDINATE_XYM  = 2000,
  COORDINATE_XYZM = 3000
};
// NOLINTEND(performance-enum-size)

/**
 * @brief OGC/SFA based Geometry abstract class
 */
class SFCGAL_API Geometry {
public:
  /**
   * @brief Default constructor.
   */
  Geometry();

  /**
   * @brief Copy constructor.
   */
  Geometry(const Geometry &) = default;

  /**
   * @brief Copy assignemnt operator.
   */
  auto
  operator=(const Geometry &other) -> Geometry & = default;

  /**
   * @brief Destructor.
   */
  virtual ~Geometry() = default;

  /**
   * @brief Get a deep copy of the geometry
   */
  [[nodiscard]] auto
  clone() const -> std::unique_ptr<Geometry>
  {
    return std::unique_ptr<Geometry>(this->cloneImpl());
  }

  /**
   * @brief [OGC/SFA]returns the geometry type
   * @warning use CamelCase (LineString, not LINESTRING)
   */
  [[nodiscard]] virtual auto
  geometryType() const -> std::string = 0;
  /**
   * @brief Returns a code corresponding to the type
   * @warning not standard
   */
  [[nodiscard]] virtual auto
  geometryTypeId() const -> GeometryType = 0;

  /**
   * [OGC/SFA]Dimension of the Geometry ( 0 : punctual, 1 : curve, ...)
   * @warning empty geometries provide the dimension corresponding to the object
   */
  [[nodiscard]] virtual auto
  dimension() const -> int = 0;
  /**
   * [OGC/SFA]returns the dimension of the coordinates
   * @pre suppose no mix of 2D/3D coordinates
   */
  [[nodiscard]] virtual auto
  coordinateDimension() const -> int = 0;
  /**
   * [OGC/SFA]test if geometry is empty
   */
  [[nodiscard]] virtual auto
  isEmpty() const -> bool = 0;

  /**
   * [OGC/SFA]test if geometry is 3d
   * @pre suppose no mix of 2D/3D coordinates
   */
  [[nodiscard]] virtual auto
  is3D() const -> bool = 0;
  /**
   * [OGC/SFA]test if geometry is measured (has an m)
   * @pre suppose no mix of M/!M points
   */
  [[nodiscard]] virtual auto
  isMeasured() const -> bool = 0;

  /**
   * @brief Drops the z coordinate of the geometry
   * @return TRUE if a Z value was present and has been removed
   * @pre suppose no mix of 2D/3D coordinates
   */
  virtual auto
  dropZ() -> bool = 0;

  /**
   * @brief Drops the m coordinate of the geometry
   * @return TRUE if a M value was present and has been removed
   * @pre suppose no mix of M/!M points
   */
  virtual auto
  dropM() -> bool = 0;

  /**
   * @brief Swaps the x and y coordinates of the geometry
   */
  virtual auto
  swapXY() -> void = 0;

  /**
   * @brief Determines the coordinate dimension of a geometry
   *
   * @return CoordinateType The coordinate dimension (XY, XYZ, XYM, XYZM)
   */
  [[nodiscard]] auto
  getCoordinateType() const -> CoordinateType
  {
    bool hasZ = is3D();
    bool hasM = isMeasured();

    if (hasZ) {
      return hasM ? CoordinateType::COORDINATE_XYZM
                  : CoordinateType::COORDINATE_XYZ;
    }
    return hasM ? CoordinateType::COORDINATE_XYM
                : CoordinateType::COORDINATE_XY;
  }

  // virtual bool isSimple() const = 0 ;

  /**
   * Force the state of the validity flag. The validity flag allows to bypass
   * validity checks If the flag is true, it means the geometry is considered
   * valid If the flag is false, it means the validity state of the geometry is
   * unknown The flag is only changed for this geometry and not the internal
   * geometries.
   * @see propagateValidityFlag
   */
  void
  forceValidityFlag(bool validity);

  /** Returns the validity flag */
  [[nodiscard]] auto
  hasValidityFlag() const -> bool;

  /**
   * [OGC/SFA]returns the WKT string
   * @param numDecimals extension specify fix precision output
   */
  [[nodiscard]] auto
  asText(const int &numDecimals = -1) const -> std::string;

  /**
   * [OGC/SFA]returns the WKB string
   */
  [[nodiscard]] auto
  asWkb(boost::endian::order wkbOrder = boost::endian::order::native,
        bool                 asHex    = false) const -> std::string;
  /**
   * [OGC/SFA]Returns a polygon representing the BBOX of the geometry
   * @todo In order to adapt to 3D, would be better to define an "Envelope
   * type", otherway would lead to Polygon and PolyhedralSurface
   */
  // std::unique_ptr< Geometry > envelope() const = 0 ;
  [[nodiscard]] auto
  envelope() const -> Envelope;

  /**
   * @brief [OGC/SFA]Returns the boundary of the geometry
   */
  [[nodiscard]] virtual auto
  boundary() const -> std::unique_ptr<Geometry>;

  /**
   * @brief Computes the distance to an other geometry
   */
  [[nodiscard]] auto
  distance(const Geometry &other) const -> double;
  /**
   * @brief Computes the 3D distance to an other geometry
   */
  [[nodiscard]] auto
  distance3D(const Geometry &other) const -> double;

  //-- helpers

  /**
   * @brief round the geometry with a corresponding scale factor
   * @param scale the scale factor (1 corresponds to the nearest integer, 1000
   * to a 0.001 tolerance)
   */
  void
  round(const long &scale = 1);

  /**
   * Equality operator
   * @todo only compare coordinate points
   * @pre the two geometries must be valid
   * @param other the other geometry to compare with
   * @param tolerance allowed tolerance
   */
  [[nodiscard]] virtual auto
  almostEqual(const Geometry &other, double tolerance) const -> bool;

  /**
   * @brief [OGC/SFA]Gets the number of geometries in a collection of geometries
   * @warning 1 for Point, LineString, Polygon, PolyhedralSurface, Triangle,
   * TriangulatedSurface
   */
  [[nodiscard]] virtual auto
  numGeometries() const -> size_t;
  /**
   * @brief [OGC/SFA]Returns the n-th geometry
   * @warning *this for Point, LineString, Polygon, PolyhedralSurface, Triangle,
   * TriangulatedSurface
   * @param idx The zero-based index of the geometry to retrieve. Must be less
   * than the total number of geometries contained in the collection.
   */
  [[nodiscard]] virtual auto
  geometryN(size_t const &idx) const -> const Geometry &;
  /**
   * @brief [OGC/SFA]Returns the n-th geometry
   * @warning *this for Point, LineString, Polygon, PolyhedralSurface, Triangle,
   * TriangulatedSurface
   * @param idx The zero-based index of the geometry to retrieve. Must be less
   * than the total number of geometries contained in the collection.
   */
  virtual auto
  geometryN(size_t const &idx) -> Geometry &;

  /**
   * @brief [OGC/SFA]Sets the n-th geometry, starting at zero
   * @warning Does nothing for Point, LineString, Polygon, Triangle
   * TriangulatedSurface
   *
   * @param geometry A unique pointer to the new Geometry object that will
   * replace the existing geometry at index `n`. Ownership of the geometry is
   * transferred to the collection.
   * @param idx The zero-based index of the geometry to set. Must be less than
   * the total number of geometries contained in the collection.
   */
  virtual void
  setGeometryN(const Geometry &geometry, size_t const &idx);

  /**
   * @brief [OGC/SFA]Sets the n-th geometry, starting at zero
   * The ownership of the geometry is taken. The caller is not responsible
   * anymore of its deallocation.
   * @warning *this for GeometryCollection, PolyhedralSurface,
   * TriangulatedSurface
   *
   * @param geometry A unique pointer to the new Geometry object that will
   * replace the existing geometry at index `n`. Ownership of the geometry is
   * transferred to the collection.
   * @param idx The zero-based index of the geometry to set. Must be less than
   * the total number of geometries contained in the collection.
   *
   * @deprecated The unique_ptr version should be used instead
   */
  virtual void
  setGeometryN(Geometry *geometry, size_t const &idx);

  /**
   * @brief [OGC/SFA]Sets the n-th geometry, starting at zero
   * @warning *this for GeometryCollection, PolyhedralSurface,
   * TriangulatedSurface
   *
   * @param geometry A unique pointer to the new Geometry object that will
   * replace the existing geometry at index `n`. Ownership of the geometry is
   * transferred to the collection.
   * @param idx The zero-based index of the geometry to set. Must be less than
   * the total number of geometries contained in the collection.
   */
  virtual void
  setGeometryN(std::unique_ptr<Geometry> geometry, size_t const &idx);

  /**
   * @brief Tests if geometry is of "Derived" type given as template parameter
   * @warning not optimized (slow with dynamic_cast)
   */
  template <typename Derived>
  [[nodiscard]] auto
  is() const -> bool
  {
    return dynamic_cast<Derived const *>(this) != NULL;
  }

  /**
   * @brief Downcast to a "Derived" class
   * @warning performs check if boost assertions are enabled
   * @pre The cast must be doable
   */
  template <typename Derived>
  auto
  as() const -> const Derived &
  {
    BOOST_ASSERT(is<Derived>());
    return *static_cast<Derived const *>(this);
  }
  /**
   * @brief Downcast to a "Derived" class
   * @warning performs check if boost assertions are enabled
   * @pre The cast must be doable
   */
  template <typename Derived>
  auto
  as() -> Derived &
  {
    BOOST_ASSERT(is<Derived>());
    return *static_cast<Derived *>(this);
  }

  /**
   * @brief [visitor]dispatch visitor
   */
  virtual void
  accept(GeometryVisitor &visitor) = 0;
  /**
   * @brief [visitor]dispatch visitor
   */
  virtual void
  accept(ConstGeometryVisitor &visitor) const = 0;

  /**
   * Serializer
   */
  template <class Archive>
  void
  serialize(Archive & /*ar*/, const unsigned int /*version*/)
  {
  }

  /**
   * @brief Computes centroid of this geometry
   */
  [[nodiscard]] auto
  centroid() const -> Point;

  /**
   * @brief Computes 3D centroid of this geometry
   */
  [[nodiscard]] auto
  centroid3D() const -> Point;

protected:
  bool validityFlag_ = false;

private:
  [[nodiscard]] virtual auto
  cloneImpl() const -> Geometry * = 0;
};

/**
 * @brief Equality operator for geometries.
 *
 * Compares two Geometry objects for equality.
 *
 * @param[in] geomA The first Geometry object to compare.
 * @param[in] geomB The second Geometry object to compare.
 * @return true if the two Geometry objects have identical coordinate points,
 * false otherwise
 *
 * @pre Both Geometry objects must be valid.
 * @todo Extend comparison to include more than just coordinate points if
 * needed.
 */
SFCGAL_API auto
operator==(const Geometry &geomA, const Geometry &geomB) -> bool;

/**
 * @brief Downcasts a std::unique_ptr Geomety to a derived class.
 *
 * @tparam Derived The derived class to cast to.
 *
 * @param geometry A unique_ptr to the Geometry object. Ownership is released
 * and transferred to the resulting unique_ptr of class @p Derived.
 *
 * @return A std::unique_ptr of class @p Derived that takes ownership of the
 * casted pointer.
 */
template <typename Derived>
SFCGAL_API auto inline geom_unique_ptr_as(std::unique_ptr<Geometry> &&geometry)
    -> std::unique_ptr<Derived>
{
  BOOST_ASSERT(geometry->is<Derived>());
  return std::unique_ptr<Derived>(static_cast<Derived *>(geometry.release()));
}

/// @{
/// @privatesection
#ifndef DOXYGEN_SHOULD_SKIP_THIS

/**
 * @brief Base class that implements covariant cloning with CRTP pattern
 */
template <typename Derived, typename Base>
class SFCGAL_API GeometryImpl : public Base {
  GeometryImpl() = default;

public:
  [[nodiscard]] auto
  clone() const -> std::unique_ptr<Derived>
  {
    return std::unique_ptr<Derived>(static_cast<Derived *>(this->cloneImpl()));
  }

private:
  [[nodiscard]] auto
  cloneImpl() const -> GeometryImpl * override
  {
    return new Derived(*static_cast<const Derived *>(this));
  }
  friend Derived;
};

template <typename T>
class AbstractMethod {};

template <typename Derived, typename Base>
class SFCGAL_API GeometryImpl<AbstractMethod<Derived>, Base> : public Base {
public:
  virtual ~GeometryImpl() = default;

private:
  virtual auto
  cloneImpl() const -> GeometryImpl * = 0;
};

#endif // ifndef DOXYGEN_SHOULD_SKIP_THIS
/// @} end of private section
/// @publicsection

} // namespace SFCGAL

#endif
