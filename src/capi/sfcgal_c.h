// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_CAPI_H_
#define SFCGAL_CAPI_H_

#include "SFCGAL/config.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__has_c_attribute)
  #if __has_c_attribute(deprecated)
    #define SFCGAL_DEPRECATED(msg) [[deprecated(msg)]]
  #elif defined(__clang__) || defined(__GNUC__)
    #define SFCGAL_DEPRECATED(msg) __attribute__((deprecated(msg)))
  #elif defined(_MSC_VER)
    #define SFCGAL_DEPRECATED(msg) __declspec(deprecated(msg))
  #else
    #define SFCGAL_DEPRECATED(msg)
  #endif
#elif defined(__has_cpp_attribute)
  #if __has_cpp_attribute(deprecated)
    #define SFCGAL_DEPRECATED(msg) [[deprecated(msg)]]
  #elif defined(__clang__) || defined(__GNUC__)
    #define SFCGAL_DEPRECATED(msg) __attribute__((deprecated(msg)))
  #elif defined(_MSC_VER)
    #define SFCGAL_DEPRECATED(msg) __declspec(deprecated(msg))
  #else
    #define SFCGAL_DEPRECATED(msg)
  #endif
#elif defined(__clang__) || defined(__GNUC__)
  #define SFCGAL_DEPRECATED(msg) __attribute__((deprecated(msg)))
#elif defined(_MSC_VER)
  #define SFCGAL_DEPRECATED(msg) __declspec(deprecated(msg))
#else
  #pragma message("WARNING: deprecated if not supported by this compiler")
  #define SFCGAL_DEPRECATED(msg)
#endif

// TODO : return of errors ! => error handler

/*--------------------------------------------------------------------------------------*
 *
 * Support for SFCGAL::Geometry class hierarchy
 *
 *--------------------------------------------------------------------------------------*/

/**
 * SRID type
 * @ingroup capi
 */
typedef uint32_t srid_t;

/**
 * sfcgal_geometry_t is an opaque pointer type that is used to represent a
 * pointer to SFCGAL::Geometry
 * @ingroup capi
 */
typedef void sfcgal_geometry_t;

/**
 * Geometric types
 * @ingroup capi
 */
typedef enum {
  //     TYPE_GEOMETRY            = 0, //abstract
  SFCGAL_TYPE_POINT              = 1,
  SFCGAL_TYPE_LINESTRING         = 2,
  SFCGAL_TYPE_POLYGON            = 3,
  SFCGAL_TYPE_MULTIPOINT         = 4,
  SFCGAL_TYPE_MULTILINESTRING    = 5,
  SFCGAL_TYPE_MULTIPOLYGON       = 6,
  SFCGAL_TYPE_GEOMETRYCOLLECTION = 7,
  //     TYPE_CIRCULARSTRING      = 8,
  //     TYPE_COMPOUNDCURVE       = 9,
  //     TYPE_CURVEPOLYGON        = 10,
  //     TYPE_MULTICURVE          = 11, //abstract
  //     TYPE_MULTISURFACE        = 12, //abstract
  //     TYPE_CURVE               = 13, //abstract
  //     TYPE_SURFACE             = 14, //abstract
  SFCGAL_TYPE_POLYHEDRALSURFACE   = 15,
  SFCGAL_TYPE_TRIANGULATEDSURFACE = 16,
  SFCGAL_TYPE_TRIANGLE            = 17,

  SFCGAL_TYPE_NURBSCURVE = 21,

  //-- not official codes
  SFCGAL_TYPE_SOLID      = 101,
  SFCGAL_TYPE_MULTISOLID = 102
} sfcgal_geometry_type_t;

/**
 * Set the geometry validation mode
 * @param enabled 1 for enabled, 0 for disabled.
 * @ingroup capi
 * @note obsolete
 * @deprecated It is not used anymore. It has no effect.
 */
SFCGAL_DEPRECATED("It is not used anymore. It has no effect.")
SFCGAL_API void
sfcgal_set_geometry_validation(int enabled);

/**
 * Returns the type of a given geometry
 * @param geom the input geometry
 * @return the type of geom
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_type_t
sfcgal_geometry_type_id(const sfcgal_geometry_t *geom);

/**
 * Returns the type of a given geometry as a string
 * @param geom the input geometry
 * @param[out] type the output buffer
 * @param[out] typeLen the size of the buffer
 * @post type is returned allocated and must be freed by the caller with
 * sfcgal_free_buffer()
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_geometry_type(const sfcgal_geometry_t *geom, char **type,
                     size_t *typeLen);

/**
 * Returns the dimension of a given geometry ( 0 : punctual, 1 : curve, ...)
 * @param geom the input geometry
 * @return the dimension as an integer
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_dimension(const sfcgal_geometry_t *geom);

/**
 * Tests if the given geometry is valid or not
 * @param geom the input geometry
 * @return 1 if the geometry is valid. 0 otherwise.
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_is_valid(const sfcgal_geometry_t *geom);

/**
 * Tests if the given geometry is valid or not
 * And return details in case of invalidity
 * @param[in] geom the input geometry
 * @param[out] invalidity_reason input/output parameter. If non null, a
 * null-terminated string could be allocated and contain reason of the
 * invalidity
 * @param[out] invalidity_location input/output parameter. If non null, a
 * geometry could be allocated and contain the location of the invalidity
 * @return 1 if the geometry is valid. 0 otherwise.
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_is_valid_detail(const sfcgal_geometry_t *geom,
                                char                   **invalidity_reason,
                                sfcgal_geometry_t      **invalidity_location);

/**
 * @deprecated Use sfcgal_geometry_is_valid_detail() instead.
 * @param[in] geom the input geometry
 * @param[out] invalidity_reason invalidity reason
 * @param[out] invalidity_location invalidity location
 * @return 1 if geom is simple. 0 otherwise.
 * @ingroup capi
 */
SFCGAL_DEPRECATED("Use sfcgal_geometry_is_valid_detail instead")
SFCGAL_API int
sfcgal_geometry_is_complexity_detail(const sfcgal_geometry_t *geom,
                                     char                   **invalidity_reason,
                                     sfcgal_geometry_t **invalidity_location);

/**
 * Tests if the given geometry is simple or not
 * @param geom the input geometry
 * @return 1 if geom is simple. 0 otherwise.
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_is_simple(const sfcgal_geometry_t *geom);

/**
 * Tests if the given geometry is simple or not
 * And return details in case of complexity
 * @param geom the input geometry
 * @param complexity_reason input/output parameter. If non null, a
 * null-terminated string could be allocated and contain reason of the
 * complexity
 * @return 1 if geom is simple. 0 otherwise.
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_is_simple_detail(const sfcgal_geometry_t *geom,
                                 char                   **complexity_reason);

/**
 * Tests if the given geometry is 3D or not
 * @param geom the input geometry
 * @return 1 if geom is 3d. 0 otherwise.
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_is_3d(const sfcgal_geometry_t *geom);

/**
 * Tests if the given geometry is measured (has an m) or not
 * @param geom the input geometry
 * @return 1 if geom is measured. 0 otherwise.
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_is_measured(const sfcgal_geometry_t *geom);

/**
 * Tests if the given geometry is empty or not
 * @param geom the input geometry
 * @return 1 if geom is empty. 0 otherwise.
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_is_empty(const sfcgal_geometry_t *geom);

/**
 * Tests if the given geometry is closed or not
 * @param geom the input geometry
 * @return 1 if geom is measured. 0 otherwise.
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_is_closed(const sfcgal_geometry_t *geom);

/**
 * Drops the z coordinate of the geometry
 * @param geom the input geometry
 * @return 1 if a Z value was present and has been removed. 0 otherwise.
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_drop_z(sfcgal_geometry_t *geom);

/**
 * Drops the m coordinate of the geometry
 * @param geom the input geometry
 * @return 1 if a M value was present and has been removed. 0 otherwise.
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_drop_m(sfcgal_geometry_t *geom);

/**
 * Adds a z-dimension to the geometry, initialized to a preset value.
 * Existing Z values remains unchanged.
 * This has no effect on empty geometries.
 * @param geom the input geometry
 * @param defaultZ z-value to use
 * @return 1 if a Z value was added. 0 otherwise.
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_force_z(sfcgal_geometry_t *geom, double defaultZ);

/**
 * Adds a m-dimension to the geometry, initialized to a preset value.
 * Existing M values remains unchanged.
 * This has no effect on empty geometries.
 * @param geom the input geometry
 * @param defaultM m-value to use
 * @return 1 if a M value was added. 0 otherwise.
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_force_m(sfcgal_geometry_t *geom, double defaultM);

/**
 * Swaps the x and y coordinates of the geometry
 * @param geom the input geometry
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_geometry_swap_xy(sfcgal_geometry_t *geom);

/**
 * Returns a deep clone of the given geometry
 * @param geom the input geometry
 * @post returns a pointer to an allocated geometry that must be deallocated by
 * sfcgal_geometry_delete()
 * @return The cloned geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_clone(const sfcgal_geometry_t *geom);

/**
 * Deletes a given geometry
 * @param geom the input geometry
 * @pre the given pointer must have been previously allocated by a creation
 * function
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_geometry_delete(sfcgal_geometry_t *geom);

/**
 * Returns the number of geometries of the given geometry collection
 * @param geometryCollection the input geometry
 * @pre geometry must be a SFCGAL::GeometryCollection.
 * Otherwise, 1 is returned.  For empty geometries 0 is
 * returned.
 * @return the the number of geometries of geometryCollection
 * @ingroup capi
 */
SFCGAL_API size_t
sfcgal_geometry_num_geometries(const sfcgal_geometry_t *geometryCollection);

/**
 * Returns a WKT representation of the given geometry using CGAL exact integer
 * fractions as coordinate values
 * @param geom the input geometry
 * @param[out] buffer the output buffer
 * @param[out] len the size of the @p buffer
 * @post @p buffer is returned allocated and must be freed by the caller with
 *       sfcgal_free_buffer()
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_geometry_as_text(const sfcgal_geometry_t *geom, char **buffer,
                        size_t *len);

/**
 * Returns a WKT representation of the given geometry using floating point
 * coordinate values. Floating point precision can be set via the numDecimals
 * parameter. Setting numDecimals to -1 yields the same result as
 * sfcgal_geometry_as_text()
 * @param geom the input geometry
 * @param numDecimals decimal precision
 * @param[out] buffer the output buffer
 * @param[out] len the size of the @p buffer
 * @post @p buffer is returned allocated and must be freed by the caller with
 *       sfcgal_free_buffer()
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_geometry_as_text_decim(const sfcgal_geometry_t *geom, int numDecimals,
                              char **buffer, size_t *len);

/**
 * Returns a WKB representation of the given geometry
 * @param geom the input geometry
 * @param[out] buffer the output buffer
 * @param[out] len the size of the @p buffer
 * @post @p buffer is returned allocated and must be freed by the caller with
 *       sfcgal_free_buffer()
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_geometry_as_wkb(const sfcgal_geometry_t *geom, char **buffer,
                       size_t *len);

/**
 * Returns a WKB representation as hexadecimal of the given geometry
 * @param geom the input geometry
 * @param[out] buffer the output buffer
 * @param[out] len the size of the @p buffer
 * @post @p buffer is returned allocated and must be freed by the caller with
 *       sfcgal_free_buffer()
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_geometry_as_hexwkb(const sfcgal_geometry_t *geom, char **buffer,
                          size_t *len);

/**
 * Creates a VTK string of the given geometry
 * @param geom the input geometry
 * @param[out] buffer the output buffer
 * @param[out] len the size of the @p buffer
 * @post @p buffer is returned allocated and must be freed by the caller with
 *       sfcgal_free_buffer()
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_geometry_as_vtk(const sfcgal_geometry_t *geom, char **buffer,
                       size_t *len);

/**
 * Creates a VTK file of the given geometry
 * @param geom the input geometry
 * @param filename VTK filename
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_geometry_as_vtk_file(const sfcgal_geometry_t *geom,
                            const char              *filename);

/**
 * Creates a OBJ file of the given geometry
 * The generated OBJ file uses the Z-up convention.
 * @param geom the input geometry
 * @param filename OBJ filename
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_geometry_as_obj_file(const sfcgal_geometry_t *geom,
                            const char              *filename);

/**
 * Creates a OBJ string of the given geometry
 * The generated OBJ string uses the Z-up convention.
 * @param[in] geom the input geometry
 * @param[out] buffer the OBJ buffer
 * @param[out] len the size of @p buffer
 * @post @p buffer is returned allocated and must be freed by the caller with
 *       sfcgal_free_buffer()
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_geometry_as_obj(const sfcgal_geometry_t *geom, char **buffer,
                       size_t *len);

/**
 * Creates a STL string of the given geometry
 * The generated STL string uses the Z-up convention.
 * @param[in] geom the input geometry
 * @param[out] buffer the STL buffer
 * @param[out] len the size of @p buffer
 * @post @p buffer is returned allocated and must be freed by the caller with
 *       sfcgal_free_buffer()
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_geometry_as_stl(const sfcgal_geometry_t *geom, char **buffer,
                       size_t *len);

/**
 * Creates a STL file of the given geometry
 * The generated STL file uses the Z-up convention.
 * @param geom the input geometry
 * @param filename STL filename
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_geometry_as_stl_file(const sfcgal_geometry_t *geom,
                            const char              *filename);

/**
 * Creates an empty point
 * @return an empty SFCGAL::Point
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_point_create();

/**
 * Creates a point from two X and Y coordinates
 * @param x X coordinate
 * @param y Y coordinate
 * @return a SFCGAL::Point with X, Y coordinates
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_point_create_from_xy(double x, double y);

/**
 * Creates a point from three X, Y and Z coordinates
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @return a SFCGAL::Point with X, Y, Z coordinates
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_point_create_from_xyz(double x, double y, double z);

/**
 * Creates a point from three X, Y and M coordinates
 * @param x X coordinate
 * @param y Y coordinate
 * @param m M coordinate
 * @return a SFCGAL::Point with X, Y, M coordinates
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_point_create_from_xym(double x, double y, double m);

/**
 * Creates a point from four X, Y, Z and M coordinates
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @param m M coordinate
 * @return a SFCGAL::Point with X, Y, Z, M coordinates
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_point_create_from_xyzm(double x, double y, double z, double m);

/**
 * Returns the X coordinate of the given SFCGAL::Point
 * @param geom the input geometry
 * @pre the given geometry must be a SFCGAL::Point
 * @pre the given point must not be empty
 * @return the X coordinate of the point
 * @ingroup capi
 */
SFCGAL_API double
sfcgal_point_x(const sfcgal_geometry_t *geom);

/**
 * Returns the Y coordinate of the given SFCGAL::Point
 * @param geom the input geometry
 * @pre the given geometry must be a SFCGAL::Point
 * @pre the given point must not be empty
 * @return the Y coordinate of the point
 * @ingroup capi
 */
SFCGAL_API double
sfcgal_point_y(const sfcgal_geometry_t *geom);

/**
 * Returns the Z coordinate of the given SFCGAL::Point
 * @param geom the input geometry
 * @pre the given geometry must be a SFCGAL::Point
 * @pre the given point must not be empty
 * @post the Z coordinate can value NaN if the given point is 2D only
 * @return the Z coordinate of the point
 * @ingroup capi
 */
SFCGAL_API double
sfcgal_point_z(const sfcgal_geometry_t *geom);

/**
 * Returns the M coordinate of the given SFCGAL::Point
 * @param geom the input geometry
 * @pre the given geometry must be a SFCGAL::Point
 * @pre the given point must not be empty
 * @post the M coordinate can value NaN if the given point has no m
 * @return the M coordinate of the point
 * @ingroup capi
 */
SFCGAL_API double
sfcgal_point_m(const sfcgal_geometry_t *geom);

/**
 * Creates an empty SFCGAL::LineString
 * @return an empty SFCGAL::LineString
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_linestring_create();

/**
 * Returns the number of points of the given SFCGAL::LineString
 * @param linestring the input geometry
 * @pre linestring must be a SFCGAL::LineString
 * @return the number of points of linestring
 * @ingroup capi
 */
SFCGAL_API size_t
sfcgal_linestring_num_points(const sfcgal_geometry_t *linestring);

/**
 * Returns the ith point of a given SFCGAL::LineString
 * @param linestring the input SFCGAL::LineString
 * @param i is the point index in the SFCGAL::LineString
 * @pre linestring must be a SFCGAL::LineString
 * @pre i >= and i < sfcgal_linestring_num_points
 * @post the returned SFCGAL::Point is not writable and must not be deallocated
 * by the caller
 * @return the ith point of linestring
 * @ingroup capi
 */
SFCGAL_API const sfcgal_geometry_t *
sfcgal_linestring_point_n(const sfcgal_geometry_t *linestring, size_t i);

/**
 * Adds a point to a SFCGAL::LineString
 * @param linestring is the SFCGAL::LineString where the SFCGAL::Point has to be
 * added to
 * @param point is the SFCGAL::Point to add to the given SFCGAL::LineString
 * @pre i >= and i < sfcgal_linestring_num_points
 * @post the ownership of SFCGAL::Point is taken by the function
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_linestring_add_point(sfcgal_geometry_t *linestring,
                            sfcgal_geometry_t *point);

/**
 * Creates an empty SFCGAL::Triangle
 * @return an empty SFCGAL::Triangle
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_triangle_create();

/**
 * Creates a SFCGAL::Triangle from three given SFCGAL::Point
 * @param pta first input point
 * @param ptb second input point
 * @param ptc third input point
 * @pre pta must be a SFCGAL::Point
 * @pre ptb must be a SFCGAL::Point
 * @pre ptc must be a SFCGAL::Point
 * @post the ownership of the three points are not taken. The caller is still
 * responsible of their deallocation
 * @return a SFCGAL::Triangle from SFCGAL::Point pta, ptb and ptc
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_triangle_create_from_points(const sfcgal_geometry_t *pta,
                                   const sfcgal_geometry_t *ptb,
                                   const sfcgal_geometry_t *ptc);

/**
 * Returns one the SFCGAL::Triangle's vertex as a SFCGAL::Point
 * @param triangle the input geometry
 * @param i index of the vertex in the triangle
 * @pre triangle must be a SFCGAL::Triangle
 * @pre i >= 0 and i < 3
 * @post returns a pointer to one of the vertices as a SFCGAL::Point. This
 * pointer is not writable and must not be deallocated by the caller
 * @return One the SFCGAL::Triangle's vertex as a SFCGAL::Point
 * @ingroup capi
 */
SFCGAL_API const sfcgal_geometry_t *
sfcgal_triangle_vertex(const sfcgal_geometry_t *triangle, int i);

/**
 * Sets one vertex of a SFCGAL::Triangle
 * @param triangle the input geometry
 * @param i index of the vertex in the triangle
 * @param vertex new point
 * @pre triangle must be a SFCGAL::Triangle
 * @pre vertex must be a SFCGAL::Point
 * @pre i >= 0 and i < 3
 * @post the ownership of the vertex is not taken. The caller is still
 * responsible of its deallocation.
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_triangle_set_vertex(sfcgal_geometry_t *triangle, int i,
                           const sfcgal_geometry_t *vertex);

/**
 * Sets one vertex of a SFCGAL::Triangle from two coordinates
 * @param triangle the input geometry
 * @param i index of the vertex in the triangle
 * @param x new point x
 * @param y new point y
 * @pre triangle must be a SFCGAL::Triangle
 * @pre i >= 0 and i < 3
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_triangle_set_vertex_from_xy(sfcgal_geometry_t *triangle, int i, double x,
                                   double y);

/**
 * Sets one vertex of a SFCGAL::Triangle from three coordinates
 * @param triangle the input geometry
 * @param i index of the vertex in the triangle
 * @param x new point x
 * @param y new point y
 * @param z new point z
 * @pre triangle must be a SFCGAL::Triangle
 * @pre i >= 0 and i < 3
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_triangle_set_vertex_from_xyz(sfcgal_geometry_t *triangle, int i,
                                    double x, double y, double z);

/**
 * Creates an empty SFCGAL::Polygon
 * @return an empty SFCGAL::Polygon
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_polygon_create();

/**
 * Creates an empty SFCGAL::Polygon from an exterior ring
 * @param ring the input geometry
 * @pre ring must be a SFCGAL::LineString
 * @post the ownership of the given ring is taken. The caller is not responsible
 * anymore of its deallocation
 * @return a SFCGAL::Polygon with ring as exterior ring
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_polygon_create_from_exterior_ring(sfcgal_geometry_t *ring);

/**
 * Returns the exterior ring of a given SFCGAL::Polygon
 * @param polygon the input geometry
 * @pre polygon must be a SFCGAL::Polygon
 * @pre polygon must not be empty
 * @post the returned ring is a SFCGAL::LineString, is not writable and must not
 * be deallocated by the caller
 * @return The exterior ing of polygon
 * @ingroup capi
 */
SFCGAL_API const sfcgal_geometry_t *
sfcgal_polygon_exterior_ring(const sfcgal_geometry_t *polygon);

/**
 * Returns the number of interior rings of a given SFCGAL::Polygon
 * @param polygon the input geometry
 * @pre polygon must be a SFCGAL::Polygon
 * @return The number of interior rings of polygon
 * @ingroup capi
 */
SFCGAL_API size_t
sfcgal_polygon_num_interior_rings(const sfcgal_geometry_t *polygon);

/**
 * Returns the ith interior ring of a given SFCGAL::Polygon
 * @param polygon the input geometry
 * @param i index of the ring in the polygon
 * @pre polygon must be a SFCGAL::Polygon
 * @pre i >= 0 and i < sfcgal_polygon_num_interior_rings
 * @post the returned ring is a SFCGAL::LineString, is not writable and must not
 * be deallocated by the caller
 * @return The ith interior ring of polygon
 * @ingroup capi
 */
SFCGAL_API const sfcgal_geometry_t *
sfcgal_polygon_interior_ring_n(const sfcgal_geometry_t *polygon, size_t i);

/**
 * Adds an interior ring to a given SFCGAL::Polygon
 * @param polygon the input geometry
 * @param ring the new ring
 * @pre polygon must be a SFCGAL::Polygon
 * @pre ring must be a SFCGAL::LineString
 * @post the ownership of the given ring is taken. The caller is not responsible
 * anymore of its deallocation
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_polygon_add_interior_ring(sfcgal_geometry_t *polygon,
                                 sfcgal_geometry_t *ring);

/**
 * Creates an empty SFCGAL::GeometryCollection
 * @return An empty SFCGAL::GeometryCollection
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_collection_create();

/**
 * Returns the number of geometries of a given SFCGAL::GeometryCollection
 * @param collection the input geometry
 * @pre collection is a SFCGAL::GeometryCollection
 * @return the number of geometries of a given SFCGAL::GeometryCollection
 * @deprecated Use sfcgal_geometry_num_geometries() instead
 * @ingroup capi
 */
SFCGAL_DEPRECATED("Use sfcgal_geometry_num_geometries instead.")
SFCGAL_API size_t
sfcgal_geometry_collection_num_geometries(const sfcgal_geometry_t *collection);

/**
 * Returns the ith geometry of a SFCGAL::GeometryCollection
 * @param collection the input geometry
 * @param i index of the geometry in the collection
 * @pre collection is a SFCGAL::GeometryCollection
 * @pre i >= 0 and i < sfcgal_geometry_collection_num_geometries
 * @post the returned Geometry is not writable and must not be deallocated by
 * the caller
 * @return the ith geometry of the collection
 * @ingroup capi
 */
SFCGAL_API const sfcgal_geometry_t *
sfcgal_geometry_collection_geometry_n(const sfcgal_geometry_t *collection,
                                      size_t                   i);

/**
 * Set the ith geometry of a given SFCGAL::GeometryCollection
 * @param collection the input geometry
 * @param geometry the new geometry
 * @param i index of the geometry in the collection
 * @pre collection is a SFCGAL::GeometryCollection
 * @pre i >= 0 and i < sfcgal_geometry_num_geometries( collection )
 * @post The ownership of the geometry is taken. The caller is not responsible
 * anymore of its deallocation.
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_geometry_collection_set_geometry_n(sfcgal_geometry_t *collection,
                                          sfcgal_geometry_t *geometry,
                                          size_t             i);

/**
 * Adds a Geometry to a given SFCGAL::GeometryCollection
 * @param collection the input geometry
 * @param geometry the new geometry
 * @pre collection must be a SFCGAL::GeometryCollection
 * @post the ownership of the given geometry is taken. The caller is not
 * responsible anymore of its deallocation
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_geometry_collection_add_geometry(sfcgal_geometry_t *collection,
                                        sfcgal_geometry_t *geometry);

/**
 * Creates an empty SFCGAL::MultiPoint
 * @return An empty SFCGAL::MultiPoint
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_multi_point_create();

/**
 * Creates an empty SFCGAL::MultiLineString
 * @return An empty SFCGAL::MultiLineString
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_multi_linestring_create();

/**
 * Creates an empty SFCGAL::MultiPolygon
 * @return An empty SFCGAL::MultiPolygon
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_multi_polygon_create();

/**
 * Creates an empty SFCGAL::MultiSolid
 * @return An empty SFCGAL::MultiSolid
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_multi_solid_create();

/**
 * Creates an empty SFCGAL::PolyhedralSurface
 * @return An empty SFCGAL::PolyhedralSurface
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_polyhedral_surface_create();

/**
 * Returns the number of patches of a given SFCGAL::PolyhedralSurface
 * @param polyhedral the input geometry
 * @pre polyhedral must be a SFCGAL::PolyhedralSurface
 * @return The number of patches of polyhedral
 * @ingroup capi
 */
SFCGAL_API size_t
sfcgal_polyhedral_surface_num_patches(const sfcgal_geometry_t *polyhedral);

/**
 * Returns the number of polygons of a given SFCGAL::PolyhedralSurface
 * @param polyhedral the input geometry
 * @pre polyhedral must be a SFCGAL::PolyhedralSurface
 * @return The number of patches of polyhedral
 * @deprecated Use sfcgal_polyhedral_surface_num_patches() instead.
 * @ingroup capi
 */
SFCGAL_DEPRECATED("Use sfcgal_polyhedral_surface_num_patches instead.")
SFCGAL_API size_t
sfcgal_polyhedral_surface_num_polygons(const sfcgal_geometry_t *polyhedral);

/**
 * Returns the ith patch of a given SFCGAL::PolyhedralSurface
 * @param polyhedral the input geometry
 * @param i index of the patch in the polyhedral
 * @pre polyhedral must be a SFCGAL::PolyhedralSurface
 * @pre i >= 0 and i < sfcgal_polyhedral_surface_num_patches(polyhedral)
 * @post the returned SFCGAL::Polygon is not writable and must not be
 * deallocated by the caller
 * @return The ith patch of geom
 * @ingroup capi
 */
SFCGAL_API const sfcgal_geometry_t *
sfcgal_polyhedral_surface_patch_n(const sfcgal_geometry_t *polyhedral,
                                  size_t                   i);

/**
 * Returns the ith polygon of a given SFCGAL::PolyhedralSurface
 * @param polyhedral the input geometry
 * @param i index of the patch in the polyhedral
 * @pre polyhedral must be a SFCGAL::PolyhedralSurface
 * @pre i >= 0 and i < sfcgal_polyhedral_surface_num_patches(polyhedral)
 * @post the returned SFCGAL::Polygon is not writable and must not be
 * deallocated by the caller
 * @return The ith patch of geom
 * @deprecated Use sfcgal_polyhedral_surface_patch_n() instead
 * @ingroup capi
 */
SFCGAL_DEPRECATED("Use sfcgal_polyhedral_surface_patch_n instead.")
SFCGAL_API const sfcgal_geometry_t *
sfcgal_polyhedral_surface_polygon_n(const sfcgal_geometry_t *polyhedral,
                                    size_t                   i);

/**
 * Adds a patch to a given SFCGAL::PolyhedralSurface
 * @param polyhedral the input geometry
 * @param patch the new patch
 * @pre polyhedral must be a SFCGAL::PolyhedralSurface
 * @pre patch must be a SFCGAL::Polygon
 * @post the ownership of the SFCGAL::Polygon is taken. The caller is not
 * responsible anymore of its deallocation
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_polyhedral_surface_add_patch(sfcgal_geometry_t *polyhedral,
                                    sfcgal_geometry_t *patch);

/**
 * Set the ith patch of a given SFCGAL::PolyhedralSurface
 * @param polyhedral the input geometry
 * @param patch the new patch
 * @param i index of the patch in the polyhedral
 * @pre polyhedral must be a SFCGAL::PolyhedralSurface.
 * @pre patch must be a SFCGAL::Polygon.
 * @pre i >= 0 and i < sfcgal_polyhedral_surface_num_patches(polyhedral)
 * @post The ownership of the polygon is taken. The caller is not responsible
 * anymore of its deallocation.
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_polyhedral_surface_set_patch_n(sfcgal_geometry_t *polyhedral,
                                      sfcgal_geometry_t *patch, size_t i);

/**
 * Adds a SFCGAL::Polygon to a given SFCGAL::PolyhedralSurface
 * @param polyhedral the input geometry
 * @param polygon the new polygon
 * @pre polyhedral must be a SFCGAL::PolyhedralSurface
 * @pre polygon must be a SFCGAL::Polygon
 * @post the ownership of the SFCGAL::Polygon is taken. The caller is not
 * responsible anymore of its deallocation
 * @deprecated Use sfcgal_polyhedral_surface_add_patch() instead.
 * @ingroup capi
 */
SFCGAL_DEPRECATED("Use sfcgal_polyhedral_surface_add_patch instead.")
SFCGAL_API void
sfcgal_polyhedral_surface_add_polygon(sfcgal_geometry_t *polyhedral,
                                      sfcgal_geometry_t *polygon);

/**
 * Creates an empty SFCGAL::TriangulatedSurface
 * @return An empty SFCGAL::TriangulatedSurface
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_triangulated_surface_create();

/**
 * Returns the number of patches of a given SFCGAL::TriangulatedSurface
 * @param tin the input geometry
 * @pre tin must be a SFCGAL::TriangulatedSurface
 * @return the number of patches of tin
 * @ingroup capi
 */
SFCGAL_API size_t
sfcgal_triangulated_surface_num_patches(const sfcgal_geometry_t *tin);

/**
 * Returns the number of triangles of a given SFCGAL::TriangulatedSurface
 * @param tin the input geometry
 * @pre tin must be a SFCGAL::TriangulatedSurface
 * @return The number of patches of tin
 * @deprecated Use sfcgal_triangulated_surface_num_patches() instead.
 * @ingroup capi
 */
SFCGAL_DEPRECATED("Use sfcgal_triangulated_surface_num_patches instead.")
SFCGAL_API size_t
sfcgal_triangulated_surface_num_triangles(const sfcgal_geometry_t *tin);

/**
 * Returns the ith patch of a given
 * SFCGAL::TriangulatedSurface
 * @param tin the input geometry
 * @param i index of the patch in the tin
 * @pre tin must be a SFCGAL::TriangulatedSurface
 * @pre i >= 0 and i < sfcgal_triangulated_surface_num_patches( tin )
 * @post the returned SFCGAL::Triangle is not writable and must not be
 * deallocated by the caller
 * @return The ith patch of tin
 * @ingroup capi
 */
SFCGAL_API const sfcgal_geometry_t *
sfcgal_triangulated_surface_patch_n(const sfcgal_geometry_t *tin, size_t i);

/**
 * Returns the ith SFCGAL::Triangle of a given SFCGAL::TriangulatedSurface
 * @param tin the input geometry
 * @param i index of the patch in the tin
 * @pre tin must be a SFCGAL::TriangulatedSurface
 * @pre i >= 0 and i < sfcgal_triangulated_surface_num_patches( tin )
 * @post the returned SFCGAL::Triangle is not writable and must not be
 * deallocated by the caller
 * @return The ith patch of tin
 * @deprecated Use sfcgal_triangulated_surface_patch_n() instead.
 * @ingroup capi
 */
SFCGAL_DEPRECATED("Use sfcgal_triangulated_surface_patch_n instead.")
SFCGAL_API const sfcgal_geometry_t *
sfcgal_triangulated_surface_triangle_n(const sfcgal_geometry_t *tin, size_t i);

/**
 * Set the ith patch of a given SFCGAL::TriangulatedSurface
 * @param tin the input geometry
 * @param patch the new patch
 * @param i index of the patch in the tin
 * @pre tin must be a SFCGAL::TriangulatedSurface
 * @pre patch must be a SFCGAL::Triangle.
 * @pre i >= 0 and i < sfcgal_triangulated_surface_num_patches( tin )
 * @post The ownership of the triangle is taken. The caller is not responsible
 * anymore of its deallocation.
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_triangulated_surface_set_patch_n(sfcgal_geometry_t *tin,
                                        sfcgal_geometry_t *patch, size_t i);

/**
 * Adds a patch to a given SFCGAL::TriangulatedSurface
 * @param tin the input geometry
 * @param patch the new patch
 * @pre tin must be a SFCGAL::TriangulatedSurface
 * @pre patch must be a SFCGAL::Triangle
 * @post the ownership of the SFCGAL::Triangle is taken. The caller is not
 * responsible anymore of its deallocation
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_triangulated_surface_add_patch(sfcgal_geometry_t *tin,
                                      sfcgal_geometry_t *patch);

/**
 * Adds a SFCGAL::Triangle to a given SFCGAL::TriangulatedSurface
 * @param tin the input geometry
 * @param triangle the new triangle
 * @pre tin must be a SFCGAL::TriangulatedSurface
 * @pre triangle must be a SFCGAL::Triangle
 * @post the ownership of the SFCGAL::Triangle is taken. The caller is not
 * responsible anymore of its deallocation
 * @deprecated Use sfcgal_triangulated_surface_add_patch() instead.
 * @ingroup capi
 */
SFCGAL_DEPRECATED("Use sfcgal_triangulated_surface_add_patch instead.")
SFCGAL_API void
sfcgal_triangulated_surface_add_triangle(sfcgal_geometry_t *tin,
                                         sfcgal_geometry_t *triangle);

/**
 * Creates an empty SFCGAL::Solid
 * @return A SFCGAL::Solid
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_solid_create();

/**
 * Creates a SFCGAL::Solid from an exterior shell
 * @param shell the input geometry
 * @pre shell must be a SFCGAL::PolyhedralSurface
 * @post the ownership of the given shell is taken. The caller is not
 * responsible anymore of its deallocation
 * @return A SFCGAL::Solid with shell as exterior shell
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_solid_create_from_exterior_shell(sfcgal_geometry_t *shell);

/**
 * Returns the number of shells of a given SFCGAL::Solid
 * @param solid the input geometry
 * @pre solid must be a SFCGAL::Solid
 * @return The number of shells of solid
 * @ingroup capi
 */
SFCGAL_API size_t
sfcgal_solid_num_shells(const sfcgal_geometry_t *solid);

/**
 * Returns the ith shell of a given SFCGAL::Solid
 * @param solid the input geometry
 * @param i index of the shell in the solid
 * @pre solid must be a SFCGAL::Solid
 * @pre i >= 0 and i < sfcgal_solid_num_shells( tin )
 * @post the returned SFCGAL::PolyhedralSurface is not writable and must not be
 * deallocated by the caller
 * @return The ith shell of solid
 * @ingroup capi
 */
SFCGAL_API const sfcgal_geometry_t *
sfcgal_solid_shell_n(const sfcgal_geometry_t *solid, size_t i);

/**
 * Adds a shell to a given SFCGAL::Solid
 * @param solid the input geometry
 * @param shell the new geometry
 * @pre solid must be a SFCGAL::Solid
 * @pre shell must be a SFCGAL::PolyhedralSurface
 * @post the ownership of the shell is taken. The caller is not responsible
 * anymore of its deallocation
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_solid_add_interior_shell(sfcgal_geometry_t *solid,
                                sfcgal_geometry_t *shell);

/**
 * Set the exterior shell of a given SFCGAL::Solid
 * @param solid the input geometry
 * @param shell the new geometry
 * @pre solid must be a SFCGAL::Solid
 * @pre shell must be a SFCGAL::PolyhedralSurface
 * @post the ownership of the shell is taken. The caller is not responsible
 * anymore of its deallocation
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_solid_set_exterior_shell(sfcgal_geometry_t *solid,
                                sfcgal_geometry_t *shell);

/**
 * Gets the validity flag of the geometry.
 * @param geom the input geometry
 * @return the validity flag of geom
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_has_validity_flag(const sfcgal_geometry_t *geom);

/**
 * Sets the validity flag of the geometry.
 * FIXME We better have geometry constructors to directly build valid geometries
 * @param geom the input geometry
 * @param valid 1 to force validity or 0 to force invalidity
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_geometry_force_valid(sfcgal_geometry_t *geom, int valid);

/*--------------------------------------------------------------------------------------*
 *
 * Support for SFCGAL::PreparedGeometry
 *
 *--------------------------------------------------------------------------------------*/

/**
 * Opaque type that represents the C++ type SFCGAL::PreparedGeometry
 * @ingroup capi
 */
typedef void sfcgal_prepared_geometry_t;

/**
 * Creates an empty SFCGAL::PreparedGeometry
 * @return an empty SFCGAL::PreparedGeometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_prepared_geometry_t *
sfcgal_prepared_geometry_create();

/**
 * Creates a SFCGAL::PreparedGeometry from a Geometry and an SRID
 * @param geometry the input geometry
 * @param srid the srid of the geometry
 * @return a SFCGAL::PreparedGeometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_prepared_geometry_t *
sfcgal_prepared_geometry_create_from_geometry(sfcgal_geometry_t *geometry,
                                              srid_t             srid);

/**
 * Deletes a given SFCGAL::PreparedGeometry
 * @param prepared the input geometry
 * @pre prepared must be a SFCGAL::PreparedGeometry
 * @post the underlying Geometry linked to the given PreparedGeometry is also
 * deleted
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_prepared_geometry_delete(sfcgal_prepared_geometry_t *prepared);

/**
 * Returns the SFCGAL::Geometry associated with a given SFCGAL::PreparedGeometry
 * @param prepared the input geometry
 * @pre prepared must be a SFCGAL::PreparedGeometry
 * @post the returned Geometry is not writable and must not be deallocated by
 * the caller
 * @return the SFCGAL::Geometry of prepared
 * @ingroup capi
 */
SFCGAL_API const sfcgal_geometry_t *
sfcgal_prepared_geometry_geometry(const sfcgal_prepared_geometry_t *prepared);

/**
 * Sets the Geometry associated with the given SFCGAL::PreparedGeometry
 * @param prepared the input geometry
 * @param geometry the geometry to set
 * @pre prepared must be a SFCGAL::PreparedGeometry
 * @post the ownership of the given geometry is taken. The caller is not
 * responsible anymore of its deallocation
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_prepared_geometry_set_geometry(sfcgal_prepared_geometry_t *prepared,
                                      sfcgal_geometry_t          *geometry);

/**
 * Returns SRID associated with a given SFCGAL::PreparedGeometry
 * @param prepared the input geometry
 * @pre prepared must be a SFCGAL::PreparedGeometry
 * @return the srid of prepared
 * @ingroup capi
 */
SFCGAL_API srid_t
sfcgal_prepared_geometry_srid(const sfcgal_prepared_geometry_t *prepared);

/**
 * Sets SRID associated with a given SFCGAL::PreparedGeometry
 * @param prepared the input geometry
 * @param srid the srid of the geometry
 * @pre prepared must be a SFCGAL::PreparedGeometry
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_prepared_geometry_set_srid(sfcgal_prepared_geometry_t *prepared,
                                  srid_t                      srid);

/**
 * Returns an EWKT representation of the given SFCGAL::PreparedGeometry
 * @param[in] prepared the input geometry
 * @param[in] num_decimals number of decimals. -2 for a variable number of
 *            decimals.
 * @param[out] buffer the output buffer
 * @param[out] len the size of the @p buffer
 * @post @p buffer is returned allocated and must be freed by the caller with
 *       sfcgal_free_buffer()
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_prepared_geometry_as_ewkt(const sfcgal_prepared_geometry_t *prepared,
                                 int num_decimals, char **buffer, size_t *len);

/*--------------------------------------------------------------------------------------*
 *
 * I/O functions
 *
 *--------------------------------------------------------------------------------------*/

/**
 * io::readWKT
 */

/**
 * Parse a Well-Known Text (WKT) representation into a Geometry.
 * @param str The Well-Known Text (WKT) string representing the geometry.
 * @param len The size of @p str
 * @post The returned geometry must be deallocated by the caller with
 * sfcgal_geometry_delete()
 * @return A SFCGAL::Geometry from the WKT
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_io_read_wkt(const char *str, size_t len);

/**
 * Parse an Extended Well-Known Text (EWKT) representation into a
 * PreparedGeometry.
 * @param str The Well-Known Text (WKT) string representing the geometry.
 * @param len The size of @p str
 * @post The returned prepared geometry must be deallocated by the caller with
 * @pre sfcgal_prepared_geometry_delete
 * @return A SFCGAL::Geometry from the EWKT
 * @ingroup capi
 */
SFCGAL_API sfcgal_prepared_geometry_t *
sfcgal_io_read_ewkt(const char *str, size_t len);

/**
 * io::readWKB
 */

/**
 * Parse a Well-Known Binary (WKB) representation into a Geometry.
 * @param str The Well-Known Text (WKB) string representing the geometry.
 * @param len The size of @p str
 * @post The returned geometry must be deallocated by the caller with
 * sfcgal_geometry_delete()
 * @return A SFCGAL::Geometry from the WKB
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_io_read_wkb(const char *str, size_t len);

/**
 * Serialization
 */

/**
 * Convert a Geometry to its binary representation.
 * @param[in] geom the prepared geometry
 * @param[out] buffer The binary buffer
 * @param[out] len The size of @p buffer
 * @post @p buffer is returned allocated and must be freed by the caller with
 *       sfcgal_free_buffer()
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_io_write_binary_prepared(const sfcgal_prepared_geometry_t *geom,
                                char **buffer, size_t *len);

/**
 * Read a SFCGAL::PreparedGeometry from a binary representation.
 * @param str The string representing the geometry.
 * @param len The size of @p str
 * @post The returned prepared geometry must be deallocated by the caller with
 * sfcgal_prepared_geometry_delete()
 * @return a SFCGAL::PreparedGeometry from its binary representation
 * @ingroup capi
 */
SFCGAL_API sfcgal_prepared_geometry_t *
sfcgal_io_read_binary_prepared(const char *str, size_t len);

/*--------------------------------------------------------------------------------------*
 *
 * Spatial processing
 *
 *--------------------------------------------------------------------------------------*/

/**
 * Tests the intersection of geom1 and geom2
 * @param geom1 the first input geometry
 * @param geom2 the second input geometry
 * @pre isValid(geom1) == true
 * @pre isValid(geom2) == true
 * @return 1 if geometries intersect in 2D. 0 otherwise.
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_intersects(const sfcgal_geometry_t *geom1,
                           const sfcgal_geometry_t *geom2);

/**
 * Tests the 3D intersection of geom1 and geom2
 * @param geom1 the first input geometry
 * @param geom2 the second input geometry
 * @pre isValid(geom1) == true
 * @pre isValid(geom2) == true
 * @return 1 if geometries intersect in 3D. 0 otherwise.
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_intersects_3d(const sfcgal_geometry_t *geom1,
                              const sfcgal_geometry_t *geom2);

/**
 * Returns the intersection of geom1 and geom2
 * @param geom1 the first input geometry
 * @param geom2 the second input geometry
 * @pre isValid(geom1) == true
 * @pre isValid(geom2) == true
 * @post isValid(return) == true
 * @return The intersection of both geometries
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_intersection(const sfcgal_geometry_t *geom1,
                             const sfcgal_geometry_t *geom2);

/**
 * Returns the 3D intersection of geom1 and geom2
 * @param geom1 the first input geometry
 * @param geom2 the second input geometry
 * @pre isValid(geom1) == true
 * @pre isValid(geom2) == true
 * @post isValid(return) == true
 * @return The 3D intersection of both geometries
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_intersection_3d(const sfcgal_geometry_t *geom1,
                                const sfcgal_geometry_t *geom2);

/**
 * Returns the difference of geom1 and geom2
 * @param geom1 the first input geometry
 * @param geom2 the second input geometry
 * @pre isValid(geom1) == true
 * @pre isValid(geom2) == true
 * @post isValid(return) == true
 * @return The difference of both geometries
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_difference(const sfcgal_geometry_t *geom1,
                           const sfcgal_geometry_t *geom2);

/**
 * Returns the 3D difference of geom1 and geom2
 * @param geom1 the first input geometry
 * @param geom2 the second input geometry
 * @pre isValid(geom1) == true
 * @pre isValid(geom2) == true
 * @post isValid(return) == true
 * @return The 3D difference of both geometries
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_difference_3d(const sfcgal_geometry_t *geom1,
                              const sfcgal_geometry_t *geom2);

/**
 * Returns the union of geom1 and geom2
 * @param geom1 the first input geometry
 * @param geom2 the second input geometry
 * @pre isValid(geom1) == true
 * @pre isValid(geom2) == true
 * @post isValid(return) == true
 * @return The union of both geometries
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_union(const sfcgal_geometry_t *geom1,
                      const sfcgal_geometry_t *geom2);

/**
 * Returns the 3D union of geom1 and geom2
 * @param geom1 the first input geometry
 * @param geom2 the second input geometry
 * @pre isValid(geom1) == true
 * @pre isValid(geom2) == true
 * @post isValid(return) == true
 * @return The 3D union of both geometries
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_union_3d(const sfcgal_geometry_t *geom1,
                         const sfcgal_geometry_t *geom2);

/**
 * Returns the convex hull of geom
 * @param geom the input geometry
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return The 2D convex hull of geom
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_convexhull(const sfcgal_geometry_t *geom);

/**
 * Returns the 3D convex hull of geom
 * @param geom the input geometry
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return The 3D convex hull of geom
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_convexhull_3d(const sfcgal_geometry_t *geom);

/**
 * Returns the volume of geom (must be a volume)
 * @param geom the input geometry
 * @pre isValid(geom) == true
 * @return The volume of geom
 * @ingroup capi
 */
SFCGAL_API double
sfcgal_geometry_volume(const sfcgal_geometry_t *geom);

/**
 * Returns the area of geom
 * @param geom the input geometry
 * @pre isValid(geom) == true
 * @return The area of geom
 * @ingroup capi
 */
SFCGAL_API double
sfcgal_geometry_area(const sfcgal_geometry_t *geom);

/**
 * Returns the 3D area of geom
 * @param geom the input geometry
 * @pre isValid(geom) == true
 * @return The 3D area of geom
 * @ingroup capi
 */
SFCGAL_API double
sfcgal_geometry_area_3d(const sfcgal_geometry_t *geom);

/**
 * Tests if the given Geometry is planar
 * @param geom the input geometry
 * @pre isValid(geom) == true
 * @return 1 if the geometry is planar. 0 otherwise.
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_is_planar(const sfcgal_geometry_t *geom);

/**
 * Returns the orientation of the given SFCGAL::Polygon
 * -1 for a counter clockwise orientation
 * 1 for a clockwise orientation
 * 0 for an invalid or undetermined orientation
 * @param geom the input geometry
 * @pre geom is a SFCGAL::Polygon
 * @pre isValid(geom) == true
 * @return The orientation of the geometry as an integer
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_orientation(const sfcgal_geometry_t *geom);

/**
 * Returns a tesselation of the given Geometry
 * @param geom the input geometry
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return The tesselation of geom
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_tesselate(const sfcgal_geometry_t *geom);

/**
 * Returns a triangulation of the given Geometry
 * @param geom the input geometry
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return The triangulation of geom
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_triangulate_2dz(const sfcgal_geometry_t *geom);

/**
 * Returns an extrusion of the given Geometry
 * @param geom the input geometry
 * @param ex extrude value on x axis
 * @param ey extrude value on y axis
 * @param ez extrude value on z axis
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return The extruded geom
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_extrude(const sfcgal_geometry_t *geom, double ex, double ey,
                        double ez);

/**
 * Convert a SFCGAL::PolyhedralSurface to a SFCGAL::Solid
 * @param geom the input geometry
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return The converted geom
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_make_solid(const sfcgal_geometry_t *geom);

/**
 * Force a Left Handed Rule on the given Geometry
 * @param geom the input geometry
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return The LHR geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_force_lhr(const sfcgal_geometry_t *geom);

/**
 * Force a Right Handed Rule on the given Geometry
 * @param geom the input geometry
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return The RHR geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_force_rhr(const sfcgal_geometry_t *geom);

/**
 * Computes the distance of the two given Geometry objects
 * @param geom1 the first input geometry
 * @param geom2 the second input geometry
 * @pre isValid(geom1) == true
 * @pre isValid(geom2) == true
 * @return The distance of the two given SFCGAL::Geometry
 * @ingroup capi
 */
SFCGAL_API double
sfcgal_geometry_distance(const sfcgal_geometry_t *geom1,
                         const sfcgal_geometry_t *geom2);

/**
 * Computes the 3D distance of the two given Geometry objects
 * @param geom1 the first input geometry
 * @param geom2 the second input geometry
 * @pre isValid(geom1) == true
 * @pre isValid(geom2) == true
 * @return The 3D distance of the two given SFCGAL::Geometry
 * @ingroup capi
 */
SFCGAL_API double
sfcgal_geometry_distance_3d(const sfcgal_geometry_t *geom1,
                            const sfcgal_geometry_t *geom2);

/**
 * Round coordinates of the given Geometry
 * @param geom the input geometry
 * @param r rounding precision
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return The rounded geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_round(const sfcgal_geometry_t *geom, int r);

/**
 * Returns the minkowski sum geom1 + geom2
 * @param geom1 the first input geometry
 * @param geom2 the second input geometry
 * @pre isValid(geom1) == true
 * @pre isValid(geom2) == true
 * @post isValid(return) == true
 * @return The minkowski sum of geom1 and geom2
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_minkowski_sum(const sfcgal_geometry_t *geom1,
                              const sfcgal_geometry_t *geom2);

/**
 * Returns the offset polygon of the given Geometry.
 * @param geom the input geometry
 * @param radius offset distance
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return The offset polygon of geom
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_offset_polygon(const sfcgal_geometry_t *geom, double radius);

/**
 * Returns the straight skeleton of the given Geometry
 * @param geom the input geometry
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return The straight skeleton of the given Geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_straight_skeleton(const sfcgal_geometry_t *geom);

/**
 * Returns the straight skeleton of the given Geometry with the distance to the
 * border as M coordinate
 * @param geom the input geometry
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return The straight skeleton with the distance to the border as M coordinate
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_straight_skeleton_distance_in_m(const sfcgal_geometry_t *geom);

/**
 * Returns the extrude straight skeleton of the given SFCGAL::Polygon
 * @param geom the input geometry
 * @param height extrusion height
 * @pre geom must be a SFCGAL::Polygon
 * @pre isValid(geom) == true
 * @pre height != 0
 * @return The extrude straight skeleton of geom
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_extrude_straight_skeleton(const sfcgal_geometry_t *geom,
                                          double                   height);

/**
 * Returns the union of the polygon z-extrusion (with respect to
 * building_height) and the extrude straight skeleton (with respect to
 * roof_height) of the given SFCGAL::Polygon
 * @param geom the input geometry
 * @param building_height extrusion height of walls
 * @param roof_height extrusion height of roof
 * @pre geom must be a SFCGAL::Polygon
 * @pre isValid(geom) == true
 * @pre roof_height != 0
 * @return The the union of the polygon z-extrusion
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_extrude_polygon_straight_skeleton(const sfcgal_geometry_t *geom,
                                                  double building_height,
                                                  double roof_height);

/**
 * Returns the approximate medial axis for the given SFCGAL::Polygon
 * Approximate medial axis is based on straight skeleton
 * @param geom the input geometry
 * @pre isValid(geom) == true
 * @return A SFCGAL::MultiLineString representing the approximate medial axis of
 * geom
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_approximate_medial_axis(const sfcgal_geometry_t *geom);

/**
 * Returns the straight skeleton partition for the given SFCGAL::Polygon
 * @param geom the input geometry
 * @param autoOrientation if 1 try to find the best orientation
 * @pre isValid(geom) == true
 * @pre geom must be a SFCGAL::Polygon, SFCGAL::Triangle or MultiPolygon
 * @post isValid(return) == true
 * @return A geometry representing the straight skeleton partition of geom
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_straight_skeleton_partition(const sfcgal_geometry_t *geom,
                                            bool autoOrientation);

/**
 * Tests the 2D coverage of geom1 and geom2
 * @param geom1 the first input geometry
 * @param geom2 the second input geometry
 * @pre isValid(geom1) == true
 * @pre isValid(geom2) == true
 * @return 1 if geom1 or geom2 2D covers the other (including if both cover each
 * other), 0 otherwise.
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_covers(const sfcgal_geometry_t *geom1,
                       const sfcgal_geometry_t *geom2);

/**
 * Tests the 3D coverage of geom1 and geom2
 * @param geom1 the first input geometry
 * @param geom2 the second input geometry
 * @pre isValid(geom1) == true
 * @pre isValid(geom2) == true
 * @return 1 if geom1 or geom2 3D covers the other (including if both cover each
 * other), 0 otherwise.
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_covers_3d(const sfcgal_geometry_t *geom1,
                          const sfcgal_geometry_t *geom2);

/**
 * Returns the substring of the given SFCGAL::LineString between fractional
 * distances
 * @param geom the input geometry
 * @param start distance from start to begin extraction
 * @param end distance to end to end extraction
 * @pre isValid(geom) == true
 * @pre geom is a Linestring
 * @pre -1 <= start <= 1
 * @pre -1 <= end <= 1
 * @post isValid(return) == true
 * @return A new SFCGAL::LineString representing the portion of the input line
 * between startFrac and endFrac. If startFrac >= endFrac, the returned
 * SFCGAL::LineString may be empty.
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_line_sub_string(const sfcgal_geometry_t *geom, double start,
                                double end);

#if !_MSC_VER
/**
 * Returns the alpha shapes of geom
 * @pre isValid(geom) == true
 * @pre alpha >= 0
 * @post isValid(return) == true
 * @param geom the input geometry
 * @param alpha It can have values from 0 to infinity. Smaller alpha values
 * produce more concave results. Alpha values greater than some data-dependent
 * value produce the convex hull of the input.
 * @param allow_holes defines whether alpha shapes are allowed to contain holes
 * @return A geometry representing the alpha shapes of the geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_alpha_shapes(const sfcgal_geometry_t *geom, double alpha,
                             bool allow_holes);

/**
 * Returns the optimal alpha shapes of geom
 * @pre isValid(geom) == true
 * @pre alpha >= 0
 * @pre nb_components >= 0
 * @post isValid(return) == true
 * @param geom the input geometry
 * @param allow_holes defines whether alpha shapes are allowed to contain holes
 * @param nb_components the number of connected components in the output
 * geometry
 * @return A geometry representing the optimal alpha shapes of geom
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_optimal_alpha_shapes(const sfcgal_geometry_t *geom,
                                     bool allow_holes, size_t nb_components);
#endif

/**
 * Returns the 3D alpha wrapping of a geometry
 * @pre isValid(geom) == true
 * @pre relative_alpha >= 0
 * @pre relative_offset >= 0
 * @post isValid(return) == true
 * @param geom input geometry
 * @param relative_alpha This parameter is used to determine which features will
 *  appear in the output. A small relative_alpha will produce an output less
 *  complex but less faithful to the input.
 * @param relative_offset  This parameter controls the tightness of the result.
 *  A large relative_offset parameter will tend to better preserve sharp
 * features as projection. If this parameter is equal to 0, it is computed from
 * the alpha parameter
 * @return A SFCGAL::PolyhedralSurface representing the 3D alpha wrapping of the
 * geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_alpha_wrapping_3d(const sfcgal_geometry_t *geom,
                                  size_t                   relative_alpha,
                                  size_t                   relative_offset);

/**
 * Returns the envelope of geom
 * @param geom input geometry
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return The 2D envelope of geom
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_envelope(const sfcgal_geometry_t *geom);

/**
 * Returns the 3d envelope of geom
 * @param geom input geometry
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return The 3D envelope of geom
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_envelope_3d(const sfcgal_geometry_t *geom);

/**
 * Returns the 2D length of geom
 * @param geom input geometry
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return The 2D length of geom
 * @ingroup capi
 */
SFCGAL_API double
sfcgal_geometry_length(const sfcgal_geometry_t *geom);

/**
 * Returns the 3D length of geom
 * @param geom input geometry
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return The 3D length of geom
 * @ingroup capi
 */
SFCGAL_API double
sfcgal_geometry_length_3d(const sfcgal_geometry_t *geom);

/**
 * Returns the boundary of geom
 * @param geom input geometry
 * @pre isValid(geom) == true
 * @post the caller is responsible for the boundary deallocation.
 * @return The boundary of geom
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_boundary(const sfcgal_geometry_t *geom);

/**
 * Returns a SFCGAL::Point representing the geometry centroid
 * @param geom input geometry
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return A SFCGAL::Point representing the geometry centroid
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_centroid(const sfcgal_geometry_t *geom);

/**
 * Returns a SFCGAL::Point representing the geometry centroid
 * @param geom input geometry
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return A SFCGAL::Point representing the 3D geometry centroid
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_centroid_3d(const sfcgal_geometry_t *geom);

/**
 * Returns true if geom1 is equals to geom2.
 *
 * For each point of geom1 there is a point in geom2.
 * @param geom1 first input geometry
 * @param geom2 second input geometry
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return 1 if geometries are almost equal. 0 otherwise.
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_is_equals(const sfcgal_geometry_t *geom1,
                          const sfcgal_geometry_t *geom2);

/**
 * Returns true if geom1 is almost equals to geom2.
 *
 * For each point of geom1 there is a point in geom2 within tolerance distance.
 * @param geom1 the first geometry
 * @param geom2 the second geometry
 * @param tolerance the tolerance
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return 1 if geometries are almost equal. 0 otherwise.
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_geometry_is_almost_equals(const sfcgal_geometry_t *geom1,
                                 const sfcgal_geometry_t *geom2,
                                 double                   tolerance);

/**
 * Returns the y monotone partition of a geometry (polygon without hole)
 * @param geom input geometry
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return The y monotone partition of geom
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_y_monotone_partition_2(const sfcgal_geometry_t *geom);

/**
 * Returns the approximal convex partition of a geometry (polygon without hole)
 * @param geom input geometry
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return The approximal convex partition of geom
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_approx_convex_partition_2(const sfcgal_geometry_t *geom);

/**
 * Returns the greene approximal convex partition of a geometry (polygon without
 * hole)
 * @param geom input geometry
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return The greene approximal convex partition of geom
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_greene_approx_convex_partition_2(const sfcgal_geometry_t *geom);

/**
 * Returns the optimal convex partition of a geometry (polygon without hole)
 * @param geom input geometry
 * @pre isValid(geom) == true
 * @post isValid(return) == true
 * @return The optimal convex partition of geom
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_optimal_convex_partition_2(const sfcgal_geometry_t *geom);

/**
 * Returns the visibility polygon of a SFCGAL::Point inside a SFCGAL::Polygon
 * @param polygon input geometry
 * @param point input geometry
 * @pre polygon is a valid geometry
 * @pre point must be inside polygon or on the boundary
 * @return The visibility polygon of point inside polygon
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_visibility_point(const sfcgal_geometry_t *polygon,
                                 const sfcgal_geometry_t *point);

/**
 * @brief build the visibility polygon of the segment [pointA ; pointB] on a
 * SFCGAL::Polygon
 * @param polygon input geometry
 * @param pointA input geometry
 * @param pointB input geometry
 * @ingroup capi
 * @pre polygon is a valid geometry
 * @pre pointA and pointB must be vertices of poly, adjacents and respect the
 * @return A visibility polygon
 * direction
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_visibility_segment(const sfcgal_geometry_t *polygon,
                                   const sfcgal_geometry_t *pointA,
                                   const sfcgal_geometry_t *pointB);

/**
 * Buffer3D types
 * @ingroup capi
 */
typedef enum {
  SFCGAL_BUFFER3D_ROUND,
  SFCGAL_BUFFER3D_CYLSPHERE,
  SFCGAL_BUFFER3D_FLAT
} sfcgal_buffer3d_type_t;

/**
 * Computes a 3D buffer around a geometry
 * @param geom The input geometry (must be a SFCGAL::Point or
 * SFCGAL::LineString)
 * @param radius The buffer radius
 * @param segments The number of segments to use for approximating curved
 * surfaces
 * @param buffer_type The type of buffer to compute (ROUND, CYLSPHERE, or FLAT)
 * @return A new geometry representing the 3D buffer
 * @pre isValid(geom) == true
 * @pre radius > 0
 * @pre segments > 3
 * @post isValid(return) == true
 * @post The returned geometry must be deallocated by the caller with
 * sfcgal_geometry_delete()
 * @return A 3D buffer as a SFCGAL::PolyhedralSurface
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_buffer3d(const sfcgal_geometry_t *geom, double radius,
                         int segments, sfcgal_buffer3d_type_t buffer_type);

/*--------------------------------------------------------------------------------------*
 *
 * Transformation
 *
 *--------------------------------------------------------------------------------------*/

/**
 * Rotates a geometry around the origin (0,0,0) by a given angle
 * @param geom The geometry to rotate
 * @param angle Rotation angle in radians
 * @return The rotated geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_rotate(const sfcgal_geometry_t *geom, double angle);

/**
 * Rotates a geometry around a specified point by a given angle
 * @param geom The geometry to rotate
 * @param angle Rotation angle in radians
 * @param cx X-coordinate of the center point
 * @param cy Y-coordinate of the center point
 * @return The rotated geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_rotate_2d(const sfcgal_geometry_t *geom, double angle,
                          double cx, double cy);

/**
 * Rotates a 3D geometry around a specified axis by a given angle
 * @param geom The geometry to rotate
 * @param angle Rotation angle in radians
 * @param ax X-coordinate of the axis vector
 * @param ay Y-coordinate of the axis vector
 * @param az Z-coordinate of the axis vector
 * @return The rotated geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_rotate_3d(const sfcgal_geometry_t *geom, double angle,
                          double ax, double ay, double az);

/**
 * Rotates a 3D geometry around a specified axis and center point by a given
 * angle
 * @param geom The geometry to rotate
 * @param angle Rotation angle in radians
 * @param ax X-coordinate of the axis vector
 * @param ay Y-coordinate of the axis vector
 * @param az Z-coordinate of the axis vector
 * @param cx X-coordinate of the center point
 * @param cy Y-coordinate of the center point
 * @param cz Z-coordinate of the center point
 * @return The rotated geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_rotate_3d_around_center(const sfcgal_geometry_t *geom,
                                        double angle, double ax, double ay,
                                        double az, double cx, double cy,
                                        double cz);

/**
 * Rotates a geometry around the X axis by a given angle
 * @param geom The geometry to rotate
 * @param angle Rotation angle in radians
 * @return The rotated geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_rotate_x(const sfcgal_geometry_t *geom, double angle);

/**
 * Rotates a geometry around the Y axis by a given angle
 * @param geom The geometry to rotate
 * @param angle Rotation angle in radians
 * @return The rotated geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_rotate_y(const sfcgal_geometry_t *geom, double angle);

/**
 * Rotates a geometry around the Z axis by a given angle
 * @param geom The geometry to rotate
 * @param angle Rotation angle in radians
 * @return The rotated geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_rotate_z(const sfcgal_geometry_t *geom, double angle);

/**
 * Scale a geometry by a given factor
 * @param geom The geometry to scale
 * @param s Scale factor
 * @return The scaled geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_scale(const sfcgal_geometry_t *geom, double s);

/**
 * Scale a geometry by different factors for each dimension
 * @param geom The geometry to scale
 * @param sx Scale factor for x dimension
 * @param sy Scale factor for y dimension
 * @param sz Scale factor for z dimension
 * @return The scaled geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_scale_3d(const sfcgal_geometry_t *geom, double sx, double sy,
                         double sz);

/**
 * Scale a geometry by different factors for each dimension around a center
 * point
 * @param geom The geometry to scale
 * @param sx Scale factor for x dimension
 * @param sy Scale factor for y dimension
 * @param sz Scale factor for z dimension
 * @param cx X-coordinate of the center point
 * @param cy Y-coordinate of the center point
 * @param cz Z-coordinate of the center point
 * @return The scaled geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_scale_3d_around_center(const sfcgal_geometry_t *geom, double sx,
                                       double sy, double sz, double cx,
                                       double cy, double cz);

/**
 * Translate a geometry by a 3D vector
 * @param geom the geometry to translate
 * @param dx x component of the translation vector
 * @param dy y component of the translation vector
 * @param dz z component of the translation vector
 * @return The translated geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_translate_3d(const sfcgal_geometry_t *geom, double dx,
                             double dy, double dz);

/**
 * Translate a geometry by a 2D vector
 * @param geom the geometry to translate
 * @param dx x component of the translation vector
 * @param dy y component of the translation vector
 * @return The translated geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_translate_2d(const sfcgal_geometry_t *geom, double dx,
                             double dy);

/**
 * Simplify a geometry
 * @param geom the geometry to simplify
 * @param threshold threshold parameter
 * @param preserveTopology preserve the topology
 * @return The simplified geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_geometry_simplify(const sfcgal_geometry_t *geom, double threshold,
                         bool preserveTopology);

/*--------------------------------------------------------------------------------------*
 *
 * Support for SFCGAL::Primitive class hierarchy
 *
 *--------------------------------------------------------------------------------------*/

/**
 * sfcgal_primitive_t is an opaque pointer type that is used to represent a
 * pointer to SFCGAL::Primitive
 * @ingroup capi
 */
typedef void sfcgal_primitive_t;

/**
 * Primitive types
 * @ingroup capi
 */
typedef enum {
  SFCGAL_TYPE_CYLINDER = 0,
  SFCGAL_TYPE_SPHERE   = 1,
  SFCGAL_TYPE_TORUS    = 2,
  SFCGAL_TYPE_BOX      = 3,
  SFCGAL_TYPE_CUBE     = 4,
  SFCGAL_TYPE_CONE     = 5
} sfcgal_primitive_type_t;

/**
 * @brief Creates a new primitive of the specified type with default parameters.
 * @param primitive_type The type of primitive to create.
 * @return A pointer to an allocated ::sfcgal_primitive_t object. The
 * returned object must be freed by calling sfcgal_primitive_delete().
 * @post The returned pointer is non-NULL if allocation succeeds.
 * @return A SFCGAL::Primitive with default parameters
 * @ingroup capi
 */
SFCGAL_API sfcgal_primitive_t *
sfcgal_primitive_create(sfcgal_primitive_type_t primitive_type);

/**
 * @brief Deletes a previously allocated primitive.
 * @param primitive Pointer to the primitive to delete.
 * @pre @p primitive must be a valid pointer returned by a creation function
 * ( e.g. sfcgal_primitive_create() ).
 * @post After this call, the pointer is invalid and must not be used.
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_primitive_delete(sfcgal_primitive_t *primitive);

/**
 * @brief Retrieves the list of primitive parameters as JSON array.
 * @param primitive Pointer to the primitive.
 * @param[out] buffer the output buffer
 * @param[out] len the size of the @p buffer
 * @post @p buffer is returned allocated and must be freed by the caller with
 *       sfcgal_free_buffer()
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_primitive_parameters(sfcgal_primitive_t *primitive, char **buffer,
                            size_t *len);

/**
 * @brief Retrieves the value of a primitive parameter as a double.
 * @param primitive Pointer to the primitive.
 * @param name Name of the parameter to retrieve.
 * @return The parameter value as a double.
 * @pre The parameter identified by @p name must exist and be of type double.
 * @return Primitive parameter as a double
 * @ingroup capi
 */
SFCGAL_API double
sfcgal_primitive_parameter_double(const sfcgal_primitive_t *primitive,
                                  const char               *name);

/**
 * @brief Sets the value of a primitive parameter as a double.
 * @param primitive Pointer to the primitive.
 * @param name Name of the parameter to set.
 * @param parameter The new parameter value.
 * @pre The parameter identified by @p name must exist and accept values of type
 * double.
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_primitive_set_parameter_double(sfcgal_primitive_t *primitive,
                                      const char *name, double parameter);

/**
 * @brief Retrieves the value of a primitive parameter as an unsigned int.
 * @param primitive Pointer to the primitive.
 * @param name Name of the parameter to retrieve.
 * @return The parameter value as an unsigned int.
 * @pre The parameter identified by @p name must exist and be of type unsigned
 * int.
 * @return Primitive parameter as an unsigned int
 * @ingroup capi
 */
SFCGAL_API unsigned int
sfcgal_primitive_parameter_int(const sfcgal_primitive_t *primitive,
                               const char               *name);

/**
 * @brief Sets the value of a primitive parameter as an unsigned int.
 * @param primitive Pointer to the primitive.
 * @param name Name of the parameter to set.
 * @param parameter The new parameter value.
 * @pre The parameter identified by @p name must exist and accept values of type
 * unsigned int.
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_primitive_set_parameter_int(sfcgal_primitive_t *primitive,
                                   const char *name, unsigned int parameter);

/**
 * @brief Retrieves the value of a primitive parameter as a 3D point: an array
 * of three doubles (x, y, z).
 * @param primitive Pointer to the primitive.
 * @param name Name of the parameter to retrieve.
 * @return A newly allocated array of three doubles representing the 3D point on
 * success, or a nullptr if the operation failed.
 * @pre The parameter identified by @p name must exist and be of type 3D point.
 * @post The returned array must be deallocated by calling
 * sfcgal_free_buffer()
 * @ingroup capi
 */
SFCGAL_API double *
sfcgal_primitive_parameter_point(const sfcgal_primitive_t *primitive,
                                 const char               *name);

/**
 * @brief Sets the value of a primitive parameter as a 3D point: an array of
 * three doubles (x, y, z).
 * @param primitive Pointer to the primitive.
 * @param name Name of the parameter to set.
 * @param point Array of three doubles representing the new 3D point.
 * @pre The parameter identified by @p name must exist and accept values of type
 * 3D point.
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_primitive_set_parameter_point(sfcgal_primitive_t *primitive,
                                     const char *name, const double *point);

/**
 * @brief Retrieves the value of a primitive parameter as a 3D vector: an array
 * of three doubles (x, y, z).
 * @param primitive Pointer to the primitive.
 * @param name Name of the parameter to retrieve.
 * @return A newly allocated array of three doubles representing the 3D vector
 * on success, or a nullptr if the operation failed.
 * @pre The parameter identified by @p name must exist and be of type 3D vector.
 * @post The returned array must be deallocated by calling
 * sfcgal_free_buffer()
 * @return Primitive parameter as a 3D vector
 * @ingroup capi
 */
SFCGAL_API double *
sfcgal_primitive_parameter_vector(const sfcgal_primitive_t *primitive,
                                  const char               *name);

/**
 * @brief Sets the value of a primitive parameter as a 3D vector: an array of
 * three doubles (x, y, z).
 * @param primitive Pointer to the primitive.
 * @param name Name of the parameter to set.
 * @param vector Array of three doubles representing the new 3D vector.
 * @pre The parameter identified by @p name must exist and accept values of type
 * 3D vector.
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_primitive_set_parameter_vector(sfcgal_primitive_t *primitive,
                                      const char *name, const double *vector);

/**
 * @brief Generates a polyhedral surface representation of the primitive
 * @param primitive Pointer to the primitive.
 * @post The returned geometry must be deallocated by the caller with
 * sfcgal_geometry_delete()
 * @return A SFCGAL::PolyhedralSurface representing the primitive
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_primitive_as_polyhedral_surface(const sfcgal_primitive_t *primitive);

/*--------------------------------------------------------------------------------------*
 *
 * Error handling
 *
 *--------------------------------------------------------------------------------------*/

/**
 * Warning and error handlers
 * @ingroup capi
 */
typedef int (*sfcgal_error_handler_t)(const char *, ...);

/**
 * Sets the error handlers. These callbacks are called on warning or error
 * @param warning_handler is the printf-styled callback function that will be
 * called when a function raises a warning. The default behaviour is to call
 * printf.
 * @param error_handler is the printf-style callback function that will be
 * called when a function generates an error. The default behaviour is to call
 * printf.
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_set_error_handlers(sfcgal_error_handler_t warning_handler,
                          sfcgal_error_handler_t error_handler);

/*--------------------------------------------------------------------------------------*
 *
 * Memory allocation
 *
 *--------------------------------------------------------------------------------------*/

typedef void *(*sfcgal_alloc_handler_t)(size_t);
typedef void (*sfcgal_free_handler_t)(void *);

/**
 * Sets the allocation handlers. These functions are called on memory allocation
 * and deallocation.
 * @param malloc_handler is the function to call for memory allocation. The
 * default behaviour is to call malloc()
 * @param free_handler is the function to call for memory deallocation. The
 * default behaviour is to call free()
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_set_alloc_handlers(sfcgal_alloc_handler_t malloc_handler,
                          sfcgal_free_handler_t  free_handler);

/**
 * Delete a buffer previously allocated and returned by SFCGAL.
 * @param buffer a buffer previously allocated and returned by SFCGAL
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_free_buffer(void *buffer);

/*--------------------------------------------------------------------------------------*
 *
 * Init
 *
 *--------------------------------------------------------------------------------------*/

/**
 * This function must be called before all the other one.
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_init();

/**
 * Get version
 * @ingroup capi
 * @return SFCGAL version as string
 */
SFCGAL_API const char *
sfcgal_version();

/**
 * Get full version (including CGAL and Boost versions)
 * @return SFCGAL version as string including CGAL and Boost versions
 * @ingroup capi
 */
SFCGAL_API const char *
sfcgal_full_version();

/*--------------------------------------------------------------------------------------*
 *
 * NURBSCurve support
 *
 *--------------------------------------------------------------------------------------*/

/**
 * Knot vector generation methods for NURBS curves
 * @ingroup capi
 */
typedef enum {
  SFCGAL_KNOT_METHOD_UNIFORM      = 0,
  SFCGAL_KNOT_METHOD_CHORD_LENGTH = 1,
  SFCGAL_KNOT_METHOD_CENTRIPETAL  = 2
} sfcgal_knot_method_t;

/**
 * End conditions for curve fitting and interpolation
 * @ingroup capi
 */
typedef enum {
  SFCGAL_END_CONDITION_CLAMPED  = 0,
  SFCGAL_END_CONDITION_NATURAL  = 1,
  SFCGAL_END_CONDITION_PERIODIC = 2,
  SFCGAL_END_CONDITION_TANGENT  = 3
} sfcgal_end_condition_t;

/**
 * Curve fitting methods
 * @ingroup capi
 */
typedef enum {
  SFCGAL_FIT_METHOD_INTERPOLATE = 0,
  SFCGAL_FIT_METHOD_APPROXIMATE = 1
} sfcgal_fit_method_t;

/**
 * Creates an empty NURBS curve
 * @return New empty NURBS curve geometry
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_nurbs_curve_create();

/**
 * Creates a NURBS curve from control points with uniform weights
 * @param points Array of control points
 * @param num_points Number of control points
 * @param degree Polynomial degree of curve
 * @param knot_method Method for generating knot vector
 * @return New NURBS curve geometry
 * @post The returned geometry must be deallocated by the caller with @pre
 * sfcgal_geometry_delete
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_nurbs_curve_create_from_points(const sfcgal_geometry_t **points,
                                      size_t num_points, unsigned int degree,
                                      sfcgal_knot_method_t knot_method);

/**
 * Creates a NURBS curve from control points and weights
 * @param points Array of control points
 * @param weights Array of weights (must match num_points)
 * @param num_points Number of control points
 * @param degree Polynomial degree of curve
 * @param knot_method Method for generating knot vector
 * @return New NURBS curve geometry
 * @post The returned geometry must be deallocated by the caller with @pre
 * sfcgal_geometry_delete
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_nurbs_curve_create_from_points_and_weights(
    const sfcgal_geometry_t **points, const double *weights, size_t num_points,
    unsigned int degree, sfcgal_knot_method_t knot_method);

/**
 * Creates a NURBS curve with full control over parameters
 * @param points Array of control points
 * @param weights Array of weights (can be NULL for uniform weights)
 * @param num_points Number of control points
 * @param degree Polynomial degree
 * @param knots Array of knot values
 * @param num_knots Number of knot values (must equal num_points + degree + 1)
 * @return New NURBS curve geometry
 * @post The returned geometry must be deallocated by the caller with @pre
 * sfcgal_geometry_delete
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_nurbs_curve_create_from_full_data(const sfcgal_geometry_t **points,
                                         const double             *weights,
                                         size_t num_points, unsigned int degree,
                                         const double *knots, size_t num_knots);

/**
 * Creates a NURBS curve from Bezier control points
 * @param points Array of Bezier control points
 * @param num_points Number of control points
 * @return New NURBS curve geometry
 * @post The returned geometry must be deallocated by the caller with @pre
 * sfcgal_geometry_delete
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_nurbs_curve_create_bezier(const sfcgal_geometry_t **points,
                                 size_t                    num_points);

/**
 * Creates a uniform B-spline curve
 * @param points Array of control points
 * @param num_points Number of control points
 * @param degree Spline degree
 * @return New NURBS curve geometry
 * @post The returned geometry must be deallocated by the caller with @pre
 * sfcgal_geometry_delete
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_nurbs_curve_create_bspline(const sfcgal_geometry_t **points,
                                  size_t num_points, unsigned int degree);

/**
 * Creates a circular arc as NURBS curve
 * @param center Arc center point
 * @param radius Arc radius (must be positive)
 * @param start_angle Start angle in radians
 * @param end_angle End angle in radians
 * @param normal Normal vector for 3D arcs (can be NULL for Z-axis)
 * @return New NURBS curve geometry representing the circular arc
 * @post The returned geometry must be deallocated by the caller with @pre
 * sfcgal_geometry_delete
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_nurbs_curve_create_circular_arc(const sfcgal_geometry_t *center,
                                       double radius, double start_angle,
                                       double                   end_angle,
                                       const sfcgal_geometry_t *normal);

/**
 * Creates an interpolating NURBS curve through points
 * @param points Array of points to interpolate through
 * @param num_points Number of points
 * @param degree Curve degree
 * @param knot_method Parameterization method
 * @param end_condition Boundary conditions
 * @return New NURBS curve geometry interpolating through the points
 * @post The returned geometry must be deallocated by the caller with @pre
 * sfcgal_geometry_delete
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_nurbs_curve_interpolate(const sfcgal_geometry_t **points,
                               size_t num_points, unsigned int degree,
                               sfcgal_knot_method_t   knot_method,
                               sfcgal_end_condition_t end_condition);

/**
 * Creates an approximating NURBS curve
 * @param points Array of points to approximate
 * @param num_points Number of points
 * @param degree Target curve degree
 * @param tolerance Maximum allowed deviation
 * @param max_control_points Maximum control points to use
 * @return New NURBS curve geometry approximating the points
 * @post The returned geometry must be deallocated by the caller with @pre
 * sfcgal_geometry_delete
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_nurbs_curve_approximate(const sfcgal_geometry_t **points,
                               size_t num_points, unsigned int degree,
                               double tolerance, size_t max_control_points);

/**
 * Unified curve fitting interface
 * @param points Array of points to fit
 * @param num_points Number of points
 * @param degree Target curve degree
 * @param fit_method Whether to interpolate exactly or approximate
 * @param knot_method Parameterization method
 * @param end_condition Boundary conditions (for interpolation)
 * @param tolerance Maximum deviation (for approximation, ignored otherwise)
 * @param max_control_points Maximum control points (for approximation)
 * @return New NURBS curve geometry fitted to the points
 * @post The returned geometry must be deallocated by the caller with @pre
 * sfcgal_geometry_delete
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_nurbs_curve_fit(const sfcgal_geometry_t **points, size_t num_points,
                       unsigned int degree, sfcgal_fit_method_t fit_method,
                       sfcgal_knot_method_t   knot_method,
                       sfcgal_end_condition_t end_condition, double tolerance,
                       size_t max_control_points);

/**
 * Returns the number of control points of a NURBS curve
 * @param curve The NURBS curve geometry
 * @return Number of control points
 * @pre curve must be a NURBS curve
 * @ingroup capi
 */
SFCGAL_API size_t
sfcgal_nurbs_curve_num_control_points(const sfcgal_geometry_t *curve);

/**
 * Returns the ith control point of a NURBS curve
 * @param curve the input NURBS curve
 * @param index control point index
 * @return The control point at the specified index
 * @pre curve must be a NURBS curve
 * @pre i >= 0 and i < sfcgal_nurbs_curve_num_control_points
 * @post the returned Point is not writable and must not be deallocated by the
 * caller
 * @ingroup capi
 */
SFCGAL_API const sfcgal_geometry_t *
sfcgal_nurbs_curve_control_point_n(const sfcgal_geometry_t *curve,
                                   size_t                   index);

/**
 * Sets the ith control point of a NURBS curve
 * @param curve the NURBS curve to modify
 * @param index control point index
 * @param point new control point
 * @pre curve must be a NURBS curve
 * @pre point must be a Point
 * @pre i >= 0 and i < sfcgal_nurbs_curve_num_control_points
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_nurbs_curve_set_control_point_n(sfcgal_geometry_t *curve, size_t index,
                                       const sfcgal_geometry_t *point);

/**
 * Returns the weight at index i of a NURBS curve
 * @param curve the input NURBS curve
 * @param index weight index
 * @return The weight value at the specified index
 * @pre curve must be a NURBS curve
 * @pre i >= 0 and i < sfcgal_nurbs_curve_num_control_points
 * @ingroup capi
 */
SFCGAL_API double
sfcgal_nurbs_curve_weight_n(const sfcgal_geometry_t *curve, size_t index);

/**
 * Sets the weight at index i of a NURBS curve
 * @param curve the NURBS curve to modify
 * @param index weight index
 * @param weight new weight value (must be positive)
 * @pre curve must be a NURBS curve
 * @pre i >= 0 and i < sfcgal_nurbs_curve_num_control_points
 * @ingroup capi
 */
SFCGAL_API void
sfcgal_nurbs_curve_set_weight_n(sfcgal_geometry_t *curve, size_t index,
                                double weight);

/**
 * Returns the degree of a NURBS curve
 * @param curve the input NURBS curve
 * @return The degree of the NURBS curve
 * @pre curve must be a NURBS curve
 * @ingroup capi
 */
SFCGAL_API unsigned int
sfcgal_nurbs_curve_degree(const sfcgal_geometry_t *curve);

/**
 * Returns the number of knots in the knot vector of a NURBS curve
 * @param curve the input NURBS curve
 * @return The number of knots in the knot vector
 * @pre curve must be a NURBS curve
 * @ingroup capi
 */
SFCGAL_API size_t
sfcgal_nurbs_curve_num_knots(const sfcgal_geometry_t *curve);

/**
 * Returns the knot value at index i of a NURBS curve
 * @param curve the input NURBS curve
 * @param index knot index
 * @return The knot value at the specified index
 * @pre curve must be a NURBS curve
 * @pre i >= 0 and i < sfcgal_nurbs_curve_num_knots
 * @ingroup capi
 */
SFCGAL_API double
sfcgal_nurbs_curve_knot_n(const sfcgal_geometry_t *curve, size_t index);

/**
 * Tests if a NURBS curve is rational (uses non-uniform weights)
 * @param curve the input NURBS curve
 * @return 1 if the curve is rational, 0 otherwise
 * @pre curve must be a NURBS curve
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_nurbs_curve_is_rational(const sfcgal_geometry_t *curve);

/**
 * Tests if a NURBS curve is in Bezier form
 * @param curve the input NURBS curve
 * @return 1 if the curve is in Bezier form, 0 otherwise
 * @pre curve must be a NURBS curve
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_nurbs_curve_is_bezier(const sfcgal_geometry_t *curve);

/**
 * Tests if a NURBS curve is a B-spline (uniform weights)
 * @param curve the input NURBS curve
 * @return 1 if the curve is a B-spline, 0 otherwise
 * @pre curve must be a NURBS curve
 * @ingroup capi
 */
SFCGAL_API int
sfcgal_nurbs_curve_is_bspline(const sfcgal_geometry_t *curve);

/**
 * Evaluates a point on the NURBS curve at given parameter
 * @param curve the input NURBS curve
 * @param parameter parameter value (between curve bounds)
 * @return The evaluated point on the curve
 * @pre curve must be a NURBS curve
 * @post The returned Point must be deallocated by the caller with @pre
 * sfcgal_geometry_delete
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_nurbs_curve_evaluate(const sfcgal_geometry_t *curve, double parameter);

/**
 * Computes the derivative of a NURBS curve at given parameter
 * @param curve the input NURBS curve
 * @param parameter parameter value
 * @param order derivative order (1 for first derivative, 2 for second, etc.)
 * @return The derivative vector at the specified parameter
 * @pre curve must be a NURBS curve
 * @post The returned Point must be deallocated by the caller with @pre
 * sfcgal_geometry_delete
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_nurbs_curve_derivative(const sfcgal_geometry_t *curve, double parameter,
                              unsigned int order);

/**
 * Converts a NURBS curve to a LineString with specified number of segments
 * @param curve the input NURBS curve
 * @param num_segments number of line segments
 * @return The converted LineString geometry
 * @pre curve must be a NURBS curve
 * @post The returned LineString must be deallocated by the caller with @pre
 * sfcgal_geometry_delete
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_nurbs_curve_to_linestring(const sfcgal_geometry_t *curve,
                                 unsigned int             num_segments);

/**
 * Converts a NURBS curve to a LineString with adaptive tessellation
 * @param curve the input NURBS curve
 * @param tolerance maximum deviation tolerance
 * @param min_segments minimum number of segments
 * @param max_segments maximum number of segments
 * @return The converted LineString geometry
 * @pre curve must be a NURBS curve
 * @post The returned LineString must be deallocated by the caller with @pre
 * sfcgal_geometry_delete
 * @ingroup capi
 */
SFCGAL_API sfcgal_geometry_t *
sfcgal_nurbs_curve_to_linestring_adaptive(const sfcgal_geometry_t *curve,
                                          double                   tolerance,
                                          unsigned int             min_segments,
                                          unsigned int max_segments);

#ifdef __cplusplus
}
#endif
#endif
