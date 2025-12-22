// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/version.h"

#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/PreparedGeometry.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/primitive3d/Box.h"
#include "SFCGAL/primitive3d/Cone.h"
#include "SFCGAL/primitive3d/Cube.h"
#include "SFCGAL/primitive3d/Cylinder.h"
#include "SFCGAL/primitive3d/Primitive.h"
#include "SFCGAL/primitive3d/Sphere.h"
#include "SFCGAL/primitive3d/Torus.h"
#include "SFCGAL/version.h"

#include "SFCGAL/capi/sfcgal_c.h"

#include "SFCGAL/detail/io/Serialization.h"
#include "SFCGAL/io/OBJ.h"
#include "SFCGAL/io/STL.h"
#include "SFCGAL/io/ewkt.h"
#include "SFCGAL/io/vtk.h"
#include "SFCGAL/io/wkb.h"
#include "SFCGAL/io/wkt.h"
#include "sfcgal_c.h"

#if !defined(_MSC_VER)
  #include "SFCGAL/algorithm/alphaShapes.h"
#endif
#include "SFCGAL/algorithm/alphaWrapping3D.h"
#include "SFCGAL/algorithm/area.h"
#include "SFCGAL/algorithm/buffer3D.h"
#include "SFCGAL/algorithm/centroid.h"
#include "SFCGAL/algorithm/convexHull.h"
#include "SFCGAL/algorithm/covers.h"
#include "SFCGAL/algorithm/difference.h"
#include "SFCGAL/algorithm/distance.h"
#include "SFCGAL/algorithm/distance3d.h"
#include "SFCGAL/algorithm/extrude.h"
#include "SFCGAL/algorithm/force3D.h"
#include "SFCGAL/algorithm/forceMeasured.h"
#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/algorithm/intersects.h"
#include "SFCGAL/algorithm/isClosed.h"
#include "SFCGAL/algorithm/isSimple.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/algorithm/length.h"
#include "SFCGAL/algorithm/lineSubstring.h"
#include "SFCGAL/algorithm/minkowskiSum.h"
#include "SFCGAL/algorithm/offset.h"
#include "SFCGAL/algorithm/partition_2.h"
#include "SFCGAL/algorithm/plane.h"
#if SFCGAL_CGAL_VERSION_MAJOR >= 6
  #include "SFCGAL/algorithm/polygonRepair.h"
#endif
#include "SFCGAL/algorithm/rotate.h"
#include "SFCGAL/algorithm/scale.h"
#include "SFCGAL/algorithm/simplification.h"
#include "SFCGAL/algorithm/straightSkeleton.h"
#include "SFCGAL/algorithm/tesselate.h"
#include "SFCGAL/algorithm/translate.h"
#include "SFCGAL/algorithm/union.h"
#include "SFCGAL/algorithm/visibility.h"
#include "SFCGAL/algorithm/volume.h"
#include "SFCGAL/detail/transform/AffineTransform3.h"
#include "SFCGAL/triangulate/triangulate2DZ.h"

#include "SFCGAL/detail/transform/ForceOrderPoints.h"
#include "SFCGAL/detail/transform/ForceZOrderPoints.h"
#include "SFCGAL/detail/transform/RoundTransform.h"
#include <cmath>

#include <nlohmann/json.hpp>

//
// Note about sfcgal_geometry_t pointers: they are basically void* pointers that
// represent pointers to a SFCGAL::Geometry. In order to support multiple
// inheritance: every input or output sfcgal_geometry_t* is a pointer to the
// *base class* SFCGAL::Geometry. If a function wants to return a sub-class, it
// must be up casted before returned (static_cast) If a function wants to use a
// sub-class from a parameter, it must also be down casted. For instance,
// static_cast<SFCGAL::Point*>(reinterpret_cast<SFCGAL::Geometry*>(p))
//
// SFCGAL::PreparedGeometry has no vtable and can thus be manipuled through
// reinterpret_cast without problem

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection
#ifndef DOXYGEN_SHOULD_SKIP_THIS

static sfcgal_error_handler_t __sfcgal_warning_handler = printf;
static sfcgal_error_handler_t __sfcgal_error_handler   = printf;

  #define SFCGAL_WARNING __sfcgal_warning_handler
  #define SFCGAL_ERROR __sfcgal_error_handler

  #define SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(call)                         \
    try {                                                                      \
      call                                                                     \
    } catch (std::exception & e) {                                             \
      SFCGAL_ERROR("%s", e.what());                                            \
      return 0;                                                                \
    }

  #define SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(call)                  \
    try {                                                                      \
      call                                                                     \
    } catch (std::exception & e) {                                             \
      SFCGAL_ERROR("%s", e.what());                                            \
    }

  // Functions that take two geometries and return a scalar
  //
  // name: C function name
  // ret_type: C function return type
  // sfcgal_function: C++ SFCGAL method to call
  // cpp_type: C++ return type (might be different than ret_type)
  // fail_value: returned value on failure
  #define SFCGAL_GEOMETRY_FUNCTION_BINARY_SCALAR(                              \
      name, sfcgal_function, ret_type, cpp_type, fail_value)                   \
    extern "C" ret_type sfcgal_geometry_##name(const sfcgal_geometry_t *ga,    \
                                               const sfcgal_geometry_t *gb)    \
    {                                                                          \
      cpp_type result;                                                         \
      try {                                                                    \
        result = sfcgal_function(*static_cast<const SFCGAL::Geometry *>(ga),   \
                                 *static_cast<const SFCGAL::Geometry *>(gb));  \
      } catch (std::exception & e) {                                           \
        SFCGAL_WARNING("During " #name "(A,B) :");                             \
        SFCGAL_WARNING(                                                        \
            "  with A: %s",                                                    \
            static_cast<const SFCGAL::Geometry *>(ga)->asText().c_str());      \
        SFCGAL_WARNING(                                                        \
            "   and B: %s",                                                    \
            static_cast<const SFCGAL::Geometry *>(gb)->asText().c_str());      \
        SFCGAL_ERROR("%s", e.what());                                          \
        return fail_value;                                                     \
      }                                                                        \
      return result;                                                           \
    }

  #define SFCGAL_GEOMETRY_FUNCTION_BINARY_PREDICATE(name, sfcgal_function)     \
    SFCGAL_GEOMETRY_FUNCTION_BINARY_SCALAR(name, sfcgal_function, int, bool, -1)

  #define SFCGAL_GEOMETRY_FUNCTION_BINARY_MEASURE(name, sfcgal_function)       \
    SFCGAL_GEOMETRY_FUNCTION_BINARY_SCALAR(name, sfcgal_function, double,      \
                                           double, -1.0)

  #define SFCGAL_GEOMETRY_FUNCTION_BINARY_CONSTRUCTION(name, sfcgal_function)  \
    extern "C" sfcgal_geometry_t *sfcgal_geometry_##name(                      \
        const sfcgal_geometry_t *ga, const sfcgal_geometry_t *gb)              \
    {                                                                          \
      std::unique_ptr<SFCGAL::Geometry> result;                                \
      try {                                                                    \
        result = sfcgal_function(*static_cast<const SFCGAL::Geometry *>(ga),   \
                                 *static_cast<const SFCGAL::Geometry *>(gb));  \
      } catch (std::exception & e) {                                           \
        SFCGAL_WARNING("During " #name "(A,B) :");                             \
        SFCGAL_WARNING(                                                        \
            "  with A: %s",                                                    \
            static_cast<const SFCGAL::Geometry *>(ga)->asText().c_str());      \
        SFCGAL_WARNING(                                                        \
            "   and B: %s",                                                    \
            static_cast<const SFCGAL::Geometry *>(gb)->asText().c_str());      \
        SFCGAL_ERROR("%s", e.what());                                          \
        return 0;                                                              \
      }                                                                        \
      return result.release();                                                 \
    }

  #define SFCGAL_GEOMETRY_FUNCTION_UNARY_CONSTRUCTION(name, sfcgal_function)   \
    extern "C" sfcgal_geometry_t *sfcgal_geometry_##name(                      \
        const sfcgal_geometry_t *ga)                                           \
    {                                                                          \
      std::unique_ptr<SFCGAL::Geometry> result;                                \
      try {                                                                    \
        result = sfcgal_function(*static_cast<const SFCGAL::Geometry *>(ga));  \
      } catch (std::exception & e) {                                           \
        SFCGAL_WARNING("During " #name "(A) :");                               \
        SFCGAL_WARNING(                                                        \
            "  with A: %s",                                                    \
            static_cast<const SFCGAL::Geometry *>(ga)->asText().c_str());      \
        SFCGAL_ERROR("%s", e.what());                                          \
        return 0;                                                              \
      }                                                                        \
      return result.release();                                                 \
    }

  #define SFCGAL_GEOMETRY_FUNCTION_UNARY_MEASURE(name, sfcgal_function)        \
    extern "C" double sfcgal_geometry_##name(const sfcgal_geometry_t *ga)      \
    {                                                                          \
      double result;                                                           \
      try {                                                                    \
        result = sfcgal_function(*static_cast<const SFCGAL::Geometry *>(ga));  \
      } catch (std::exception & e) {                                           \
        SFCGAL_WARNING("During " #name "(A) :");                               \
        SFCGAL_WARNING(                                                        \
            "  with A: %s",                                                    \
            static_cast<const SFCGAL::Geometry *>(ga)->asText().c_str());      \
        SFCGAL_ERROR("%s", e.what());                                          \
        return -1.0;                                                           \
      }                                                                        \
      return result;                                                           \
    }

template <class T>
inline auto
down_cast(sfcgal_geometry_t *p) -> T *
{
  T *result = dynamic_cast<T *>(static_cast<SFCGAL::Geometry *>(p));

  if (!result) {
    BOOST_THROW_EXCEPTION(SFCGAL::Exception("wrong geometry type"));
  }

  return result;
}

template <class T>
inline auto
down_const_cast(const sfcgal_geometry_t *p) -> const T *
{
  const T *q =
      dynamic_cast<const T *>(static_cast<const SFCGAL::Geometry *>(p));

  if (!q) {
    BOOST_THROW_EXCEPTION(SFCGAL::Exception("wrong geometry type"));
  }

  return q;
}

static sfcgal_alloc_handler_t sfcgal_alloc_handler = malloc;
static sfcgal_free_handler_t  sfcgal_free_handler  = free;

inline auto
alloc_and_copy(const std::string &str, char **buffer, size_t *len) -> void
{
  *len    = str.size();
  *buffer = (char *)sfcgal_alloc_handler(*len + 1);
  if (*buffer != nullptr) {
    memset(*buffer, 0, *len + 1);
    memcpy(*buffer, str.data(), *len);
  } else {
    *len = 0;
  }
}
#endif // ifndef DOXYGEN_SHOULD_SKIP_THIS
/// @} end of private section
// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

extern "C" void
sfcgal_set_error_handlers(sfcgal_error_handler_t warning_handler,
                          sfcgal_error_handler_t error_handler)
{
  __sfcgal_warning_handler = warning_handler;
  __sfcgal_error_handler   = error_handler;
}

extern "C" void
sfcgal_set_alloc_handlers(sfcgal_alloc_handler_t alloc_handler,
                          sfcgal_free_handler_t  free_handler)
{
  sfcgal_alloc_handler = alloc_handler;
  sfcgal_free_handler  = free_handler;
}

extern "C" void
sfcgal_free_buffer(void *buffer)
{
  if (buffer != nullptr) {
    sfcgal_free_handler(buffer);
  }
}

extern "C" void
sfcgal_init(void)
{
  set_error_behaviour(CGAL::THROW_EXCEPTION);
  set_warning_behaviour(CGAL::THROW_EXCEPTION);
}

extern "C" auto
sfcgal_version(void) -> const char *
{
  return SFCGAL_Version();
}

extern "C" auto
sfcgal_full_version(void) -> const char *
{
  return SFCGAL_Full_Version();
}

extern "C" void
sfcgal_set_geometry_validation(int /*enabled*/)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      BOOST_THROW_EXCEPTION(SFCGAL::Exception("Not implemented")););
}

extern "C" auto
sfcgal_geometry_type_id(const sfcgal_geometry_t *geom) -> sfcgal_geometry_type_t
{

  try {
    return (sfcgal_geometry_type_t) reinterpret_cast<const SFCGAL::Geometry *>(
               geom)
        ->geometryTypeId();
  } catch (std::exception &e) {
    SFCGAL_ERROR("%s", e.what());
    return SFCGAL_TYPE_POINT; // to avoid warning
  }
}

extern "C" auto
sfcgal_geometry_type(const sfcgal_geometry_t *geom, char **type,
                     size_t *typeLen) -> void
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      std::string wkt =
          reinterpret_cast<const SFCGAL::Geometry *>(geom)->geometryType();
      alloc_and_copy(wkt, type, typeLen);)
}

extern "C" auto
sfcgal_geometry_dimension(const sfcgal_geometry_t *geom) -> int
{

  try {
    return reinterpret_cast<const SFCGAL::Geometry *>(geom)->dimension();
  } catch (std::exception &e) {
    SFCGAL_ERROR("%s", e.what());
    return 0; // to avoid warning
  }
}

extern "C" auto
sfcgal_geometry_is_valid(const sfcgal_geometry_t *geom) -> int
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return (int)bool(SFCGAL::algorithm::isValid(
          *reinterpret_cast<const SFCGAL::Geometry *>(geom)));)
}

// deprecated!
extern "C" auto
sfcgal_geometry_is_complexity_detail(const sfcgal_geometry_t *geom,
                                     char                   **invalidity_reason,
                                     sfcgal_geometry_t **invalidity_location)
    -> int
{
  return sfcgal_geometry_is_valid_detail(geom, invalidity_reason,
                                         invalidity_location);
}

extern "C" auto
sfcgal_geometry_is_valid_detail(const sfcgal_geometry_t *geom,
                                char                   **invalidity_reason,
                                sfcgal_geometry_t **invalidity_location) -> int
{
  // invalidity location is not supported for now
  if (invalidity_location != nullptr) {
    *invalidity_location = nullptr;
  }
  // set to null for now
  if (invalidity_reason != nullptr) {
    *invalidity_reason = nullptr;
  }

  const auto *g = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  if (g->hasValidityFlag()) {
    return 1;
  }
  bool is_valid = false;
  try {
    SFCGAL::Validity const validity = SFCGAL::algorithm::isValid(*g);
    is_valid                        = validity;
    if (!is_valid && (invalidity_reason != nullptr)) {
      *invalidity_reason = strdup(validity.reason().c_str());
    }
  } catch (SFCGAL::Exception &e) {
    if (invalidity_reason != nullptr) {
      *invalidity_reason = strdup(e.what());
    }
  }
  return static_cast<int>(is_valid);
}

extern "C" auto
sfcgal_geometry_is_simple(const sfcgal_geometry_t *geom) -> int
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return (int)bool(SFCGAL::algorithm::isSimple(
          *reinterpret_cast<const SFCGAL::Geometry *>(geom)));)
}

extern "C" auto
sfcgal_geometry_is_simple_detail(const sfcgal_geometry_t *geom,
                                 char **complexity_reason) -> int
{
  // set to null for now
  if (complexity_reason != nullptr) {
    *complexity_reason = nullptr;
  }

  const auto *g         = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  bool        is_simple = false;
  try {
    SFCGAL::Simplicity const simplicity = SFCGAL::algorithm::isSimple(*g);
    is_simple                           = simplicity;
    if (!is_simple && (complexity_reason != nullptr)) {
      *complexity_reason = strdup(simplicity.reason().c_str());
    }
  } catch (SFCGAL::Exception &e) {
    if (complexity_reason != nullptr) {
      *complexity_reason = strdup(e.what());
    }
  }
  return static_cast<int>(is_simple);
}

extern "C" auto
sfcgal_geometry_is_3d(const sfcgal_geometry_t *geom) -> int
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return (int)reinterpret_cast<const SFCGAL::Geometry *>(geom)->is3D();)
}

extern "C" auto
sfcgal_geometry_is_measured(const sfcgal_geometry_t *geom) -> int
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return (int)reinterpret_cast<const SFCGAL::Geometry *>(geom)
          ->isMeasured();)
}

extern "C" auto
sfcgal_geometry_is_empty(const sfcgal_geometry_t *geom) -> int
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return (int)reinterpret_cast<const SFCGAL::Geometry *>(geom)->isEmpty();)
}

extern "C" auto
sfcgal_geometry_is_closed(const sfcgal_geometry_t *geom) -> int
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return (int)bool(SFCGAL::algorithm::isClosed(
          *reinterpret_cast<const SFCGAL::Geometry *>(geom)));)
}

extern "C" auto
sfcgal_geometry_drop_z(sfcgal_geometry_t *geom) -> int
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return (int)reinterpret_cast<SFCGAL::Geometry *>(geom)->dropZ();)
}

extern "C" auto
sfcgal_geometry_drop_m(sfcgal_geometry_t *geom) -> int
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return (int)reinterpret_cast<SFCGAL::Geometry *>(geom)->dropM();)
}

extern "C" auto
sfcgal_geometry_force_z(sfcgal_geometry_t *geom, double defaultZ) -> int
{
  auto *sfcgalGeom = reinterpret_cast<SFCGAL::Geometry *>(geom);
  // If the geometry is already 3D or empty, don't do anything
  if (sfcgalGeom->isEmpty() || sfcgalGeom->is3D()) {
    return 0;
  }

  try {
    SFCGAL::algorithm::force3D(*sfcgalGeom, defaultZ);
  } catch (std::exception &e) {
    SFCGAL_ERROR("%s", e.what());
    return 0;
  }

  return 1;
}

extern "C" auto
sfcgal_geometry_force_m(sfcgal_geometry_t *geom, double defaultM) -> int
{
  auto *sfcgalGeom = reinterpret_cast<SFCGAL::Geometry *>(geom);
  // If the geometry is already measured or empty, don't do anything
  if (sfcgalGeom->isEmpty() || sfcgalGeom->isMeasured()) {
    return 0;
  }

  try {
    SFCGAL::algorithm::forceMeasured(*sfcgalGeom, defaultM);
  } catch (std::exception &e) {
    SFCGAL_ERROR("%s", e.what());
    return 0;
  }

  return 1;
}

extern "C" auto
sfcgal_geometry_swap_xy(sfcgal_geometry_t *geom) -> void
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      return reinterpret_cast<SFCGAL::Geometry *>(geom)->swapXY();)
}

extern "C" auto
sfcgal_geometry_clone(const sfcgal_geometry_t *geom) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return reinterpret_cast<const SFCGAL::Geometry *>(geom)
          ->clone()
          .release();)
}

extern "C" void
sfcgal_geometry_delete(sfcgal_geometry_t *geom)
{
  delete reinterpret_cast<SFCGAL::Geometry *>(geom);
}

extern "C" auto
sfcgal_geometry_num_geometries(const sfcgal_geometry_t *geometryCollection)
    -> size_t
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return down_const_cast<SFCGAL::Geometry>(geometryCollection)
          ->numGeometries();)
}

extern "C" void
sfcgal_geometry_as_text(const sfcgal_geometry_t *pgeom, char **buffer,
                        size_t *len)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      std::string wkt =
          reinterpret_cast<const SFCGAL::Geometry *>(pgeom)->asText();
      alloc_and_copy(wkt, buffer, len);)
}

extern "C" void
sfcgal_geometry_as_text_decim(const sfcgal_geometry_t *pgeom, int numDecimals,
                              char **buffer, size_t *len)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      std::string wkt =
          reinterpret_cast<const SFCGAL::Geometry *>(pgeom)->asText(
              numDecimals);
      alloc_and_copy(wkt, buffer, len);)
}

extern "C" void
sfcgal_geometry_as_wkb(const sfcgal_geometry_t *pgeom, char **buffer,
                       size_t *len)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      std::string wkb =
          reinterpret_cast<const SFCGAL::Geometry *>(pgeom)->asWkb(
              boost::endian::order::native, false);
      alloc_and_copy(wkb, buffer, len);)
}

extern "C" void
sfcgal_geometry_as_hexwkb(const sfcgal_geometry_t *pgeom, char **buffer,
                          size_t *len)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      std::string wkb =
          reinterpret_cast<const SFCGAL::Geometry *>(pgeom)->asWkb(
              boost::endian::order::native, true);
      alloc_and_copy(wkb, buffer, len);)
}

extern "C" auto
sfcgal_geometry_as_vtk_file(const sfcgal_geometry_t *pgeom,
                            const char              *filename) -> void
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      SFCGAL::io::VTK::save(*reinterpret_cast<const SFCGAL::Geometry *>(pgeom),
                            filename);)
}

extern "C" auto
sfcgal_geometry_as_vtk(const sfcgal_geometry_t *pgeom, char **buffer,
                       size_t *len) -> void
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      std::string obj = SFCGAL::io::VTK::saveToString(
          *reinterpret_cast<const SFCGAL::Geometry *>(pgeom));
      alloc_and_copy(obj, buffer, len);)
}

extern "C" auto
sfcgal_geometry_as_stl_file(const sfcgal_geometry_t *pgeom,
                            const char              *filename) -> void
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      SFCGAL::io::STL::save(*reinterpret_cast<const SFCGAL::Geometry *>(pgeom),
                            filename);)
}

extern "C" auto
sfcgal_geometry_as_stl(const sfcgal_geometry_t *pgeom, char **buffer,
                       size_t *len) -> void
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      std::string obj = SFCGAL::io::STL::saveToString(
          *reinterpret_cast<const SFCGAL::Geometry *>(pgeom));
      alloc_and_copy(obj, buffer, len);)
}

extern "C" void
sfcgal_geometry_as_obj_file(const sfcgal_geometry_t *pgeom,
                            const char              *filename)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      SFCGAL::io::OBJ::save(*reinterpret_cast<const SFCGAL::Geometry *>(pgeom),
                            filename);)
}

extern "C" void
sfcgal_geometry_as_obj(const sfcgal_geometry_t *pgeom, char **buffer,
                       size_t *len)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      std::string obj = SFCGAL::io::OBJ::saveToString(
          *reinterpret_cast<const SFCGAL::Geometry *>(pgeom));
      alloc_and_copy(obj, buffer, len);)
}

/**
 * Point
 */
extern "C" auto
sfcgal_point_create(void) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<SFCGAL::Geometry *>(new SFCGAL::Point());)
}

extern "C" auto
sfcgal_point_create_from_xy(double x, double y) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<SFCGAL::Geometry *>(new SFCGAL::Point(x, y));)
}

extern "C" auto
sfcgal_point_create_from_xym(double x, double y, double m)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      auto g = std::make_unique<SFCGAL::Point>(x, y); g->setM(m);
      return static_cast<SFCGAL::Geometry *>(g.release());)
}

extern "C" auto
sfcgal_point_create_from_xyz(double x, double y, double z)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<SFCGAL::Geometry *>(new SFCGAL::Point(x, y, z));)
}

extern "C" auto
sfcgal_point_create_from_xyzm(double x, double y, double z, double m)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<SFCGAL::Geometry *>(new SFCGAL::Point(x, y, z, m));)
}

extern "C" auto
sfcgal_point_x(const sfcgal_geometry_t *geom) -> double
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return CGAL::to_double(down_const_cast<SFCGAL::Point>(geom)->x());)
}

extern "C" auto
sfcgal_point_y(const sfcgal_geometry_t *geom) -> double
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return CGAL::to_double(down_const_cast<SFCGAL::Point>(geom)->y());)
}

extern "C" auto
sfcgal_point_z(const sfcgal_geometry_t *geom) -> double
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return CGAL::to_double(down_const_cast<SFCGAL::Point>(geom)->z());)
}

extern "C" auto
sfcgal_point_m(const sfcgal_geometry_t *geom) -> double
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return CGAL::to_double(down_const_cast<SFCGAL::Point>(geom)->m());)
}

/**
 * LineString
 */
extern "C" auto
sfcgal_linestring_create(void) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<SFCGAL::Geometry *>(new SFCGAL::LineString());)
}

extern "C" auto
sfcgal_linestring_num_points(const sfcgal_geometry_t *geom) -> size_t
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return down_const_cast<SFCGAL::LineString>(geom)->numPoints();)
}

extern "C" auto
sfcgal_linestring_point_n(const sfcgal_geometry_t *geom, size_t i)
    -> const sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<const SFCGAL::Geometry *>(
                 &(down_const_cast<SFCGAL::LineString>(geom)->pointN(i)));)
}

extern "C" void
sfcgal_linestring_add_point(sfcgal_geometry_t *geom, sfcgal_geometry_t *point)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      down_cast<SFCGAL::LineString>(geom)->addPoint(
          down_cast<SFCGAL::Point>(point));)
}

extern "C" void
sfcgal_linestring_closes(sfcgal_geometry_t *linestring)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      down_cast<SFCGAL::LineString>(linestring)->closes();)
}

/**
 * Triangle
 */
extern "C" auto
sfcgal_triangle_create(void) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<SFCGAL::Geometry *>(new SFCGAL::Triangle());)
}

extern "C" auto
sfcgal_triangle_create_from_points(const sfcgal_geometry_t *pa,
                                   const sfcgal_geometry_t *pb,
                                   const sfcgal_geometry_t *pc)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<SFCGAL::Geometry *>(
                 new SFCGAL::Triangle(*down_const_cast<SFCGAL::Point>(pa),
                                      *down_const_cast<SFCGAL::Point>(pb),
                                      *down_const_cast<SFCGAL::Point>(pc)));)
}

extern "C" auto
sfcgal_triangle_vertex(const sfcgal_geometry_t *geom, int i)
    -> const sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<const SFCGAL::Geometry *>(
                 &down_const_cast<SFCGAL::Triangle>(geom)->vertex(i));)
}

extern "C" void
sfcgal_triangle_set_vertex(sfcgal_geometry_t *geom, int i,
                           const sfcgal_geometry_t *point)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      down_cast<SFCGAL::Triangle>(geom)->vertex(i) =
          *down_const_cast<const SFCGAL::Point>(point);)
}

extern "C" void
sfcgal_triangle_set_vertex_from_xy(sfcgal_geometry_t *geom, int i, double x,
                                   double y)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      down_cast<SFCGAL::Triangle>(geom)->vertex(i) = SFCGAL::Point(x, y);)
}

extern "C" void
sfcgal_triangle_set_vertex_from_xyz(sfcgal_geometry_t *geom, int i, double x,
                                    double y, double z)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      down_cast<SFCGAL::Triangle>(geom)->vertex(i) = SFCGAL::Point(x, y, z);)
}

/**
 * Polygon
 */
extern "C" auto
sfcgal_polygon_create(void) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<SFCGAL::Geometry *>(new SFCGAL::Polygon());)
}

extern "C" auto
sfcgal_polygon_create_from_exterior_ring(sfcgal_geometry_t *ring)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<SFCGAL::Geometry *>(
                 new SFCGAL::Polygon(down_cast<SFCGAL::LineString>(ring)));)
}

extern "C" void
sfcgal_polygon_set_exterior_ring(sfcgal_geometry_t *polygon,
                                 sfcgal_geometry_t *ring)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      down_cast<SFCGAL::Polygon>(polygon)->setExteriorRing(
          down_cast<SFCGAL::LineString>(ring));)
}

extern "C" auto
sfcgal_polygon_exterior_ring(const sfcgal_geometry_t *geom)
    -> const sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<const SFCGAL::Geometry *>(
                 &down_const_cast<SFCGAL::Polygon>(geom)->exteriorRing());)
}

extern "C" auto
sfcgal_polygon_num_interior_rings(const sfcgal_geometry_t *geom) -> size_t
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return down_const_cast<SFCGAL::Polygon>(geom)->numInteriorRings();)
}

extern "C" auto
sfcgal_polygon_interior_ring_n(const sfcgal_geometry_t *geom, size_t i)
    -> const sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<const SFCGAL::Geometry *>(
                 &down_const_cast<SFCGAL::Polygon>(geom)->interiorRingN(i));)
}

extern "C" void
sfcgal_polygon_add_interior_ring(sfcgal_geometry_t *geom,
                                 sfcgal_geometry_t *ring)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      down_cast<SFCGAL::Polygon>(geom)->addRing(
          down_cast<SFCGAL::LineString>(ring));)
}

/**
 * Geometry collection
 */

extern "C" auto
sfcgal_geometry_collection_create(void) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<SFCGAL::Geometry *>(new SFCGAL::GeometryCollection());)
}

extern "C" auto
sfcgal_geometry_collection_num_geometries(const sfcgal_geometry_t *geom)
    -> size_t
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return down_const_cast<SFCGAL::GeometryCollection>(geom)
          ->numGeometries();)
}

extern "C" auto
sfcgal_geometry_collection_geometry_n(const sfcgal_geometry_t *geom, size_t i)
    -> const sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const auto *g = down_const_cast<SFCGAL::GeometryCollection>(geom);
      return static_cast<const SFCGAL::Geometry *>(&g->geometryN(i));)
}

extern "C" auto
sfcgal_geometry_collection_set_geometry_n(sfcgal_geometry_t *collection,
                                          sfcgal_geometry_t *geometry, size_t i)
    -> void
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      return down_cast<SFCGAL::GeometryCollection>(collection)
          ->setGeometryN(reinterpret_cast<SFCGAL::Geometry *>(geometry), i);)
}

extern "C" void
sfcgal_geometry_collection_add_geometry(sfcgal_geometry_t *geom,
                                        sfcgal_geometry_t *ngeom)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      down_cast<SFCGAL::GeometryCollection>(geom)->addGeometry(
          reinterpret_cast<SFCGAL::Geometry *>(ngeom));)
}

/**
 * Multi-*
 */
extern "C" auto
sfcgal_multi_point_create(void) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<SFCGAL::Geometry *>(new SFCGAL::MultiPoint());)
}

extern "C" auto
sfcgal_multi_linestring_create(void) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<SFCGAL::Geometry *>(new SFCGAL::MultiLineString());)
}

extern "C" auto
sfcgal_multi_polygon_create(void) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<SFCGAL::Geometry *>(new SFCGAL::MultiPolygon());)
}

extern "C" auto
sfcgal_multi_solid_create(void) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<SFCGAL::Geometry *>(new SFCGAL::MultiSolid());)
}

/**
 * Polyhedral surface
 */

extern "C" auto
sfcgal_polyhedral_surface_create(void) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<SFCGAL::Geometry *>(new SFCGAL::PolyhedralSurface());)
}

extern "C" auto
sfcgal_polyhedral_surface_num_patches(const sfcgal_geometry_t *polyhedral)
    -> size_t
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return down_const_cast<SFCGAL::PolyhedralSurface>(polyhedral)
          ->numPatches();)
}

extern "C" auto
sfcgal_polyhedral_surface_num_polygons(const sfcgal_geometry_t *geom) -> size_t
{
  return sfcgal_polyhedral_surface_num_patches(geom);
}

extern "C" auto
sfcgal_polyhedral_surface_patch_n(const sfcgal_geometry_t *geom, size_t i)
    -> const sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<const SFCGAL::Geometry *>(
                 &down_const_cast<SFCGAL::PolyhedralSurface>(geom)->patchN(i));)
}

extern "C" auto
sfcgal_polyhedral_surface_polygon_n(const sfcgal_geometry_t *geom, size_t i)
    -> const sfcgal_geometry_t *
{
  return sfcgal_polyhedral_surface_patch_n(geom, i);
}

extern "C" void
sfcgal_polyhedral_surface_add_patch(sfcgal_geometry_t *geom,
                                    sfcgal_geometry_t *patch)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      return down_cast<SFCGAL::PolyhedralSurface>(geom)->addPatch(
          down_cast<SFCGAL::Polygon>(patch));)
}

extern "C" void
sfcgal_polyhedral_surface_set_patch_n(sfcgal_geometry_t *polyhedral,
                                      sfcgal_geometry_t *patch, size_t i)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      return down_cast<SFCGAL::PolyhedralSurface>(polyhedral)
          ->setPatchN(down_cast<SFCGAL::Polygon>(patch), i);)
}

extern "C" void
sfcgal_polyhedral_surface_add_polygon(sfcgal_geometry_t *geom,
                                      sfcgal_geometry_t *poly)
{
  sfcgal_polyhedral_surface_add_patch(geom, poly);
}

/**
 * Triangulated surface
 */

extern "C" auto
sfcgal_triangulated_surface_create(void) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<SFCGAL::Geometry *>(
                 new SFCGAL::TriangulatedSurface());)
}

extern "C" auto
sfcgal_triangulated_surface_num_patches(const sfcgal_geometry_t *tin) -> size_t
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return down_const_cast<SFCGAL::TriangulatedSurface>(tin)->numPatches();)
}

extern "C" auto
sfcgal_triangulated_surface_num_triangles(const sfcgal_geometry_t *geom)
    -> size_t
{
  return sfcgal_triangulated_surface_num_patches(geom);
}

extern "C" auto
sfcgal_triangulated_surface_patch_n(const sfcgal_geometry_t *tin, size_t i)
    -> const sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<const SFCGAL::Geometry *>(
                 &down_const_cast<SFCGAL::TriangulatedSurface>(tin)->patchN(
                     i));)
}

extern "C" auto
sfcgal_triangulated_surface_triangle_n(const sfcgal_geometry_t *geom, size_t i)
    -> const sfcgal_geometry_t *
{
  return sfcgal_triangulated_surface_patch_n(geom, i);
}

extern "C" void
sfcgal_triangulated_surface_set_patch_n(sfcgal_geometry_t *tin,
                                        sfcgal_geometry_t *patch, size_t i)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      return down_cast<SFCGAL::TriangulatedSurface>(tin)->setPatchN(
          down_cast<SFCGAL::Triangle>(patch), i);)
}

extern "C" void
sfcgal_triangulated_surface_add_patch(sfcgal_geometry_t *tin,
                                      sfcgal_geometry_t *patch)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      down_cast<SFCGAL::TriangulatedSurface>(tin)->addPatch(
          down_cast<SFCGAL::Triangle>(patch));)
}

extern "C" void
sfcgal_triangulated_surface_add_triangle(sfcgal_geometry_t *geom,
                                         sfcgal_geometry_t *triangle)
{
  sfcgal_triangulated_surface_add_patch(geom, triangle);
}

/**
 * Solid
 */

extern "C" auto
sfcgal_solid_create(void) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<SFCGAL::Geometry *>(new SFCGAL::Solid());)
}

extern "C" auto
sfcgal_solid_create_from_exterior_shell(sfcgal_geometry_t *shell)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<SFCGAL::Geometry *>(new SFCGAL::Solid(
          down_cast<SFCGAL::PolyhedralSurface>(shell)));)
}

extern "C" auto
sfcgal_solid_num_shells(const sfcgal_geometry_t *geom) -> size_t
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return down_const_cast<SFCGAL::Solid>(geom)->numShells();)
}

extern "C" auto
sfcgal_solid_shell_n(const sfcgal_geometry_t *geom, size_t i)
    -> const sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<const SFCGAL::Geometry *>(
                 &down_const_cast<SFCGAL::Solid>(geom)->shellN(i));)
}

extern "C" void
sfcgal_solid_add_interior_shell(sfcgal_geometry_t *geom,
                                sfcgal_geometry_t *shell)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      down_cast<SFCGAL::Solid>(geom)->addInteriorShell(
          down_cast<SFCGAL::PolyhedralSurface>(shell));)
}

extern "C" void
sfcgal_solid_set_exterior_shell(sfcgal_geometry_t *geom,
                                sfcgal_geometry_t *shell)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      down_cast<SFCGAL::Solid>(geom)->setExteriorShell(
          down_cast<SFCGAL::PolyhedralSurface>(shell));)
}

extern "C" auto
sfcgal_prepared_geometry_create(void) -> sfcgal_prepared_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(return new SFCGAL::PreparedGeometry();)
}

extern "C" auto
sfcgal_prepared_geometry_create_from_geometry(sfcgal_geometry_t *geom,
                                              srid_t             srid)
    -> sfcgal_prepared_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return new SFCGAL::PreparedGeometry(
                 reinterpret_cast<SFCGAL::Geometry *>(geom), srid);)
}

extern "C" void
sfcgal_prepared_geometry_delete(sfcgal_prepared_geometry_t *pgeom)
{
  delete reinterpret_cast<SFCGAL::PreparedGeometry *>(pgeom);
}

extern "C" auto
sfcgal_prepared_geometry_geometry(const sfcgal_prepared_geometry_t *pgeom)
    -> const sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return &reinterpret_cast<const SFCGAL::PreparedGeometry *>(pgeom)
                  ->geometry();)
}

extern "C" void
sfcgal_prepared_geometry_set_geometry(sfcgal_prepared_geometry_t *pgeom,
                                      sfcgal_geometry_t          *geom)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      reinterpret_cast<SFCGAL::PreparedGeometry *>(pgeom)->resetGeometry(
          reinterpret_cast<SFCGAL::Geometry *>(geom));)
}

extern "C" auto
sfcgal_prepared_geometry_srid(const sfcgal_prepared_geometry_t *pgeom) -> srid_t
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return reinterpret_cast<const SFCGAL::PreparedGeometry *>(pgeom)->SRID();)
}

extern "C" void
sfcgal_prepared_geometry_set_srid(sfcgal_prepared_geometry_t *pgeom,
                                  srid_t                      srid)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      reinterpret_cast<SFCGAL::PreparedGeometry *>(pgeom)->SRID() = srid;)
}

extern "C" void
sfcgal_prepared_geometry_as_ewkt(const sfcgal_prepared_geometry_t *pgeom,
                                 int num_decimals, char **buffer, size_t *len)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      std::string ewkt =
          reinterpret_cast<const SFCGAL::PreparedGeometry *>(pgeom)->asEWKT(
              num_decimals);
      alloc_and_copy(ewkt, buffer, len);)
}

extern "C" auto
sfcgal_io_read_wkt(const char *str, size_t len) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return SFCGAL::io::readWkt(str, len).release();)
}

extern "C" auto
sfcgal_io_read_wkb(const char *str, size_t len) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      if (len > 2 && str[0] == '0' &&
          (str[1] == '0' || str[1] == '1')) return SFCGAL::io::readWkb(str, len,
                                                                       true)
          .release();
      return SFCGAL::io::readWkb(str, len, false).release();)
}

extern "C" void
sfcgal_io_write_binary_prepared(const sfcgal_prepared_geometry_t *geom,
                                char **buffer, size_t *len)
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      const auto *g = reinterpret_cast<const SFCGAL::PreparedGeometry *>(geom);
      std::string str = SFCGAL::io::writeBinaryPrepared(*g);
      alloc_and_copy(str, buffer, len);)
}

extern "C" auto
sfcgal_io_read_binary_prepared(const char *str, size_t len)
    -> sfcgal_prepared_geometry_t *
{
  std::string const                         sstr(str, len);
  std::unique_ptr<SFCGAL::PreparedGeometry> g;

  try {
    g = SFCGAL::io::readBinaryPrepared(sstr);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During read_binary_prepared");
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return g.release();
}

extern "C" auto
sfcgal_io_read_ewkt(const char *str, size_t len) -> sfcgal_prepared_geometry_t *
{
  std::unique_ptr<SFCGAL::PreparedGeometry> g;

  try {
    g = SFCGAL::io::readEwkt(str, len);
  } catch (std::exception &e) {
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  SFCGAL::PreparedGeometry *pg = g.release();
  return pg;
}

/// @{
/// @privatesection
#ifndef DOXYGEN_SHOULD_SKIP_THIS

SFCGAL_GEOMETRY_FUNCTION_BINARY_PREDICATE(covers, SFCGAL::algorithm::covers)
SFCGAL_GEOMETRY_FUNCTION_BINARY_PREDICATE(covers_3d,
                                          SFCGAL::algorithm::covers3D)

SFCGAL_GEOMETRY_FUNCTION_BINARY_PREDICATE(intersects,
                                          SFCGAL::algorithm::intersects)
SFCGAL_GEOMETRY_FUNCTION_BINARY_PREDICATE(intersects_3d,
                                          SFCGAL::algorithm::intersects3D)

SFCGAL_GEOMETRY_FUNCTION_BINARY_MEASURE(distance, SFCGAL::algorithm::distance)
SFCGAL_GEOMETRY_FUNCTION_BINARY_MEASURE(distance_3d,
                                        SFCGAL::algorithm::distance3D)

SFCGAL_GEOMETRY_FUNCTION_BINARY_CONSTRUCTION(intersection,
                                             SFCGAL::algorithm::intersection)
SFCGAL_GEOMETRY_FUNCTION_BINARY_CONSTRUCTION(intersection_3d,
                                             SFCGAL::algorithm::intersection3D)
SFCGAL_GEOMETRY_FUNCTION_BINARY_CONSTRUCTION(difference,
                                             SFCGAL::algorithm::difference)
SFCGAL_GEOMETRY_FUNCTION_BINARY_CONSTRUCTION(difference_3d,
                                             SFCGAL::algorithm::difference3D)
SFCGAL_GEOMETRY_FUNCTION_BINARY_CONSTRUCTION(union, SFCGAL::algorithm::union_)
SFCGAL_GEOMETRY_FUNCTION_BINARY_CONSTRUCTION(union_3d,
                                             SFCGAL::algorithm::union3D)

SFCGAL_GEOMETRY_FUNCTION_UNARY_CONSTRUCTION(convexhull,
                                            SFCGAL::algorithm::convexHull)
SFCGAL_GEOMETRY_FUNCTION_UNARY_CONSTRUCTION(convexhull_3d,
                                            SFCGAL::algorithm::convexHull3D)
SFCGAL_GEOMETRY_FUNCTION_UNARY_CONSTRUCTION(straight_skeleton,
                                            SFCGAL::algorithm::straightSkeleton)
SFCGAL_GEOMETRY_FUNCTION_UNARY_CONSTRUCTION(
    approximate_medial_axis, SFCGAL::algorithm::approximateMedialAxis)
SFCGAL_GEOMETRY_FUNCTION_UNARY_CONSTRUCTION(tesselate,
                                            SFCGAL::algorithm::tesselate)

SFCGAL_GEOMETRY_FUNCTION_UNARY_MEASURE(area, SFCGAL::algorithm::area)
SFCGAL_GEOMETRY_FUNCTION_UNARY_MEASURE(area_3d, SFCGAL::algorithm::area3D)

#endif // ifndef DOXYGEN_SHOULD_SKIP_THIS
/// @} end of private section
/// @publicsection

extern "C" auto
sfcgal_geometry_volume(const sfcgal_geometry_t *geom) -> double
{
  double r = std::numeric_limits<double>::quiet_NaN();

  try {
    r = CGAL::to_double(SFCGAL::algorithm::volume(
        *static_cast<const SFCGAL::Geometry *>(geom)));
  } catch (std::exception &e) {
    SFCGAL_WARNING("During volume(A) :");
    SFCGAL_WARNING(
        "  with A: %s",
        static_cast<const SFCGAL::Geometry *>(geom)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return -1.0;
  }

  return r;
}

extern "C" auto
sfcgal_geometry_is_planar(const sfcgal_geometry_t *ga) -> int
{
  const auto *g = static_cast<const SFCGAL::Geometry *>(ga);

  if (g->geometryTypeId() != SFCGAL::TYPE_POLYGON) {
    SFCGAL_ERROR("is_planar() only applies to polygons");
    return -1;
  }

  bool r = false;

  try {
    r = SFCGAL::algorithm::isPlane3D<SFCGAL::Kernel>(
        g->as<const SFCGAL::Polygon>(), 1e-9);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During is_planar(A) :");
    SFCGAL_WARNING("  with A: %s",
                   static_cast<const SFCGAL::Geometry *>(ga)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return -1.0;
  }

  return r ? 1 : 0;
}

/*
 * Get geometry orientation.
 * Returns:
 * -1 for a counter clock wise orientation,
 * 1 for a clock wise orientation,
 * 0 for invalid or undetermined orientation
 */
extern "C" auto
sfcgal_geometry_orientation(const sfcgal_geometry_t *ga) -> int
{
  const auto *g = reinterpret_cast<const SFCGAL::Geometry *>(ga);

  if (g->geometryTypeId() != SFCGAL::TYPE_POLYGON) {
    SFCGAL_ERROR("orientation() only applies to polygons");
    return 0;
  }

  if (g->isEmpty()) {
    // Empty polygon has invalid/undetermined orientation
    return 0;
  }

  bool r = false;

  try {
    r = g->as<const SFCGAL::Polygon>().isCounterClockWiseOriented();
  } catch (std::exception &e) {
    SFCGAL_WARNING("During orientation(A) :");
    SFCGAL_WARNING("  with A: %s",
                   static_cast<const SFCGAL::Geometry *>(ga)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return -1;
  }

  return r ? -1 : 1;
}

extern "C" auto
sfcgal_geometry_make_solid(const sfcgal_geometry_t *ga) -> sfcgal_geometry_t *
{
  const auto *g = reinterpret_cast<const SFCGAL::Geometry *>(ga);

  if (g->geometryTypeId() != SFCGAL::TYPE_POLYHEDRALSURFACE) {
    SFCGAL_ERROR("make_solid() only applies to polyhedral surfaces");
    return nullptr;
  }

  return static_cast<SFCGAL::Geometry *>(
      new SFCGAL::Solid(g->as<const SFCGAL::PolyhedralSurface>()));
}

extern "C" auto
sfcgal_geometry_force_lhr(const sfcgal_geometry_t *ga) -> sfcgal_geometry_t *
{
  const auto *g = reinterpret_cast<const SFCGAL::Geometry *>(ga);
  std::unique_ptr<SFCGAL::Geometry>   gb = g->clone();
  SFCGAL::transform::ForceOrderPoints force(/* ccw */ true);

  try {
    gb->accept(force);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During force_lhr(A) :");
    SFCGAL_WARNING("  with A: %s",
                   static_cast<const SFCGAL::Geometry *>(ga)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return gb.release();
}

extern "C" auto
sfcgal_geometry_force_rhr(const sfcgal_geometry_t *ga) -> sfcgal_geometry_t *
{
  const auto *g = reinterpret_cast<const SFCGAL::Geometry *>(ga);
  std::unique_ptr<SFCGAL::Geometry>   gb = g->clone();
  SFCGAL::transform::ForceOrderPoints force(/* ccw */ false);

  try {
    gb->accept(force);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During force_rhr(A) :");
    SFCGAL_WARNING("  with A: %s",
                   static_cast<const SFCGAL::Geometry *>(ga)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return gb.release();
}

extern "C" auto
sfcgal_geometry_triangulate_2dz(const sfcgal_geometry_t *ga)
    -> sfcgal_geometry_t *
{
  const auto *g    = reinterpret_cast<const SFCGAL::Geometry *>(ga);
  auto       *surf = new SFCGAL::TriangulatedSurface;

  try {
    SFCGAL::triangulate::ConstraintDelaunayTriangulation cdt;
    SFCGAL::triangulate::triangulate2DZ(*g, cdt);
    cdt.getTriangles(*surf);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During triangulate_2d(A) :");
    SFCGAL_WARNING("  with A: %s",
                   static_cast<const SFCGAL::Geometry *>(ga)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return static_cast<SFCGAL::Geometry *>(surf);
}

extern "C" auto
sfcgal_geometry_extrude(const sfcgal_geometry_t *ga, double x, double y,
                        double z) -> sfcgal_geometry_t *
{
  const auto *g = reinterpret_cast<const SFCGAL::Geometry *>(ga);
  std::unique_ptr<SFCGAL::Geometry>    gb = g->clone();
  SFCGAL::transform::ForceZOrderPoints forceZ;
  std::unique_ptr<SFCGAL::Geometry>    result;

  try {
    gb->accept(forceZ);
    result = SFCGAL::algorithm::extrude(*gb, x, y, z);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During extrude(A, %g, %g, %g) :", x, y, z);
    SFCGAL_WARNING("  with A: %s",
                   static_cast<const SFCGAL::Geometry *>(ga)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return result.release();
}

extern "C" auto
sfcgal_geometry_round(const sfcgal_geometry_t *ga, int scale)
    -> sfcgal_geometry_t *
{
  const auto *g = reinterpret_cast<const SFCGAL::Geometry *>(ga);
  std::unique_ptr<SFCGAL::Geometry> gb = g->clone();
  //	SFCGAL_WARNING( "geom: %s %s", gb->asText().c_str(), typeid(g).name() );

  SFCGAL::transform::RoundTransform roundT(scale);

  try {
    gb->accept(roundT);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During round(A):");
    SFCGAL_WARNING("  with A: %s",
                   static_cast<const SFCGAL::Geometry *>(ga)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  //	SFCGAL_WARNING( "processed geom: %s", gb->asText().c_str() );
  return gb.release();
}

extern "C" auto
sfcgal_geometry_minkowski_sum(const sfcgal_geometry_t *ga,
                              const sfcgal_geometry_t *gb)
    -> sfcgal_geometry_t *
{
  const auto *g1 = reinterpret_cast<const SFCGAL::Geometry *>(ga);
  const auto *g2 = reinterpret_cast<const SFCGAL::Geometry *>(gb);

  if (g2->geometryTypeId() != SFCGAL::TYPE_POLYGON) {
    SFCGAL_ERROR("minkowski_sum(): the second argument must be a polygon");
    return nullptr;
  }

  std::unique_ptr<SFCGAL::Geometry> sum;

  try {
    sum = SFCGAL::algorithm::minkowskiSum(*g1, g2->as<const SFCGAL::Polygon>());
  } catch (std::exception &e) {
    SFCGAL_WARNING("During minkowski_sum(A,B):");
    SFCGAL_WARNING("  with A: %s",
                   static_cast<const SFCGAL::Geometry *>(ga)->asText().c_str());
    SFCGAL_WARNING("   and B: %s",
                   static_cast<const SFCGAL::Geometry *>(gb)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return sum.release();
}

extern "C" auto
sfcgal_geometry_offset_polygon(const sfcgal_geometry_t *ga, double offset)
    -> sfcgal_geometry_t *
{
  const auto *g1 = reinterpret_cast<const SFCGAL::Geometry *>(ga);
  std::unique_ptr<SFCGAL::MultiPolygon> mp;

  try {
    mp = SFCGAL::algorithm::offset(*g1, offset);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During offset(A,%g):", offset);
    SFCGAL_WARNING("  with A: %s",
                   static_cast<const SFCGAL::Geometry *>(ga)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return mp.release();
}

extern "C" void
sfcgal_geometry_force_valid(sfcgal_geometry_t *geom, int valid)
{
  auto *g1 = reinterpret_cast<SFCGAL::Geometry *>(geom);
  SFCGAL::algorithm::propagateValidityFlag(*g1, valid != 0);
}

extern "C" auto
sfcgal_geometry_has_validity_flag(const sfcgal_geometry_t *geom) -> int
{
  const auto *g1 = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  return g1->hasValidityFlag() ? 1 : 0;
}

extern "C" auto
sfcgal_geometry_extrude_straight_skeleton(const sfcgal_geometry_t *geom,
                                          double height) -> sfcgal_geometry_t *
{
  const auto *g1 = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  std::unique_ptr<SFCGAL::PolyhedralSurface> polys;

  try {
    polys = SFCGAL::algorithm::extrudeStraightSkeleton(*g1, height);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During straight_extrude_skeleton_distance(A):");
    SFCGAL_WARNING(
        "  with A: %s",
        static_cast<const SFCGAL::Geometry *>(geom)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return polys.release();
}

extern "C" auto
sfcgal_geometry_extrude_polygon_straight_skeleton(const sfcgal_geometry_t *geom,
                                                  double building_height,
                                                  double roof_height)
    -> sfcgal_geometry_t *
{
  const auto *g1 = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  std::unique_ptr<SFCGAL::Geometry> polys;

  try {
    polys = SFCGAL::algorithm::extrudeStraightSkeleton(*g1, building_height,
                                                       roof_height);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During straight_extrude_skeleton_distance(A):");
    SFCGAL_WARNING(
        "  with A: %s",
        static_cast<const SFCGAL::Geometry *>(geom)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return polys.release();
}

extern "C" auto
sfcgal_geometry_straight_skeleton_distance_in_m(const sfcgal_geometry_t *geom)
    -> sfcgal_geometry_t *
{
  const auto *g1 = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  std::unique_ptr<SFCGAL::MultiLineString> mls;

  try {
    mls = SFCGAL::algorithm::straightSkeleton(*g1, /*autoOrientation*/ true,
                                              /*innerOnly*/ false,
                                              /*outputDistanceInM*/ true);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During straight_skeleton_distance_in_m(A):");
    SFCGAL_WARNING(
        "  with A: %s",
        static_cast<const SFCGAL::Geometry *>(geom)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return mls.release();
}

extern "C" auto
sfcgal_geometry_line_sub_string(const sfcgal_geometry_t *geom, double start,
                                double end) -> sfcgal_geometry_t *
{
  const auto *g1 = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  if (g1->geometryTypeId() != SFCGAL::TYPE_LINESTRING) {
    SFCGAL_ERROR("line_sub_string(): the first argument must be a lineString");
    return nullptr;
  }
  std::unique_ptr<SFCGAL::LineString> ls;
  try {
    ls = SFCGAL::algorithm::lineSubstring(g1->as<const SFCGAL::LineString>(),
                                          start, end);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During line_sub_string(A, %g, %g):", start, end);
    SFCGAL_WARNING(
        "  with A: %s",
        static_cast<const SFCGAL::Geometry *>(geom)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return ls.release();
}

#if !_MSC_VER
extern "C" auto
sfcgal_geometry_alpha_shapes(const sfcgal_geometry_t *geom, double alpha,
                             bool allow_holes) -> sfcgal_geometry_t *
{
  const auto *g1 = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  std::unique_ptr<SFCGAL::Geometry> result;

  try {
    result = SFCGAL::algorithm::alphaShapes(g1->as<const SFCGAL::Geometry>(),
                                            alpha, allow_holes);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During alphaShapes(A,%g):", alpha);
    SFCGAL_WARNING(
        "  with A: %s",
        static_cast<const SFCGAL::Geometry *>(geom)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return result.release();
}

extern "C" auto
sfcgal_geometry_optimal_alpha_shapes(const sfcgal_geometry_t *geom,
                                     bool allow_holes, size_t nb_components)
    -> sfcgal_geometry_t *
{
  const auto *g1 = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  std::unique_ptr<SFCGAL::Geometry> result;

  try {
    result = SFCGAL::algorithm::optimal_alpha_shapes(
        g1->as<const SFCGAL::Geometry>(), allow_holes, nb_components);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During optimal_alpha_shapes(A, %g %g):",
                   static_cast<int>(allow_holes), nb_components);
    SFCGAL_WARNING(
        "  with A: %s",
        static_cast<const SFCGAL::Geometry *>(geom)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return result.release();
}
#endif

extern "C" auto
sfcgal_geometry_alpha_wrapping_3d(const sfcgal_geometry_t *geom,
                                  size_t relativeAlpha, size_t relativeOffset)
    -> sfcgal_geometry_t *
{
  const auto *inputGeom = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  std::unique_ptr<SFCGAL::Geometry> result;

  try {
    result = SFCGAL::algorithm::alphaWrapping3D(
        inputGeom->as<const SFCGAL::Geometry>(), relativeAlpha, relativeOffset);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During alpha_wrapping_3d(A, %g %g):", relativeAlpha,
                   relativeOffset);
    SFCGAL_WARNING(
        "  with A: %s",
        static_cast<const SFCGAL::Geometry *>(geom)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return result.release();
}

extern "C" auto
sfcgal_y_monotone_partition_2(const sfcgal_geometry_t *geom)
    -> sfcgal_geometry_t *
{
  const auto *g1 = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  std::unique_ptr<SFCGAL::Geometry> result;

  try {
    result = SFCGAL::algorithm::partition_2(g1->as<const SFCGAL::Geometry>(),
                                            SFCGAL::algorithm::y_monotone);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During y_monotone_partition_2(A):");
    SFCGAL_WARNING(
        "  with A: %s",
        static_cast<const SFCGAL::Geometry *>(geom)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return result.release();
}

extern "C" auto
sfcgal_approx_convex_partition_2(const sfcgal_geometry_t *geom)
    -> sfcgal_geometry_t *
{
  const auto *g1 = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  std::unique_ptr<SFCGAL::Geometry> result;

  try {
    result = SFCGAL::algorithm::partition_2(g1->as<const SFCGAL::Geometry>(),
                                            SFCGAL::algorithm::approx_convex);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During approx_convex_partition_2(A):");
    SFCGAL_WARNING(
        "  with A: %s",
        static_cast<const SFCGAL::Geometry *>(geom)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return result.release();
}

extern "C" auto
sfcgal_greene_approx_convex_partition_2(const sfcgal_geometry_t *geom)
    -> sfcgal_geometry_t *
{
  const auto *g1 = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  std::unique_ptr<SFCGAL::Geometry> result;

  try {
    result =
        SFCGAL::algorithm::partition_2(g1->as<const SFCGAL::Geometry>(),
                                       SFCGAL::algorithm::greene_approx_convex);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During greene_approx_convex_partition_2(A):");
    SFCGAL_WARNING(
        "  with A: %s",
        static_cast<const SFCGAL::Geometry *>(geom)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return result.release();
}
extern "C" auto
sfcgal_optimal_convex_partition_2(const sfcgal_geometry_t *geom)
    -> sfcgal_geometry_t *
{
  const auto *g1 = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  std::unique_ptr<SFCGAL::Geometry> result;

  try {
    result = SFCGAL::algorithm::partition_2(g1->as<const SFCGAL::Geometry>(),
                                            SFCGAL::algorithm::optimal_convex);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During optimal_convex_partition_2(A):");
    SFCGAL_WARNING(
        "  with A: %s",
        static_cast<const SFCGAL::Geometry *>(geom)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return result.release();
}

extern "C" auto
sfcgal_geometry_visibility_point(const sfcgal_geometry_t *polygon,
                                 const sfcgal_geometry_t *point)
    -> sfcgal_geometry_t *
{

  const auto *poly = reinterpret_cast<const SFCGAL::Geometry *>(polygon);
  const auto *pt   = reinterpret_cast<const SFCGAL::Geometry *>(point);
  std::unique_ptr<SFCGAL::Geometry> result;

  if (poly->geometryTypeId() != SFCGAL::TYPE_POLYGON) {
    SFCGAL_ERROR("visibility() only applies to polygons");
    return result.release();
  }

  if (pt->geometryTypeId() != SFCGAL::TYPE_POINT) {
    SFCGAL_ERROR("second argument must be a point");
    return result.release();
  }

  try {
    result = SFCGAL::algorithm::visibility(poly->as<const SFCGAL::Polygon>(),
                                           pt->as<const SFCGAL::Point>());
  } catch (std::exception &e) {
    SFCGAL_WARNING("During visibility(A, B) :");
    SFCGAL_WARNING("  with A: %s", poly->asText().c_str());
    SFCGAL_WARNING("  and B: %s", pt->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return result.release();
  }

  return result.release();
}

extern "C" auto
sfcgal_geometry_visibility_segment(const sfcgal_geometry_t *polygon,
                                   const sfcgal_geometry_t *pointA,
                                   const sfcgal_geometry_t *pointB)
    -> sfcgal_geometry_t *
{
  const auto *poly = reinterpret_cast<const SFCGAL::Geometry *>(polygon);
  const auto *ptA  = reinterpret_cast<const SFCGAL::Geometry *>(pointA);
  const auto *ptB  = reinterpret_cast<const SFCGAL::Geometry *>(pointB);
  std::unique_ptr<SFCGAL::Geometry> result;

  if (poly->geometryTypeId() != SFCGAL::TYPE_POLYGON) {
    SFCGAL_ERROR("visibility() only applies to polygons");
    return result.release();
  }

  if ((ptA->geometryTypeId() != SFCGAL::TYPE_POINT) ||
      (ptB->geometryTypeId() != SFCGAL::TYPE_POINT)) {
    SFCGAL_ERROR("second and third argument must be a point");
    return result.release();
  }

  try {
    result = SFCGAL::algorithm::visibility(poly->as<const SFCGAL::Polygon>(),
                                           ptA->as<const SFCGAL::Point>(),
                                           ptB->as<const SFCGAL::Point>());
  } catch (std::exception &e) {
    SFCGAL_WARNING("During visibility(A, B, C) :");
    SFCGAL_WARNING("  with A: %s", poly->asText().c_str());
    SFCGAL_WARNING("  and B: %s", ptA->asText().c_str());
    SFCGAL_WARNING("  and C: %s", ptB->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return result.release();
  }

  return result.release();
}

extern "C" auto

sfcgal_geometry_translate_2d(const sfcgal_geometry_t *geom, double dx,
                             double dy) -> sfcgal_geometry_t *
{
  const auto *g = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  std::unique_ptr<SFCGAL::Geometry> gb = g->clone();
  try {

    SFCGAL::algorithm::translate(*gb, SFCGAL::Kernel::Vector_2(dx, dy));
  } catch (std::exception &e) {
    SFCGAL_WARNING("During translate(A, %g, %g):", dx, dy);
    SFCGAL_WARNING("  with A: %s", gb->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return gb.release();
}

extern "C" auto
sfcgal_geometry_translate_3d(const sfcgal_geometry_t *geom, double dx,
                             double dy, double dz) -> sfcgal_geometry_t *
{
  const auto *g = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  std::unique_ptr<SFCGAL::Geometry> gb = g->clone();

  try {
    SFCGAL::algorithm::translate(*gb, SFCGAL::Kernel::Vector_3(dx, dy, dz));
  } catch (std::exception &e) {
    SFCGAL_WARNING("During translate(A, %g, %g, %g):", dx, dy, dz);
    SFCGAL_WARNING("  with A: %s", gb->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return gb.release();
}

extern "C" auto
sfcgal_geometry_transform(const sfcgal_geometry_t *geom, const float *matrix4x4)
    -> sfcgal_geometry_t *
{
  const auto *g = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  std::unique_ptr<SFCGAL::Geometry> gb = g->clone();

  // row by row
  CGAL::Aff_transformation_3<SFCGAL::Kernel> transf(
      matrix4x4[(0 * 4) + 0], matrix4x4[(1 * 4) + 0], matrix4x4[(2 * 4) + 0],
      matrix4x4[(3 * 4) + 0], //
      matrix4x4[(0 * 4) + 1], matrix4x4[(1 * 4) + 1], matrix4x4[(2 * 4) + 1],
      matrix4x4[(3 * 4) + 1], //
      matrix4x4[(0 * 4) + 2], matrix4x4[(1 * 4) + 2], matrix4x4[(2 * 4) + 2],
      matrix4x4[(3 * 4) + 2], //
      matrix4x4[(3 * 4) + 3]);

  SFCGAL::transform::AffineTransform3 visitor(transf);
  gb->accept(visitor);

  return gb.release();
}

extern "C" auto
sfcgal_geometry_scale(const sfcgal_geometry_t *geom, double s)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const SFCGAL::Geometry &g =
          *reinterpret_cast<const SFCGAL::Geometry *>(geom);
      std::unique_ptr<SFCGAL::Geometry> result = g.clone();
      SFCGAL::algorithm::scale(*result, s); return result.release();)
}

extern "C" auto
sfcgal_geometry_scale_3d(const sfcgal_geometry_t *geom, double sx, double sy,
                         double sz) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const SFCGAL::Geometry &g =
          *reinterpret_cast<const SFCGAL::Geometry *>(geom);
      std::unique_ptr<SFCGAL::Geometry> result = g.clone();
      SFCGAL::algorithm::scale(*result, sx, sy, sz); return result.release();)
}

extern "C" auto
sfcgal_geometry_scale_3d_around_center(const sfcgal_geometry_t *geom, double sx,
                                       double sy, double sz, double cx,
                                       double cy, double cz)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const SFCGAL::Geometry &g =
          *reinterpret_cast<const SFCGAL::Geometry *>(geom);
      std::unique_ptr<SFCGAL::Geometry> result = g.clone();
      SFCGAL::algorithm::scale(*result, sx, sy, sz, cx, cy, cz);
      return result.release();)
}

extern "C" auto
sfcgal_geometry_rotate(const sfcgal_geometry_t *geom, double angle)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const SFCGAL::Geometry &g =
          *reinterpret_cast<const SFCGAL::Geometry *>(geom);
      std::unique_ptr<SFCGAL::Geometry> result = g.clone();
      SFCGAL::algorithm::rotate(*result, angle); return result.release();)
}

extern "C" auto
sfcgal_geometry_rotate_2d(const sfcgal_geometry_t *geom, double angle,
                          double cx, double cy) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const SFCGAL::Geometry &g =
          *reinterpret_cast<const SFCGAL::Geometry *>(geom);
      std::unique_ptr<SFCGAL::Geometry> result = g.clone();
      SFCGAL::algorithm::rotate(*result, angle, SFCGAL::Point(cx, cy));
      return result.release();)
}

extern "C" auto
sfcgal_geometry_rotate_3d(const sfcgal_geometry_t *geom, double angle,
                          double ax, double ay, double az)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const SFCGAL::Geometry &g =
          *reinterpret_cast<const SFCGAL::Geometry *>(geom);
      std::unique_ptr<SFCGAL::Geometry> result = g.clone();
      SFCGAL::algorithm::rotate(*result, angle,
                                SFCGAL::Kernel::Vector_3(ax, ay, az));
      return result.release();)
}

extern "C" auto
sfcgal_geometry_rotate_3d_around_center(const sfcgal_geometry_t *geom,
                                        double angle, double ax, double ay,
                                        double az, double cx, double cy,
                                        double cz) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const SFCGAL::Geometry &g =
          *reinterpret_cast<const SFCGAL::Geometry *>(geom);
      std::unique_ptr<SFCGAL::Geometry> result = g.clone();
      SFCGAL::algorithm::rotate(*result, angle,
                                SFCGAL::Kernel::Vector_3(ax, ay, az),
                                SFCGAL::Point(cx, cy, cz));
      return result.release();)
}

extern "C" auto
sfcgal_geometry_rotate_x(const sfcgal_geometry_t *geom, double angle)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const SFCGAL::Geometry &g =
          *reinterpret_cast<const SFCGAL::Geometry *>(geom);
      std::unique_ptr<SFCGAL::Geometry> result = g.clone();
      SFCGAL::algorithm::rotateX(*result, angle); return result.release();)
}

extern "C" auto
sfcgal_geometry_rotate_y(const sfcgal_geometry_t *geom, double angle)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const SFCGAL::Geometry &g =
          *reinterpret_cast<const SFCGAL::Geometry *>(geom);
      std::unique_ptr<SFCGAL::Geometry> result = g.clone();
      SFCGAL::algorithm::rotateY(*result, angle); return result.release();)
}

extern "C" auto
sfcgal_geometry_rotate_z(const sfcgal_geometry_t *geom, double angle)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const SFCGAL::Geometry &g =
          *reinterpret_cast<const SFCGAL::Geometry *>(geom);
      std::unique_ptr<SFCGAL::Geometry> result = g.clone();
      SFCGAL::algorithm::rotateZ(*result, angle); return result.release();)
}

extern "C" auto
sfcgal_geometry_straight_skeleton_partition(const sfcgal_geometry_t *geom,
                                            bool autoOrientation)
    -> sfcgal_geometry_t *
{
  std::unique_ptr<SFCGAL::Geometry> result;
  try {
    result = SFCGAL::algorithm::straightSkeletonPartition(
        *(const SFCGAL::Geometry *)(geom));
  } catch (std::exception &e) {
    SFCGAL_WARNING("During straight_skeleton_partition (A, %g) :",
                   static_cast<int>(autoOrientation));
    SFCGAL_WARNING(
        "  with A: %s",
        static_cast<const SFCGAL::Geometry *>(geom)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }
  return result.release();
}

extern "C" auto
sfcgal_geometry_buffer3d(const sfcgal_geometry_t *geom, double radius,
                         int segments, sfcgal_buffer3d_type_t buffer_type)
    -> sfcgal_geometry_t *
{
  std::unique_ptr<SFCGAL::Geometry> result;
  try {
    SFCGAL::algorithm::Buffer3D::BufferType type;

    switch (buffer_type) {
    case SFCGAL_BUFFER3D_ROUND:
      type = SFCGAL::algorithm::Buffer3D::ROUND;
      break;
    case SFCGAL_BUFFER3D_CYLSPHERE:
      type = SFCGAL::algorithm::Buffer3D::CYLSPHERE;
      break;
    case SFCGAL_BUFFER3D_FLAT:
      type = SFCGAL::algorithm::Buffer3D::FLAT;
      break;
    default:
      SFCGAL_ERROR("Invalid buffer type");
      return nullptr;
    }

    SFCGAL::algorithm::Buffer3D buffer3d(*(const SFCGAL::Geometry *)(geom),
                                         radius, segments);
    result = buffer3d.compute(type);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During buffer3d (A, %g, %d, %d) :", radius, segments,
                   buffer_type);
    SFCGAL_WARNING(
        "  with A: %s",
        static_cast<const SFCGAL::Geometry *>(geom)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }
  return result.release();
}

extern "C" auto
sfcgal_geometry_envelope(const sfcgal_geometry_t *geom) -> sfcgal_geometry_t *
{
  const auto      *geometry = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  SFCGAL::Envelope result;

  try {
    result = geometry->envelope();
  } catch (std::exception &e) {
    SFCGAL_WARNING("During envelope(A):");
    SFCGAL_WARNING("  with A: %s", geometry->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return result.toPolygon().release();
}

extern "C" auto
sfcgal_geometry_envelope_3d(const sfcgal_geometry_t *geom)
    -> sfcgal_geometry_t *
{
  const auto      *geometry = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  SFCGAL::Envelope result;

  try {
    result = geometry->envelope();
  } catch (std::exception &e) {
    SFCGAL_WARNING("During envelope_3d(A):");
    SFCGAL_WARNING("  with A: %s", geometry->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  if (result.is3D()) {
    return result.toShell().release();
  }

  return result.toPolygon().release();
}

extern "C" auto
sfcgal_geometry_length(const sfcgal_geometry_t *geom) -> double
{
  const auto *geometry = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  double      result;

  try {
    result = SFCGAL::algorithm::length(*geometry);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During length(A):");
    SFCGAL_WARNING("  with A: %s", geometry->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return std::numeric_limits<double>::quiet_NaN();
  }

  return result;
}

extern "C" auto
sfcgal_geometry_length_3d(const sfcgal_geometry_t *geom) -> double
{
  const auto *geometry = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  double      result;

  try {
    result = SFCGAL::algorithm::length3D(*geometry);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During length(A):");
    SFCGAL_WARNING("  with A: %s", geometry->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return std::numeric_limits<double>::quiet_NaN();
  }

  return result;
}

SFCGAL_API auto
sfcgal_geometry_boundary(const sfcgal_geometry_t *geom) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return reinterpret_cast<const SFCGAL::Geometry *>(geom)
          ->boundary()
          .release();)
}

extern "C" auto
sfcgal_geometry_centroid(const sfcgal_geometry_t *geom) -> sfcgal_geometry_t *
{
  const auto *geometry = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  std::unique_ptr<SFCGAL::Point> result;

  try {
    result = SFCGAL::algorithm::centroid(*geometry);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During centroid(A):");
    SFCGAL_WARNING("  with A: %s", geometry->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  std::unique_ptr<SFCGAL::Geometry> out = result->clone();
  return out.release();
}

extern "C" auto
sfcgal_geometry_centroid_3d(const sfcgal_geometry_t *geom)
    -> sfcgal_geometry_t *
{
  const auto *geometry = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  std::unique_ptr<SFCGAL::Point> result;

  try {
    result = SFCGAL::algorithm::centroid3D(*geometry);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During centroid3D(A):");
    SFCGAL_WARNING("  with A: %s", geometry->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  std::unique_ptr<SFCGAL::Geometry> out = result->clone();
  return out.release();
}

extern "C" auto
sfcgal_geometry_is_equals(const sfcgal_geometry_t *ga,
                          const sfcgal_geometry_t *gb) -> int
{
  return sfcgal_geometry_is_almost_equals(ga, gb, -1.0);
}

extern "C" auto
sfcgal_geometry_is_almost_equals(const sfcgal_geometry_t *ga,
                                 const sfcgal_geometry_t *gb, double tolerance)
    -> int
{
  const auto *g1 = reinterpret_cast<const SFCGAL::Geometry *>(ga);
  const auto *g2 = reinterpret_cast<const SFCGAL::Geometry *>(gb);

  bool result;
  try {
    result = g1->almostEqual(*g2, tolerance);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During is_almost_equals(A, B, %g):", tolerance);
    SFCGAL_WARNING("  with A: %s", g1->asText().c_str());
    SFCGAL_WARNING("  with B: %s", g2->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    result = false;
  }

  return static_cast<int>(result);
}

extern "C" auto
sfcgal_geometry_simplify(const sfcgal_geometry_t *geom, double threshold,
                         bool preserveTopology) -> sfcgal_geometry_t *
{
  const auto *geometry = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  std::unique_ptr<SFCGAL::Geometry> result;

  try {
    result =
        SFCGAL::algorithm::simplify(*geometry, threshold, preserveTopology);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During simplify(A, %g, %d):", threshold,
                   static_cast<int>(preserveTopology));
    SFCGAL_WARNING("  with A: %s", geometry->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  std::unique_ptr<SFCGAL::Geometry> out = result->clone();
  return out.release();
}

extern "C" auto
sfcgal_primitive_create(sfcgal_primitive_type_t primitive_type)
    -> sfcgal_primitive_t *
{
  sfcgal_primitive_t *primitive = nullptr;

  switch (primitive_type) {
  case SFCGAL_TYPE_CYLINDER:
    primitive = new SFCGAL::Cylinder();
    break;
  case SFCGAL_TYPE_SPHERE:
    primitive = new SFCGAL::Sphere();
    break;
  case SFCGAL_TYPE_TORUS:
    primitive = new SFCGAL::Torus();
    break;
  case SFCGAL_TYPE_BOX:
    primitive = new SFCGAL::Box();
    break;
  case SFCGAL_TYPE_CUBE:
    primitive = new SFCGAL::Cube();
    break;
  case SFCGAL_TYPE_CONE:
    primitive = new SFCGAL::Cone();
    break;
  }

  if (primitive == nullptr) {
    SFCGAL_ERROR("Invalid primitive type");
  }

  return primitive;
}

extern "C" auto
sfcgal_primitive_delete(sfcgal_primitive_t *primitive) -> void
{
  delete reinterpret_cast<SFCGAL::Primitive *>(primitive);
}

extern "C" auto
sfcgal_primitive_clone(const sfcgal_primitive_t *prim) -> sfcgal_primitive_t *
{
  const auto *primitive1Cast =
      reinterpret_cast<const SFCGAL::Primitive *>(prim);
  sfcgal_primitive_t *out = sfcgal_primitive_create(
      static_cast<sfcgal_primitive_type_t>(primitive1Cast->primitiveTypeId()));

  auto *outCast = reinterpret_cast<SFCGAL::Primitive *>(out);
  outCast->operator=(*primitive1Cast);

  return out;
}

extern "C" auto
sfcgal_primitive_is_almost_equals(const sfcgal_primitive_t *prim1,
                                  const sfcgal_primitive_t *prim2,
                                  double                    tolerance) -> int
{
  const auto *prim1Cast = reinterpret_cast<const SFCGAL::Primitive *>(prim1);
  const auto *prim2Cast = reinterpret_cast<const SFCGAL::Primitive *>(prim2);

  bool result;
  try {
    result = prim1Cast->almostEqual(*prim2Cast, tolerance);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During primitive_is_almost_equals(A, B, %g):", tolerance);
    SFCGAL_WARNING(
        "  with A: %s",
        (prim1Cast == nullptr ? "null" : prim1Cast->toString().c_str()));
    SFCGAL_WARNING(
        "  with B: %s",
        (prim2Cast == nullptr ? "null" : prim2Cast->toString().c_str()));
    SFCGAL_ERROR("%s", e.what());
    result = false;
  }

  return static_cast<int>(result);
}

extern "C" auto
sfcgal_primitive_area(const sfcgal_primitive_t *prim, bool withDiscretization)
    -> double
{
  const auto *primCast = reinterpret_cast<const SFCGAL::Primitive *>(prim);

  double result;
  try {
    result = primCast->area3D(withDiscretization);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During primitive_area(A, %d):",
                   static_cast<int>(withDiscretization));
    SFCGAL_WARNING(
        "  with A: %s",
        (primCast == nullptr ? "null" : primCast->toString().c_str()));
    SFCGAL_ERROR("%s", e.what());
    result = -1.0;
  }

  return result;
}

extern "C" auto
sfcgal_primitive_volume(const sfcgal_primitive_t *prim, bool withDiscretization)
    -> double
{
  const auto *primCast = reinterpret_cast<const SFCGAL::Primitive *>(prim);

  double result;
  try {
    result = primCast->volume(withDiscretization);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During primitive_volume(A, %d):",
                   static_cast<int>(withDiscretization));
    SFCGAL_WARNING(
        "  with A: %s",
        (primCast == nullptr ? "null" : primCast->toString().c_str()));
    SFCGAL_ERROR("%s", e.what());
    result = -1.0;
  }

  return result;
}

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection
#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace SFCGAL {

// Temp struct dedicated to JSON serialize/deserialize
struct PrimitiveParameterDesc {
  std::string        name;
  std::string        type;
  PrimitiveParameter value;
};

// deserialize json to PrimitiveParameterDesc
void
from_json(const nlohmann::json &jsonObject, PrimitiveParameterDesc &desc)
{
  jsonObject.at("name").get_to(desc.name);
  jsonObject.at("type").get_to(desc.type);
  if (jsonObject.contains("value")) {
    const nlohmann::json &value = jsonObject.at("value");
    if (desc.type == "int") {
      desc.value = value.get<unsigned int>();
    } else if (desc.type == "double") {
      desc.value = value.get<double>();
    } else if (desc.type == "point3") {
      std::vector<double> vect;
      vect = value.get<std::vector<double>>();
      Kernel::Point_3 point(
          vect[0], vect[1], //
          (vect.size() > 2 ? vect[2]
                           : std::numeric_limits<double>::quiet_NaN()), //
          (vect.size() > 3 ? vect[3]
                           : std::numeric_limits<double>::quiet_NaN()));
      desc.value = point;
    } else if (desc.type == "vector3") {
      std::vector<double> vect;
      vect = value.get<std::vector<double>>();
      Kernel::Vector_3 point(vect[0], vect[1], //
                             (vect.size() > 2
                                  ? vect[2]
                                  : std::numeric_limits<double>::quiet_NaN()));
      desc.value = point;
    } else {
      throw nlohmann::json::type_error::create(
          306, (boost::format("Unknown type '%1%'.") % desc.type).str(),
          nullptr);
    }
  }
}

// serialize PrimitiveParameterDesc to json
void
to_json(nlohmann::json &jsonObject, const PrimitiveParameterDesc &desc)
{
  jsonObject["name"] = desc.name;
  if (std::holds_alternative<SFCGAL::Kernel::FT>(desc.value)) {
    jsonObject["type"] = "double";
  } else if (std::holds_alternative<unsigned int>(desc.value)) {
    jsonObject["type"] = "int";
  } else if (std::holds_alternative<SFCGAL::Kernel::Point_3>(desc.value)) {
    jsonObject["type"] = "point3";
  } else if (std::holds_alternative<SFCGAL::Kernel::Vector_3>(desc.value)) {
    jsonObject["type"] = "vector3";
  } else {
    jsonObject["type"] = "unknown";
  }

  if (std::holds_alternative<SFCGAL::Kernel::FT>(desc.value)) {
    jsonObject["value"] =
        CGAL::to_double(std::get<SFCGAL::Kernel::FT>(desc.value));
  } else if (std::holds_alternative<unsigned int>(desc.value)) {
    jsonObject["value"] = std::get<unsigned int>(desc.value);
  } else if (std::holds_alternative<SFCGAL::Kernel::Point_3>(desc.value)) {
    SFCGAL::Kernel::Point_3 point =
        std::get<SFCGAL::Kernel::Point_3>(desc.value);
    jsonObject["value"] = std::vector<double>{CGAL::to_double(point.x()),
                                              CGAL::to_double(point.y()),
                                              CGAL::to_double(point.z())};
  } else if (std::holds_alternative<SFCGAL::Kernel::Vector_3>(desc.value)) {
    SFCGAL::Kernel::Vector_3 point =
        std::get<SFCGAL::Kernel::Vector_3>(desc.value);
    jsonObject["value"] = std::vector<double>{CGAL::to_double(point.x()),
                                              CGAL::to_double(point.y()),
                                              CGAL::to_double(point.z())};
  }
}

} // namespace SFCGAL
#endif // ifndef DOXYGEN_SHOULD_SKIP_THIS
/// @} end of private section
/// @publicsection

extern "C" auto
sfcgal_primitive_parameters(const sfcgal_primitive_t *primitive, char **buffer,
                            size_t *len) -> void
{
  auto do_get_parameters = [](const sfcgal_primitive_t *primitive,
                              char **buffer, size_t *len) {
    const auto *primitiveCast =
        reinterpret_cast<const SFCGAL::Primitive *>(primitive);

    std::vector<nlohmann::json> params;
    for (const auto &[key, value] : primitiveCast->parameters()) {
      SFCGAL::PrimitiveParameterDesc desc;
      desc.name  = key;
      desc.value = value;
      nlohmann::json param(desc);
      param.erase("value"); // remove useless value in this case
      params.emplace_back(param);
    }

    nlohmann::json jsonOut(params);
    std::string    paramsStr = jsonOut.dump();
    alloc_and_copy(paramsStr, buffer, len);
  };

  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      do_get_parameters(primitive, buffer, len);)
}

extern "C" auto
sfcgal_primitive_parameter(const sfcgal_primitive_t *primitive,
                           const char *name, char **buffer, size_t *len) -> void
{
  auto do_get_parameter = [](const sfcgal_primitive_t *primitive,
                             const char *name, char **buffer, size_t *len) {
    const auto *primitiveCast =
        reinterpret_cast<const SFCGAL::Primitive *>(primitive);

    const SFCGAL::PrimitiveParameter value = primitiveCast->parameter(name);
    SFCGAL::PrimitiveParameterDesc   desc;
    desc.name  = std::string(name);
    desc.value = value;

    nlohmann::json param(desc);
    std::string    paramsStr = param.dump();
    alloc_and_copy(paramsStr, buffer, len);
  };

  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      do_get_parameter(primitive, name, buffer, len);)
}

extern "C" auto
sfcgal_primitive_set_parameter(sfcgal_primitive_t *primitive, const char *name,
                               const char *jsonChars) -> void
{
  auto do_set_parameter = [](sfcgal_primitive_t *primitive, const char *name,
                             const char *jsonChars) {
    try {
      auto *primitiveCast = reinterpret_cast<SFCGAL::Primitive *>(primitive);

      nlohmann::json                 jsonObj = nlohmann::json::parse(jsonChars);
      SFCGAL::PrimitiveParameterDesc newValue =
          jsonObj.template get<SFCGAL::PrimitiveParameterDesc>();

      primitiveCast->setParameter(name, newValue.value);
    } catch (std::bad_alloc const &e) {
      SFCGAL_ERROR(e.what());
    } catch (...) {
      SFCGAL_ERROR("Json exception");
    }
  };

  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      do_set_parameter(primitive, name, jsonChars);)
}

extern "C" auto
sfcgal_primitive_parameter_double(const sfcgal_primitive_t *primitive,
                                  const char               *name) -> double
{
  double returnValue = 0.;

  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const auto *primitiveCast =
          reinterpret_cast<const SFCGAL::Primitive *>(primitive);
      SFCGAL::PrimitiveParameter parameter =
          primitiveCast->parameter(std::string(name));
      const SFCGAL::Kernel::FT *kernel_parameter =
          std::get_if<SFCGAL::Kernel::FT>(&parameter);
      if (kernel_parameter != nullptr) {
        returnValue = CGAL::to_double(*kernel_parameter);
      } else {
        SFCGAL_ERROR("Parameter %s is not a double", name);
      } return returnValue;);
}

extern "C" auto
sfcgal_primitive_set_parameter_double(sfcgal_primitive_t *primitive,
                                      const char *name, double parameter)
    -> void
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      auto *primitiveCast = reinterpret_cast<SFCGAL::Primitive *>(primitive);
      primitiveCast->setParameter(std::string(name),
                                  SFCGAL::Kernel::FT(parameter));)
}

extern "C" auto
sfcgal_primitive_parameter_int(const sfcgal_primitive_t *primitive,
                               const char               *name) -> unsigned int
{
  unsigned int returnValue = 0;

  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const auto *primitiveCast =
          reinterpret_cast<const SFCGAL::Primitive *>(primitive);
      SFCGAL::PrimitiveParameter parameter =
          primitiveCast->parameter(std::string(name));
      const unsigned int *kernel_parameter =
          std::get_if<unsigned int>(&parameter);
      if (kernel_parameter != nullptr) {
        returnValue = *kernel_parameter;
      } else {
        SFCGAL_ERROR("Parameter %s is not a unsigned int", name);
      } return returnValue;);
}

extern "C" auto
sfcgal_primitive_set_parameter_int(sfcgal_primitive_t *primitive,
                                   const char *name, unsigned int parameter)
    -> void
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      auto *primitiveCast = reinterpret_cast<SFCGAL::Primitive *>(primitive);
      primitiveCast->setParameter(std::string(name), parameter);)
}

/// @{
/// @privatesection

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template <class T>
auto
primitive_parameter_array(const sfcgal_primitive_t *primitive, const char *name,
                          const char *parameter_type) -> double *
{
  double *returnArray = nullptr;

  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const auto *primitiveCast =
          reinterpret_cast<const SFCGAL::Primitive *>(primitive);
      SFCGAL::PrimitiveParameter parameter =
          primitiveCast->parameter(std::string(name));

      std::cout << "primitive_parameter_array name:" << name << " ==> "
                << *primitiveCast << "\n";

      const T *array = std::get_if<T>(&parameter); //
      if (array != nullptr) {
        returnArray =
            static_cast<double *>(sfcgal_alloc_handler(3 * sizeof(double)));
        returnArray[0] = CGAL::to_double(array->x());
        returnArray[1] = CGAL::to_double(array->y());
        returnArray[2] = CGAL::to_double(array->z());
      } else {
        SFCGAL_ERROR("Parameter %s is not a %s", name, parameter_type);
      } //
      return returnArray; //
  );
}

#endif // ifndef DOXYGEN_SHOULD_SKIP_THIS
/// @} end of private section
/// @publicsection

template <class T>
auto
primitive_set_parameter_array(sfcgal_primitive_t *primitive, const char *name,
                              const double *array) -> void
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      auto *primitiveCast = reinterpret_cast<SFCGAL::Primitive *>(primitive);
      primitiveCast->setParameter(std::string(name),
                                  T(array[0], array[1], array[2]));)
}

extern "C" auto
sfcgal_primitive_parameter_point(const sfcgal_primitive_t *primitive,
                                 const char               *name) -> double *
{
  return primitive_parameter_array<SFCGAL::Kernel::Point_3>(primitive, name,
                                                            "point");
}

extern "C" auto
sfcgal_primitive_set_parameter_point(sfcgal_primitive_t *primitive,
                                     const char *name, const double *point)
    -> void
{
  primitive_set_parameter_array<SFCGAL::Kernel::Point_3>(primitive, name,
                                                         point);
}

extern "C" auto
sfcgal_primitive_parameter_vector(const sfcgal_primitive_t *primitive,
                                  const char               *name) -> double *
{
  return primitive_parameter_array<SFCGAL::Kernel::Vector_3>(primitive, name,
                                                             "vector");
}

extern "C" auto
sfcgal_primitive_set_parameter_vector(sfcgal_primitive_t *primitive,
                                      const char *name, const double *vector)
    -> void
{
  primitive_set_parameter_array<SFCGAL::Kernel::Vector_3>(primitive, name,
                                                          vector);
}

extern "C" auto
sfcgal_primitive_as_polyhedral_surface(const sfcgal_primitive_t *primitive)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      auto *primitiveCast =
          reinterpret_cast<const SFCGAL::Primitive *>(primitive);
      return new SFCGAL::PolyhedralSurface(
          primitiveCast->generatePolyhedralSurface());)
}

/*--------------------------------------------------------------------------------------*
 *
 * NURBSCurve support
 *
 *--------------------------------------------------------------------------------------*/

// Helper functions for enum conversion
namespace {
auto
cKnotMethodToCpp(sfcgal_knot_method_t method) -> SFCGAL::NURBSCurve::KnotMethod
{
  switch (method) {
  case SFCGAL_KNOT_METHOD_UNIFORM:
    return SFCGAL::NURBSCurve::KnotMethod::UNIFORM;
  case SFCGAL_KNOT_METHOD_CHORD_LENGTH:
    return SFCGAL::NURBSCurve::KnotMethod::CHORD_LENGTH;
  case SFCGAL_KNOT_METHOD_CENTRIPETAL:
  default:
    return SFCGAL::NURBSCurve::KnotMethod::CENTRIPETAL;
  }
}

auto
cEndConditionToCpp(sfcgal_end_condition_t condition)
    -> SFCGAL::NURBSCurve::EndCondition
{
  switch (condition) {
  case SFCGAL_END_CONDITION_CLAMPED:
    return SFCGAL::NURBSCurve::EndCondition::CLAMPED;
  case SFCGAL_END_CONDITION_NATURAL:
    return SFCGAL::NURBSCurve::EndCondition::NATURAL;
  case SFCGAL_END_CONDITION_PERIODIC:
    return SFCGAL::NURBSCurve::EndCondition::PERIODIC;
  case SFCGAL_END_CONDITION_TANGENT:
    return SFCGAL::NURBSCurve::EndCondition::TANGENT;
  default:
    return SFCGAL::NURBSCurve::EndCondition::CLAMPED;
  }
}

auto
cFitMethodToCpp(sfcgal_fit_method_t method) -> SFCGAL::NURBSCurve::FitMethod
{
  switch (method) {
  case SFCGAL_FIT_METHOD_INTERPOLATE:
    return SFCGAL::NURBSCurve::FitMethod::INTERPOLATE;
  case SFCGAL_FIT_METHOD_APPROXIMATE:
    return SFCGAL::NURBSCurve::FitMethod::APPROXIMATE;
  default:
    return SFCGAL::NURBSCurve::FitMethod::INTERPOLATE;
  }
}

// Helper function for safe geometry to Point conversion
auto
safeGeometryToPoint(const sfcgal_geometry_t *geom) -> const SFCGAL::Point *
{
  if (geom == nullptr) {
    return nullptr;
  }
  const auto *geometry = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  if (geometry->geometryTypeId() != SFCGAL::TYPE_POINT) {
    BOOST_THROW_EXCEPTION(SFCGAL::Exception("Expected Point geometry"));
  }
  return static_cast<const SFCGAL::Point *>(geometry);
}

// Helper function for safe geometry to NURBSCurve conversion
auto
safeGeometryToNURBSCurve(sfcgal_geometry_t *geom) -> SFCGAL::NURBSCurve *
{
  if (geom == nullptr) {
    BOOST_THROW_EXCEPTION(SFCGAL::Exception("Geometry cannot be null"));
  }
  auto *geometry = reinterpret_cast<SFCGAL::Geometry *>(geom);
  if (geometry->geometryTypeId() != SFCGAL::TYPE_NURBSCURVE) {
    BOOST_THROW_EXCEPTION(SFCGAL::Exception("Expected NURBSCurve geometry"));
  }
  return static_cast<SFCGAL::NURBSCurve *>(geometry);
}

// Helper function for safe geometry to const NURBSCurve conversion
auto
safeGeometryToConstNURBSCurve(const sfcgal_geometry_t *geom)
    -> const SFCGAL::NURBSCurve *
{
  if (geom == nullptr) {
    BOOST_THROW_EXCEPTION(SFCGAL::Exception("Geometry cannot be null"));
  }
  const auto *geometry = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  if (geometry->geometryTypeId() != SFCGAL::TYPE_NURBSCURVE) {
    BOOST_THROW_EXCEPTION(SFCGAL::Exception("Expected NURBSCurve geometry"));
  }
  return static_cast<const SFCGAL::NURBSCurve *>(geometry);
}

auto
convertPointArray(const sfcgal_geometry_t **points, size_t num_points)
    -> std::vector<SFCGAL::Point>
{
  if (points == nullptr) {
    BOOST_THROW_EXCEPTION(SFCGAL::Exception("Points array cannot be null"));
  }
  std::vector<SFCGAL::Point> result;
  result.reserve(num_points);
  for (size_t i = 0; i < num_points; ++i) {
    const auto *point = safeGeometryToPoint(points[i]);
    if (point == nullptr) {
      BOOST_THROW_EXCEPTION(SFCGAL::Exception("Point cannot be null"));
    }
    result.push_back(*point);
  }
  return result;
}

auto
convertWeightArray(const double *weights, size_t num_weights)
    -> std::vector<SFCGAL::NURBSCurve::FT>
{
  std::vector<SFCGAL::NURBSCurve::FT> result;
  if (weights != nullptr) {
    result.reserve(num_weights);
    for (size_t i = 0; i < num_weights; ++i) {
      result.emplace_back(weights[i]);
    }
  }
  return result;
}

auto
convertKnotArray(const double *knots, size_t num_knots)
    -> std::vector<SFCGAL::NURBSCurve::Knot>
{
  std::vector<SFCGAL::NURBSCurve::Knot> result;
  result.reserve(num_knots);
  for (size_t i = 0; i < num_knots; ++i) {
    result.emplace_back(knots[i]);
  }
  return result;
}

} // anonymous namespace

extern "C" auto
sfcgal_nurbs_curve_create(void) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      return static_cast<SFCGAL::Geometry *>(new SFCGAL::NURBSCurve());)
}

extern "C" auto
sfcgal_nurbs_curve_create_from_points(const sfcgal_geometry_t **points,
                                      size_t num_points, unsigned int degree,
                                      sfcgal_knot_method_t knot_method)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      auto pointVec = convertPointArray(points, num_points);
      auto curve    = std::make_unique<SFCGAL::NURBSCurve>(
          pointVec, degree, cKnotMethodToCpp(knot_method));
      return static_cast<SFCGAL::Geometry *>(curve.release());)
}

extern "C" auto
sfcgal_nurbs_curve_create_from_points_and_weights(
    const sfcgal_geometry_t **points, const double *weights, size_t num_points,
    unsigned int degree, sfcgal_knot_method_t knot_method)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      auto pointVec  = convertPointArray(points, num_points);
      auto weightVec = convertWeightArray(weights, num_points);
      auto curve     = std::make_unique<SFCGAL::NURBSCurve>(
          pointVec, weightVec, degree, cKnotMethodToCpp(knot_method));
      return static_cast<SFCGAL::Geometry *>(curve.release());)
}

extern "C" auto
sfcgal_nurbs_curve_create_from_full_data(const sfcgal_geometry_t **points,
                                         const double             *weights,
                                         size_t num_points, unsigned int degree,
                                         const double *knots, size_t num_knots)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      auto pointVec  = convertPointArray(points, num_points);
      auto weightVec = convertWeightArray(weights, num_points);
      auto knotVec   = convertKnotArray(knots, num_knots);
      auto curve     = std::make_unique<SFCGAL::NURBSCurve>(pointVec, weightVec,
                                                            degree, knotVec);
      return static_cast<SFCGAL::Geometry *>(curve.release());)
}

extern "C" auto
sfcgal_nurbs_curve_create_bezier(const sfcgal_geometry_t **points,
                                 size_t num_points) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      auto pointVec = convertPointArray(points, num_points);
      auto curve    = SFCGAL::NURBSCurve::fromBezier(pointVec);
      return static_cast<SFCGAL::Geometry *>(curve.release());)
}

extern "C" auto
sfcgal_nurbs_curve_create_bspline(const sfcgal_geometry_t **points,
                                  size_t num_points, unsigned int degree)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      auto pointVec = convertPointArray(points, num_points);
      auto curve    = SFCGAL::NURBSCurve::createBSpline(pointVec, degree);
      return static_cast<SFCGAL::Geometry *>(curve.release());)
}

extern "C" auto
sfcgal_nurbs_curve_create_circular_arc(const sfcgal_geometry_t *center,
                                       double radius, double start_angle,
                                       double                   end_angle,
                                       const sfcgal_geometry_t *normal)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const auto *centerPoint = safeGeometryToPoint(center);
      if (centerPoint == nullptr) {
        BOOST_THROW_EXCEPTION(SFCGAL::Exception("Center point cannot be null"));
      } const auto *normalPoint =
          normal ? safeGeometryToPoint(normal) : nullptr;

      SFCGAL::Point normalVec(0, 0, 1);
      if (normalPoint != nullptr) { normalVec = *normalPoint; }

      auto curve = SFCGAL::NURBSCurve::createCircularArc(
          *centerPoint, SFCGAL::NURBSCurve::FT(radius),
          SFCGAL::NURBSCurve::FT(start_angle),
          SFCGAL::NURBSCurve::FT(end_angle), normalVec);
      return static_cast<SFCGAL::Geometry *>(curve.release());)
}

extern "C" auto
sfcgal_nurbs_curve_interpolate(const sfcgal_geometry_t **points,
                               size_t num_points, unsigned int degree,
                               sfcgal_knot_method_t   knot_method,
                               sfcgal_end_condition_t end_condition)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      auto pointVec = convertPointArray(points, num_points);
      auto curve    = SFCGAL::NURBSCurve::interpolateCurve(
          pointVec, degree, cKnotMethodToCpp(knot_method),
          cEndConditionToCpp(end_condition));
      return static_cast<SFCGAL::Geometry *>(curve.release());)
}

extern "C" auto
sfcgal_nurbs_curve_approximate(const sfcgal_geometry_t **points,
                               size_t num_points, unsigned int degree,
                               double tolerance, size_t max_control_points)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      auto pointVec = convertPointArray(points, num_points);
      auto curve    = SFCGAL::NURBSCurve::approximateCurve(
          pointVec, degree, SFCGAL::NURBSCurve::FT(tolerance),
          max_control_points);
      return static_cast<SFCGAL::Geometry *>(curve.release());)
}

extern "C" auto
sfcgal_nurbs_curve_fit(const sfcgal_geometry_t **points, size_t num_points,
                       unsigned int degree, sfcgal_fit_method_t fit_method,
                       sfcgal_knot_method_t   knot_method,
                       sfcgal_end_condition_t end_condition, double tolerance,
                       size_t max_control_points) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      auto pointVec = convertPointArray(points, num_points);
      auto curve    = SFCGAL::NURBSCurve::fitCurve(
          pointVec, degree, cFitMethodToCpp(fit_method),
          cKnotMethodToCpp(knot_method), cEndConditionToCpp(end_condition),
          SFCGAL::NURBSCurve::FT(tolerance), max_control_points);
      return static_cast<SFCGAL::Geometry *>(curve.release());)
}

extern "C" auto
sfcgal_nurbs_curve_num_control_points(const sfcgal_geometry_t *curve) -> size_t
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const auto *nurbsCurve = safeGeometryToConstNURBSCurve(curve);
      return nurbsCurve->numControlPoints();)
}

extern "C" auto
sfcgal_nurbs_curve_control_point_n(const sfcgal_geometry_t *curve, size_t index)
    -> const sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const auto *nurbsCurve   = safeGeometryToConstNURBSCurve(curve);
      const auto &controlPoint = nurbsCurve->controlPointN(index);
      return static_cast<const SFCGAL::Geometry *>(&controlPoint);)
}

extern "C" auto
sfcgal_nurbs_curve_set_control_point_n(sfcgal_geometry_t *curve, size_t index,
                                       const sfcgal_geometry_t *point) -> void
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      auto       *nurbsCurve = safeGeometryToNURBSCurve(curve);
      const auto *pointPtr   = safeGeometryToPoint(point);
      if (pointPtr == nullptr) {
        BOOST_THROW_EXCEPTION(SFCGAL::Exception("Point cannot be null"));
      } nurbsCurve->setControlPoint(index, *pointPtr);)
}

extern "C" auto
sfcgal_nurbs_curve_weight_n(const sfcgal_geometry_t *curve, size_t index)
    -> double
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const auto *nurbsCurve = safeGeometryToConstNURBSCurve(curve);
      auto weight = nurbsCurve->weight(index); return CGAL::to_double(weight);)
}

extern "C" auto
sfcgal_nurbs_curve_set_weight_n(sfcgal_geometry_t *curve, size_t index,
                                double weight) -> void
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR_NO_RET(
      auto *nurbsCurve = safeGeometryToNURBSCurve(curve);
      nurbsCurve->setWeight(index, SFCGAL::NURBSCurve::FT(weight));)
}

extern "C" auto
sfcgal_nurbs_curve_degree(const sfcgal_geometry_t *curve) -> unsigned int
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const auto *nurbsCurve = safeGeometryToConstNURBSCurve(curve);
      return nurbsCurve->degree();)
}

extern "C" auto
sfcgal_nurbs_curve_num_knots(const sfcgal_geometry_t *curve) -> size_t
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const auto *nurbsCurve = safeGeometryToConstNURBSCurve(curve);
      return nurbsCurve->numKnots();)
}

extern "C" auto
sfcgal_nurbs_curve_knot_n(const sfcgal_geometry_t *curve, size_t index)
    -> double
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const auto *nurbsCurve = safeGeometryToConstNURBSCurve(curve);
      auto        knot = nurbsCurve->knot(index); return CGAL::to_double(knot);)
}

extern "C" auto
sfcgal_nurbs_curve_is_rational(const sfcgal_geometry_t *curve) -> int
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const auto *nurbsCurve = safeGeometryToConstNURBSCurve(curve);
      return nurbsCurve->isRational() ? 1 : 0;)
}

extern "C" auto
sfcgal_nurbs_curve_is_bezier(const sfcgal_geometry_t *curve) -> int
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const auto *nurbsCurve = safeGeometryToConstNURBSCurve(curve);
      return nurbsCurve->isBezier() ? 1 : 0;)
}

extern "C" auto
sfcgal_nurbs_curve_is_bspline(const sfcgal_geometry_t *curve) -> int
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const auto *nurbsCurve = safeGeometryToConstNURBSCurve(curve);
      return nurbsCurve->isBSpline() ? 1 : 0;)
}

extern "C" auto
sfcgal_nurbs_curve_evaluate(const sfcgal_geometry_t *curve, double parameter)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const auto *nurbsCurve = safeGeometryToConstNURBSCurve(curve);
      auto        point =
          nurbsCurve->evaluate(SFCGAL::NURBSCurve::Parameter(parameter));
      return static_cast<SFCGAL::Geometry *>(new SFCGAL::Point(point));)
}

extern "C" auto
sfcgal_nurbs_curve_derivative(const sfcgal_geometry_t *curve, double parameter,
                              unsigned int order) -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const auto *nurbsCurve = safeGeometryToConstNURBSCurve(curve);
      auto        derivative = nurbsCurve->derivative(
          SFCGAL::NURBSCurve::Parameter(parameter), order);
      return static_cast<SFCGAL::Geometry *>(new SFCGAL::Point(derivative));)
}

extern "C" auto
sfcgal_nurbs_curve_to_linestring(const sfcgal_geometry_t *curve,
                                 unsigned int             num_segments)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const auto *nurbsCurve = safeGeometryToConstNURBSCurve(curve);
      auto        lineString = nurbsCurve->toLineString(num_segments);
      return static_cast<SFCGAL::Geometry *>(lineString.release());)
}

extern "C" auto
sfcgal_nurbs_curve_to_linestring_adaptive(const sfcgal_geometry_t *curve,
                                          double                   tolerance,
                                          unsigned int             min_segments,
                                          unsigned int             max_segments)
    -> sfcgal_geometry_t *
{
  SFCGAL_GEOMETRY_CONVERT_CATCH_TO_ERROR(
      const auto *nurbsCurve = safeGeometryToConstNURBSCurve(curve);
      auto        lineString = nurbsCurve->toLineStringAdaptive(
          SFCGAL::NURBSCurve::FT(tolerance), min_segments, max_segments);
      return static_cast<SFCGAL::Geometry *>(lineString.release());)
}

#if SFCGAL_CGAL_VERSION_MAJOR >= 6
extern "C" auto
sfcgal_geometry_polygon_repair(const sfcgal_geometry_t     *geom,
                               sfcgal_polygon_repair_rule_t repair_rule)
    -> sfcgal_geometry_t *
{
  std::unique_ptr<SFCGAL::Geometry> result;
  try {
    const auto *geometry = reinterpret_cast<const SFCGAL::Geometry *>(geom);

    // Validate input geometry type
    if (geometry->geometryTypeId() != SFCGAL::TYPE_POLYGON &&
        geometry->geometryTypeId() != SFCGAL::TYPE_MULTIPOLYGON) {
      SFCGAL_ERROR("polygon_repair requires Polygon or MultiPolygon geometry");
      return nullptr;
    }

    SFCGAL::algorithm::PolygonRepairRule rule;

    switch (repair_rule) {
    case SFCGAL_POLYGON_REPAIR_EVEN_ODD:
      rule = SFCGAL::algorithm::PolygonRepairRule::EVEN_ODD_RULE;
      break;

    case SFCGAL_POLYGON_REPAIR_NON_ZERO:
  #if SFCGAL_CGAL_VERSION_MAJOR == 6 && SFCGAL_CGAL_VERSION_MINOR >= 1
      rule = SFCGAL::algorithm::PolygonRepairRule::NON_ZERO_RULE;
  #else
      SFCGAL_ERROR("NON_ZERO rule requires CGAL 6.1 or later");
      return nullptr;
  #endif
      break;

    case SFCGAL_POLYGON_REPAIR_UNION:
  #if SFCGAL_CGAL_VERSION_MAJOR == 6 && SFCGAL_CGAL_VERSION_MINOR >= 1
      rule = SFCGAL::algorithm::PolygonRepairRule::UNION_RULE;
  #else
      SFCGAL_ERROR("UNION rule requires CGAL 6.1 or later");
      return nullptr;
  #endif
      break;

    case SFCGAL_POLYGON_REPAIR_INTERSECTION:
  #if SFCGAL_CGAL_VERSION_MAJOR == 6 && SFCGAL_CGAL_VERSION_MINOR >= 1
      rule = SFCGAL::algorithm::PolygonRepairRule::INTERSECTION_RULE;
  #else
      SFCGAL_ERROR("INTERSECTION rule requires CGAL 6.1 or later");
      return nullptr;
  #endif
      break;

    default:
      SFCGAL_ERROR("Invalid polygon repair rule");
      return nullptr;
    }

    result = SFCGAL::algorithm::polygonRepair(*geometry, rule);

  } catch (std::exception &e) {
    SFCGAL_WARNING("During polygon_repair (A, %d):", repair_rule);
    SFCGAL_WARNING("  with A: %s",
                   ((const SFCGAL::Geometry *)(geom))->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return reinterpret_cast<sfcgal_geometry_t *>(result.release());
}
#endif // SFCGAL_CGAL_VERSION_MAJOR >= 6

extern "C" auto
sfcgal_geometry_projected_medial_axis(const sfcgal_geometry_t *geom)
    -> sfcgal_geometry_t *
{
  const auto *geometry = reinterpret_cast<const SFCGAL::Geometry *>(geom);
  std::unique_ptr<SFCGAL::MultiLineString> multiLineString;

  try {
    multiLineString = SFCGAL::algorithm::approximateMedialAxis(*geometry, true);
  } catch (std::exception &e) {
    SFCGAL_WARNING("During projected_medial_axis(A):");
    SFCGAL_WARNING(
        "  with A: %s",
        static_cast<const SFCGAL::Geometry *>(geom)->asText().c_str());
    SFCGAL_ERROR("%s", e.what());
    return nullptr;
  }

  return multiLineString.release();
}
