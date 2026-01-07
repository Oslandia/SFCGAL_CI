// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "operations.hpp"

#include "../constructors.hpp"
#include "SFCGAL/config.h"
#include "SFCGAL/version.h"
#include <SFCGAL/Kernel.h>

#ifndef _MSC_VER
  #include <SFCGAL/algorithm/alphaShapes.h>
#endif
#include <SFCGAL/algorithm/alphaWrapping3D.h>
#include <SFCGAL/algorithm/area.h>
#include <SFCGAL/algorithm/buffer3D.h>
#include <SFCGAL/algorithm/centroid.h>
#include <SFCGAL/algorithm/collect.h>
#include <SFCGAL/algorithm/collectionExtract.h>
#include <SFCGAL/algorithm/collectionHomogenize.h>
#include <SFCGAL/algorithm/collectionToMulti.h>
#include <SFCGAL/algorithm/connection.h>
#include <SFCGAL/algorithm/convexHull.h>
#include <SFCGAL/algorithm/covers.h>
#include <SFCGAL/algorithm/difference.h>
#include <SFCGAL/algorithm/distance.h>
#include <SFCGAL/algorithm/distance3d.h>
#include <SFCGAL/algorithm/extrude.h>
#include <SFCGAL/algorithm/force2D.h>
#include <SFCGAL/algorithm/force3D.h>
#include <SFCGAL/algorithm/forceMeasured.h>
#include <SFCGAL/algorithm/intersection.h>
#include <SFCGAL/algorithm/intersects.h>
#include <SFCGAL/algorithm/isClosed.h>
#include <SFCGAL/algorithm/isSimple.h>
#include <SFCGAL/algorithm/isValid.h>
#include <SFCGAL/algorithm/length.h>
#include <SFCGAL/algorithm/lineSubstring.h>
#include <SFCGAL/algorithm/minkowskiSum.h>
#include <SFCGAL/algorithm/minkowskiSum3D.h>
#include <SFCGAL/algorithm/normal.h>
#include <SFCGAL/algorithm/offset.h>
#include <SFCGAL/algorithm/orientation.h>
#include <SFCGAL/algorithm/partition_2.h>
#include <SFCGAL/algorithm/plane.h>
#if SFCGAL_CGAL_VERSION_MAJOR >= 6
  #include <SFCGAL/algorithm/polygonRepair.h>
#endif
#include <SFCGAL/algorithm/roofGeneration.h>
#include <SFCGAL/algorithm/rotate.h>
#include <SFCGAL/algorithm/scale.h>
#include <SFCGAL/algorithm/simplification.h>
#include <SFCGAL/algorithm/straightSkeleton.h>
#include <SFCGAL/algorithm/tesselate.h>
#include <SFCGAL/algorithm/translate.h>
#include <SFCGAL/algorithm/union.h>
#include <SFCGAL/algorithm/visibility.h>
#include <SFCGAL/algorithm/volume.h>

#include <SFCGAL/Envelope.h>
#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/MultiSolid.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>

#include <CGAL/number_utils.h>
#include <cmath>

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

/**
 * @brief Structure defining a geometry operation for sfcgalop CLI
 *
 * Contains all metadata and implementation details for a single geometric
 * operation that can be executed via the command line interface.
 */
struct Operation {
  std::string name; ///< Operation name (e.g., "area", "intersection")
  std::string
      category; ///< Category for grouping (e.g., "Metrics", "Set Operations")
  std::string description; ///< Short description for operation list
  bool        requires_b;  ///< Whether operation requires a second geometry
  std::string param_help;  ///< Detailed help text with parameters and examples
  std::string input;       ///< Input specification (A, A,B, A,params)
  std::string output; ///< Output type (G=Geometry, D=Double, B=Boolean, T=Text)
  /// Function implementing the operation
  std::function<std::optional<OperationResult>(
      const std::string &, const SFCGAL::Geometry *, const SFCGAL::Geometry *)>
      func;
};

namespace {

/**
 * @brief Parse a string to double, returning a fallback on failure.
 *
 * Attempts to convert the given string to a double using std::stod.
 * If conversion fails (invalid format, out-of-range, etc.), the provided
 * default_val is returned.
 *
 * @param str Input string to parse.
 * @param default_val Value to return if parsing fails (default: 0.0).
 * @return double Parsed double on success, otherwise default_val.
 */
auto
parse_double(const std::string &str, double default_val = 0.0) -> double
{
  try {
    return std::stod(str);
  } catch (...) {
    return default_val;
  }
}

/**
 * @brief Parse a string as a boolean value.
 *
 * Accepts various boolean representations in a case-insensitive manner:
 * - "true", "t", "1" -> true
 * - "false", "f", "0" -> false
 *
 * @param str String to parse as boolean
 * @param default_val Default value to return if parsing fails
 * @return Boolean value or default_val if parsing fails
 */
auto
parse_boolean(const std::string &str, bool default_val = false) -> bool
{
  std::string lower_str = str;
  std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(),
                 [](unsigned char chr) -> int { return std::tolower(chr); });

  if (lower_str == "true" || lower_str == "t" || lower_str == "1") {
    return true;
  }
  if (lower_str == "false" || lower_str == "f" || lower_str == "0") {
    return false;
  }
  return default_val;
}

// Forward declaration for trim function
auto
trim(const std::string &str) -> std::string;

/**
 * @brief Extract and parse a boolean parameter from the parameter map.
 *
 * This function handles both string-based boolean parameters and legacy
 * double-based parameters (where 0.0 = false, non-zero = true).
 *
 * @param params Parameter map from parse_params
 * @param key Parameter name to extract
 * @param param_str Original parameter string for string-based parsing
 * @param default_val Default value if parameter not found
 * @return Boolean value
 */

auto
parse_boolean_param(const std::map<std::string, double> &params,
                    const std::string &key, const std::string &param_str,
                    bool default_val = false) -> bool
{
  auto it = params.find(key);
  if (it == params.end()) {
    return default_val;
  }

  // Try to extract the original string value from param_str for proper boolean
  // parsing
  std::string search_key = key + "=";
  auto        key_pos    = param_str.find(search_key);
  if (key_pos != std::string::npos) {
    auto value_start = key_pos + search_key.length();
    auto value_end   = param_str.find(',', value_start);
    if (value_end == std::string::npos) {
      value_end = param_str.length();
    }
    std::string value_str =
        param_str.substr(value_start, value_end - value_start);
    value_str = trim(value_str);

    // If it looks like a string boolean, parse it as such
    if (value_str == "true" || value_str == "false" || value_str == "t" ||
        value_str == "f" || value_str == "T" || value_str == "F" ||
        value_str == "True" || value_str == "False" || value_str == "TRUE" ||
        value_str == "FALSE") {
      return parse_boolean(value_str, default_val);
    }
  }

  // Fallback to standard if it's not 0(.0) it's true
  return it->second != 0.0;
}

/**
 * @brief Parse a comma-separated list of key=value pairs into a map.
 *
 * Parses `str` for entries of the form `key=value` separated by commas and
 * returns a map from each key to its value parsed as a double. Entries that
 * do not contain an '=' character are ignored. Values are converted using
 * `parse_double`; any parse fallback behavior is handled by that function.
 *
 * @param str Input string containing comma-separated `key=value` pairs.
 * @return std::map<std::string,double> Mapping of keys to parsed double values.
 */
/**
 * Trim leading and trailing whitespace from a string
 */
auto
trim(const std::string &str) -> std::string
{
  auto start = str.begin();
  auto end   = str.end();

  // Find first non-whitespace
  while (start != end &&
         (std::isspace(static_cast<unsigned char>(*start)) != 0)) {
    ++start;
  }

  // Find last non-whitespace (from the end)
  while (start != end &&
         (std::isspace(static_cast<unsigned char>(*(end - 1))) != 0)) {
    --end;
  }

  return {start, end};
}

auto
parse_params(const std::string &str) -> std::map<std::string, double>
{
  std::map<std::string, double> params;
  std::istringstream            stream(str);
  std::string                   param;

  while (std::getline(stream, param, ',')) {
    auto eq_pos = param.find('=');
    if (eq_pos != std::string::npos) {
      auto name_raw  = param.substr(0, eq_pos);
      auto value_raw = param.substr(eq_pos + 1);

      // Trim whitespace from both key and value to prevent lookup failures
      auto name_trimmed  = trim(name_raw);
      auto value_trimmed = trim(value_raw);

      params[name_trimmed] = parse_double(value_trimmed);
    }
  }
  return params;
}

// NOLINTNEXTLINE(cert-err58-cpp)
const std::vector<Operation> operations = {
    // Metrics
    {"area", "Metrics", "Calculate the 2D area of a geometry", false,
     "No parameters required.\n\nExample:\n  sfcgalop -a \"POLYGON((0 0,3 0,3 "
     "4,0 4,0 0))\" area",
     "A", "D",
     [](const std::string &, const SFCGAL::Geometry *geom,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return SFCGAL::algorithm::area(*geom);
     }},

    {"area3d", "Metrics", "Calculate the 3D surface area of a geometry", false,
     "No parameters required.\n\nExample:\n  sfcgalop -a \"POLYGON Z((0 0 0,3 "
     "0 0,3 4 2,0 4 2,0 0 0))\" area3d",
     "A", "D",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return SFCGAL::algorithm::area3D(*geom_a);
     }},

    {"volume", "Metrics", "Calculate the 3D volume of a solid geometry", false,
     "No parameters required.\nOnly works with solid geometries (SOLID, "
     "POLYHEDRALSURFACE).\n\nExample:\n  sfcgalop -a \"SOLID(...)\" volume",
     "A", "D",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return CGAL::to_double(SFCGAL::algorithm::volume(*geom_a));
     }},

    {"length", "Metrics", "Calculate the 2D length of linear geometries", false,
     "No parameters required.\n\nExample:\n  sfcgalop -a \"LINESTRING(0 0,3 "
     "4,6 0)\" length",
     "A", "D",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return SFCGAL::algorithm::length(*geom_a);
     }},

    {"length3d", "Metrics", "Calculate the 3D length of linear geometries",
     false,
     "No parameters required.\n\nExample:\n  sfcgalop -a \"LINESTRING Z(0 0 "
     "0,3 4 2,6 0 1)\" length3d",
     "A", "D",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return SFCGAL::algorithm::length3D(*geom_a);
     }},

    {"distance", "Metrics",
     "Calculate the 2D minimum distance between two geometries", true,
     "No parameters required.\nRequires two geometries.\n\nExample:\n  "
     "sfcgalop -a \"POINT(0 0)\" -b \"POINT(3 4)\" distance",
     "A, B", "D",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::distance(*geom_a, *geom_b);
     }},

    {"distance3d", "Metrics",
     "Calculate the 3D minimum distance between two geometries", true,
     "No parameters required.\nRequires two geometries.\n\nExample:\n  "
     "sfcgalop -a \"POINT Z(0 0 0)\" -b \"POINT Z(3 4 5)\" distance3d",
     "A, B", "D",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::distance3D(*geom_a, *geom_b);
     }},

    // Predicates
    {"intersects", "Predicates", "Test if two geometries intersect in 2D", true,
     "", "A, B", "B",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::intersects(*geom_a, *geom_b);
     }},

    {"intersects3d", "Predicates", "Test if two geometries intersect in 3D",
     true, "", "A, B", "B",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::intersects3D(*geom_a, *geom_b);
     }},

    {"covers", "Predicates", "Test if geometry A completely covers geometry B",
     true, "", "A, B", "B",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::covers(*geom_a, *geom_b);
     }},

    {"is_valid", "Predicates", "Test if geometry is topologically valid", false,
     "Parameters:\n  tolerance=VALUE: Absolute tolerance for planarity checks "
     "(default: 1e-6)\n\nExample:\n  sfcgalop -a \"POLYGON Z ((0 0 0,1 0 0,1 1 "
     "0.0001,0 0 0))\" is_valid \"tolerance=1e-3\"",
     "A, params", "B",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double tolerance =
           params.count("tolerance") ? params["tolerance"] : 1e-6;
       return static_cast<bool>(SFCGAL::algorithm::isValid(*geom_a, tolerance));
     }},

    {"is_simple", "Predicates", "Test if geometry has no self-intersections",
     false, "", "A", "B",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return static_cast<bool>(SFCGAL::algorithm::isSimple(*geom_a));
     }},

    {"is_closed", "Predicates", "Test if linear geometry forms a closed ring",
     false, "", "A", "B",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return static_cast<bool>(SFCGAL::algorithm::isClosed(*geom_a));
     }},

    {"is_3d", "Predicates", "Test if geometry has Z coordinates", false, "",
     "A", "B",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return geom_a->is3D();
     }},

    {"is_measured", "Predicates",
     "Test if geometry has measure (M) coordinates", false, "", "A", "B",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return geom_a->isMeasured();
     }},

    {"is_empty", "Predicates", "Test if geometry contains no points", false, "",
     "A", "B",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return geom_a->isEmpty();
     }},

    // Set operations
    {"intersection", "Set Operations",
     "Compute the geometric intersection of two geometries", true, "", "A, B",
     "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::intersection(*geom_a, *geom_b);
     }},

    {"intersection3d", "Set Operations",
     "Compute the 3D geometric intersection of two geometries", true, "",
     "A, B", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::intersection3D(*geom_a, *geom_b);
     }},

    {"difference", "Set Operations", "Compute geometry A minus geometry B",
     true, "", "A, B", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::difference(*geom_a, *geom_b);
     }},

    {"difference3d", "Set Operations", "Compute 3D geometry A minus geometry B",
     true, "", "A, B", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::difference3D(*geom_a, *geom_b);
     }},

    {"union", "Set Operations", "Compute the geometric union of two geometries",
     true, "", "A, B", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::union_(*geom_a, *geom_b);
     }},

    {"union3d", "Set Operations",
     "Compute the 3D geometric union of two geometries", true, "", "A, B", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::union3D(*geom_a, *geom_b);
     }},

    // Construction
    {"boundary", "Construction",
     "Compute the topological boundary of a geometry", false, "", "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return geom_a->boundary();
     }},

    {"envelope", "Construction", "Compute the minimum bounding rectangle",
     false, "", "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       SFCGAL::Envelope env = geom_a->envelope();
       return env.toPolygon();
     }},

    {"convexhull", "Construction", "Compute the 2D convex hull of a geometry",
     false, "", "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return SFCGAL::algorithm::convexHull(*geom_a);
     }},

    {"convexhull3d", "Construction", "Compute the 3D convex hull of a geometry",
     false, "", "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return SFCGAL::algorithm::convexHull3D(*geom_a);
     }},

    {"centroid", "Construction", "Compute the geometric centroid of a geometry",
     false, "", "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return SFCGAL::algorithm::centroid(*geom_a);
     }},

    {"straightskeleton", "Construction",
     "Compute the straight skeleton of a polygon", false,
     "Parameters:\n  auto_orientation=BOOL: Enable automatic orientation "
     "correction (default: false)\n"
     "                         Accepts: true/false, t/f, 1/0, TRUE/FALSE "
     "(case-insensitive)\n\n"
     "Example:\n  sfcgalop -a \"POLYGON((0 0,4 0,4 "
     "4,0 4,0 0))\" straightskeleton \"auto_orientation=true\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto params = parse_params(args);
       bool autoOrientation =
           parse_boolean_param(params, "auto_orientation", args, false);
       return SFCGAL::algorithm::straightSkeleton(*geom_a, autoOrientation);
     }},

    {"medial_axis", "Construction",
     "Compute the approximate medial axis of a polygon", false,
     "Computes the approximate medial axis for a polygon as a "
     "MultiLineString.\n"
     "The medial axis represents the 'skeleton' of the polygon, where each "
     "point\n"
     "is equidistant from the nearest polygon edges.\n\n"
     "Input A: Polygon geometry\n\n"
     "Parameters:\n"
     "  project_to_edges=BOOL: If true, extend free endpoints to polygon "
     "boundary\n"
     "                         using edge midpoint method (default: false)\n\n"
     "Examples:\n"
     "  sfcgalop -a \"POLYGON((0 0,6 0,6 3,0 3,0 0))\" medial_axis\n"
     "  sfcgalop -a \"POLYGON((0 0,6 0,6 3,0 3,0 0))\" medial_axis "
     "\"project_to_edges=true\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto params = parse_params(args);
       bool projectToEdges =
           parse_boolean_param(params, "project_to_edges", args, false);
       return SFCGAL::algorithm::approximateMedialAxis(*geom_a, projectToEdges);
     }},

    {"extrude", "Construction", "Extrude a 2D geometry to create a 3D solid",
     false,
     "Parameters:\n  dx=X: X-axis extrusion distance\n  dy=Y: Y-axis extrusion "
     "distance\n  dz=Z: Z-axis extrusion distance (default: 1.0)\n\nExample:\n "
     " sfcgalop -a \"POLYGON((0 0,1 0,1 1,0 1,0 0))\" extrude "
     "\"dx=0,dy=0,dz=2\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double dx     = params["dx"];
       double dy     = params["dy"];
       double dz     = params.count("dz") ? params["dz"] : 1.0;
       return SFCGAL::algorithm::extrude(*geom_a, dx, dy, dz);
     }},

    {"extrudeUntil", "Construction",
     "Extrude a 2D geometry until another "
     "geometry (typically roof) to create a 3D solid",
     true, "", "A, B", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::extrudeUntil(geom_a->as<SFCGAL::Polygon>(),
                                              *geom_b);
     }},

    {"tesselate", "Construction", "Tesselate a geometry into triangular faces",
     false,
     "No parameters required.\n\nExample:\n  sfcgalop -a \"POLYGON((0 0,3 0,3 "
     "3,0 3,0 0))\" tesselate",
     "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return SFCGAL::algorithm::tesselate(*geom_a);
     }},

    {"triangulate", "Construction",
     "Triangulate a geometry (alias for tesselate)", false,
     "No parameters required.\n\nExample:\n  sfcgalop -a \"POLYGON((0 0,3 0,3 "
     "3,0 3,0 0))\" triangulate",
     "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       return SFCGAL::algorithm::tesselate(*geom_a);
     }},

    {"offset", "Construction", "Create an offset polygon at specified distance",
     false,
     "Parameters:\n  distance=VALUE: Offset distance (default: 1.0)\n  "
     "Positive values create outward offset\n  Negative values create inward "
     "offset\n\nExample:\n  sfcgalop -a \"POLYGON((0 0,3 0,3 3,0 3,0 0))\" "
     "offset \"0.5\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params   = parse_params(args);
       double distance = params.count("distance") ? params["distance"] : 1.0;
       return SFCGAL::algorithm::offset(*geom_a, distance);
     }},

    {"buffer3d", "Construction", "Create a 3D buffer around points and lines",
     false,
     "Parameters:\n  radius=VALUE: Buffer radius (default: 1.0)\n\nNote: "
     "Segments are fixed at 16\nOnly works with Point and LineString "
     "geometries\n\nExample:\n  sfcgalop -a \"POINT(0 0 0)\" buffer3d "
     "\"radius=2.5\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params   = parse_params(args);
       double radius   = params.count("radius") ? params["radius"] : 1.0;
       int    segments = 16; // default segments

       try {
         SFCGAL::algorithm::Buffer3D buffer(*geom_a, radius, segments);
         return buffer.compute(SFCGAL::algorithm::Buffer3D::ROUND);
       } catch (const std::exception &) {
         // Buffer3D only works for Point and LineString
         return std::nullopt;
       }
     }},

    {"minkowskisum", "Construction", "Compute Minkowski sum of two geometries",
     true, "", "A, B", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       if (geom_b->geometryTypeId() == SFCGAL::TYPE_POLYGON) {
         return SFCGAL::algorithm::minkowskiSum(*geom_a,
                                                geom_b->as<SFCGAL::Polygon>());
       }
       return std::nullopt;
     }},

    {"minkowskisum3d", "Construction",
     "Compute 3D Minkowski sum of two geometries", true, "", "A, B", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       return SFCGAL::algorithm::minkowskiSum3D(*geom_a, *geom_b);
     }},

#ifndef _MSC_VER
    {"alphashapes", "Construction", "Compute alpha shapes from point cloud",
     false,
     "Parameters:\n  alpha=VALUE: Alpha parameter controlling shape detail "
     "(default: 1.0)\n  Smaller values create more detailed "
     "shapes\n\nExample:\n  sfcgalop -a \"MULTIPOINT((0 0),(1 0),(0.5 1),(2 "
     "0.5))\" alphashapes \"alpha=0.5\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params     = parse_params(args);
       double alpha      = params.count("alpha") ? params["alpha"] : 1.0;
       bool   allowHoles = true;
       return SFCGAL::algorithm::alphaShapes(*geom_a, alpha, allowHoles);
     }},
#endif

    {"alphawrapping3d", "Construction",
     "Create 3D alpha wrapping surface from points", false,
     "Parameters:\n  alpha=VALUE: Alpha parameter (default: 1.0)\n  "
     "offset=VALUE: Offset parameter (default: 0)\n\nExample:\n  sfcgalop -a "
     "\"MULTIPOINT Z((0 0 0),(1 0 0),(0 1 0),(0 0 1))\" alphawrapping3d "
     "\"alpha=1.5,offset=0\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double alpha  = params.count("alpha") ? params["alpha"] : 1.0;
       size_t offset =
           static_cast<size_t>(params.count("offset") ? params["offset"] : 0.0);
       auto alphaInt = static_cast<size_t>(
           alpha * 100); // Convert to integer representation
       return SFCGAL::algorithm::alphaWrapping3D(*geom_a, alphaInt, offset);
     }},

    {"linesubstring", "Construction",
     "Extract substring from linestring by fraction", false,
     "Parameters:\n  start=VALUE: Start fraction (0.0 to 1.0)\n  end=VALUE: "
     "End fraction (default: 1.0)\n\nExample:\n  sfcgalop -a \"LINESTRING(0 "
     "0,10 0,10 10)\" linesubstring \"start=0.25,end=0.75\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double start  = params["start"];
       double end    = params.count("end") ? params["end"] : 1.0;
       if (geom_a->geometryTypeId() == SFCGAL::TYPE_LINESTRING) {
         const auto &lineString = geom_a->as<SFCGAL::LineString>();
         return SFCGAL::algorithm::lineSubstring(lineString, start, end);
       }
       return std::nullopt;
     }},

    // Transformations
    {"translate", "Transformations", "Translate geometry by specified offset",
     false,
     "Parameters:\n  dx=X: X-axis translation\n  dy=Y: Y-axis translation\n  "
     "dz=Z: Z-axis translation\n\nExample:\n  sfcgalop -a \"POINT(0 0)\" "
     "translate \"dx=5,dy=3,dz=1\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double dx     = params["dx"];
       double dy     = params["dy"];
       double dz     = params["dz"];
       auto   result = geom_a->clone();
       SFCGAL::algorithm::translate(*result, dx, dy, dz);
       return result;
     }},

    {"rotate", "Transformations", "Rotate geometry around specified axis",
     false,
     "Parameters:\n  angle=DEGREES: Rotation angle in degrees\n  axis=x|y|z: "
     "Rotation axis (default: z)\n\nExamples:\n  sfcgalop -a \"POINT(1 0)\" "
     "rotate \"angle=90\"           # Rotate 90° around Z-axis\n  sfcgalop -a "
     "\"POINT(1 0 0)\" rotate \"angle=90,axis=y\"  # Rotate 90° around Y-axis",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto        params    = parse_params(args);
       double      angle_deg = params.count("angle") ? params["angle"] : 0.0;
       std::string axis      = "z"; // default to Z-axis

       // Check if axis parameter is provided (axis won't be in parse_params
       // because it's not a double) Parse manually for non-numeric parameters
       if (args.find("axis=") != std::string::npos) {
         size_t axis_pos = args.find("axis=");
         size_t start    = axis_pos + 5; // length of "axis="
         size_t end      = args.find(',', start);
         if (end == std::string::npos) {
           end = args.length();
         }
         axis = args.substr(start, end - start);
       }

       // Convert degrees to radians
       double angle_rad = (angle_deg * M_PI) / 180.0;

       auto result = geom_a->clone();

       try {
         if (axis == "x" || axis == "X") {
           SFCGAL::algorithm::rotateX(*result, angle_rad);
         } else if (axis == "y" || axis == "Y") {
           SFCGAL::algorithm::rotateY(*result, angle_rad);
         } else {
           // Default to Z-axis rotation
           SFCGAL::algorithm::rotateZ(*result, angle_rad);
         }
         return result;
       } catch (const std::exception &e) {
         std::cerr << "Rotation failed: " << e.what() << "\n";
         return std::nullopt;
       }
     }},

    {"scale", "Transformations", "Scale geometry by specified factors", false,
     "Parameters:\n  s=VALUE: Uniform scaling factor\n  OR\n  sx=X: X-axis "
     "scaling factor\n  sy=Y: Y-axis scaling factor\n  sz=Z: Z-axis scaling "
     "factor\n\nExamples:\n  sfcgalop -a \"POINT(1 1)\" scale \"s=2\"\n  "
     "sfcgalop -a \"POINT(1 1)\" scale \"sx=2,sy=1,sz=1\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double scale  = params.count("s") ? params["s"] : 1.0;
       double sx     = params.count("sx") ? params["sx"] : scale;
       double sy     = params.count("sy") ? params["sy"] : scale;
       double sz     = params.count("sz") ? params["sz"] : scale;
       auto   result = geom_a->clone();
       SFCGAL::algorithm::scale(*result, sx, sy, sz);
       return result;
     }},

    {"force2d", "Transformations", "Remove Z coordinates to create 2D geometry",
     false,
     "No parameters required.\n\nExample:\n  sfcgalop -a \"POINT Z(1 2 3)\" "
     "force2d",
     "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto result = geom_a->clone();
       SFCGAL::algorithm::force2D(*result);
       return result;
     }},

    {"force3d", "Transformations", "Add Z coordinates to create 3D geometry",
     false,
     "Parameters:\n  z=VALUE: Z coordinate value to assign (default: "
     "0.0)\n\nExample:\n  sfcgalop -a \"POINT(1 2)\" force3d \"z=5.0\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double z      = params.count("z") ? params["z"] : 1.0;
       auto   result = geom_a->clone();
       SFCGAL::algorithm::force3D(*result, z);
       return result;
     }},

    {"forcemeasured", "Transformations", "Add measure coordinates to geometry",
     false,
     "Parameters:\n  m=VALUE: Measure coordinate value (default: "
     "0.0)\n\nExample:\n  sfcgalop -a \"POINT(1 2)\" forcemeasured \"m=10.0\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double m      = params.count("m") ? params["m"] : 1.0;
       auto   result = geom_a->clone();
       SFCGAL::algorithm::forceMeasured(*result, m);
       return result;
     }},

    {"simplify", "Transformations",
     "Simplify geometry by removing vertices within tolerance", false,
     "Parameters:\n  tolerance=VALUE: Distance tolerance for vertex removal "
     "(default: 0.01)\n\nExample:\n  sfcgalop -a \"LINESTRING(0 0,0.01 0.01,1 "
     "1,2 2)\" simplify \"tolerance=0.1\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params    = parse_params(args);
       double tolerance = params.count("tolerance") ? params["tolerance"] : 1.0;
       bool   preserveTopology = true;
       return SFCGAL::algorithm::simplify(*geom_a, tolerance, preserveTopology);
     }},

#if SFCGAL_CGAL_VERSION_MAJOR >= 6
    {"polygonrepair", "Transformations", "Repair invalid polygons with rules",
     false,
     "Parameters:\n  method=0|1|2|3 (default: 0)\n\nMethods:\n  0 = "
     "Even-odd\n  1 = Non-zero winding\n  2 = Union of all polygons\n  3 = "
     "Intersection of all polygons\n\n"
     "Example:\n  sfcgalop -a \"POLYGON((0 0, 2 2, 2 0, 0 2, 0 0)) "
     "polygonrepair \"method=1\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto params = parse_params(args);
       int  method =
           static_cast<int>(params.count("method") ? params["method"] : 0);

       SFCGAL::algorithm::PolygonRepairRule rule;
       switch (method) {
       case 1:
         rule = SFCGAL::algorithm::PolygonRepairRule::NON_ZERO_RULE;
         break;
       case 2:
         rule = SFCGAL::algorithm::PolygonRepairRule::UNION_RULE;
         break;
       case 3:
         rule = SFCGAL::algorithm::PolygonRepairRule::INTERSECTION_RULE;
         break;
       default:
         rule = SFCGAL::algorithm::PolygonRepairRule::EVEN_ODD_RULE;
       }

       return SFCGAL::algorithm::polygonRepair(*geom_a, rule);
     }},
#endif // SFCGAL_CGAL_VERSION_MAJOR >= 6

    // Collection operations
    {"collect", "Collections", "Combine two geometries into a collection", true,
     "", "A, B", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *geom_b) -> std::optional<OperationResult> {
       if (!geom_b) {
         return std::nullopt;
       }
       auto collection = std::make_unique<SFCGAL::GeometryCollection>();
       collection->addGeometry(*geom_a);
       collection->addGeometry(*geom_b);
       return collection;
     }},

    {"collection_extract", "Collections",
     "Extract polygons from geometry collection", false, "", "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto copy = geom_a->clone();
       return SFCGAL::algorithm::collectionExtractPolygons(std::move(copy));
     }},

    {"collection_homogenize", "Collections",
     "Convert collection to appropriate multi-type", false, "", "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto copy = geom_a->clone();
       return SFCGAL::algorithm::collectionHomogenize(std::move(copy));
     }},

    {"collection_to_multi", "Collections",
     "Convert collection to multi-geometry type", false, "", "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto copy = geom_a->clone();
       return SFCGAL::algorithm::collectionToMulti(std::move(copy));
     }},

    // Analysis
    {"orientation", "Analysis",
     "Determine polygon ring orientation (clockwise/counter-clockwise)", false,
     "", "A", "D",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       if (geom_a->geometryTypeId() == SFCGAL::TYPE_POLYGON) {
         const auto &polygon = geom_a->as<SFCGAL::Polygon>();
         // Check if polygon exterior ring is counter-clockwise
         bool ccw = SFCGAL::algorithm::isCounterClockWiseOriented(
             polygon.exteriorRing());
         return static_cast<double>(ccw ? 1 : -1);
       }
       return std::nullopt;
     }},

    {"visibility", "Analysis",
     "Compute visibility polygon from a point in polygon", false,
     "Parameters:\n  x=X_COORD: X coordinate of viewpoint\n  y=Y_COORD: Y "
     "coordinate of viewpoint\n\nExample:\n  sfcgalop -a \"POLYGON((0 0,10 "
     "0,10 10,0 10,0 0))\" visibility \"x=5,y=5\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       if (geom_a->geometryTypeId() == SFCGAL::TYPE_POLYGON) {
         const auto   &polygon = geom_a->as<SFCGAL::Polygon>();
         auto          params  = parse_params(args);
         SFCGAL::Point point(params["x"], params["y"]);
         return SFCGAL::algorithm::visibility(polygon, point);
       }
       return std::nullopt;
     }},

    {"partition", "Analysis", "Partition polygon into simpler pieces", false,
     "Parameters:\n  method=0|1|2|3 (default: 0)\n\nMethods:\n  0 = "
     "y_monotone: Creates y-monotone polygons\n  1 = approx_convex: "
     "Approximate convex partition\n  2 = greene_approx_convex: Greene's "
     "approximation algorithm\n  3 = optimal_convex: Optimal convex "
     "partition\n\nExample:\n  sfcgalop -a \"POLYGON((0 0,5 5,10 0,5 -2,0 "
     "0))\" partition \"method=1\"",
     "A, params", "G",
     [](const std::string &args, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto params = parse_params(args);
       int  method =
           static_cast<int>(params.count("method") ? params["method"] : 0);

       SFCGAL::algorithm::PartitionAlgorithm alg;
       switch (method) {
       case 1:
         alg = SFCGAL::algorithm::approx_convex;
         break;
       case 2:
         alg = SFCGAL::algorithm::greene_approx_convex;
         break;
       case 3:
         alg = SFCGAL::algorithm::optimal_convex;
         break;
       default:
         alg = SFCGAL::algorithm::y_monotone;
       }

       return SFCGAL::algorithm::partition_2(*geom_a, alg);
     }},

    {"normal", "Analysis", "Compute surface normal vector for polygon/triangle",
     false, "", "A", "G",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       if (geom_a->geometryTypeId() == SFCGAL::TYPE_POLYGON ||
           geom_a->geometryTypeId() == SFCGAL::TYPE_TRIANGLE) {
         // Compute normal using first 3 points of polygon
         const auto &polygon =
             (geom_a->geometryTypeId() == SFCGAL::TYPE_POLYGON)
                 ? geom_a->as<SFCGAL::Polygon>()
                 : SFCGAL::Polygon(geom_a->as<SFCGAL::Triangle>());

         if (polygon.exteriorRing().numPoints() >= 3) {
           const auto &point0 = polygon.exteriorRing().pointN(0);
           const auto &point1 = polygon.exteriorRing().pointN(1);
           const auto &point2 = polygon.exteriorRing().pointN(2);

           // Calculate normal as cross product
           double v1x = CGAL::to_double(point1.x() - point0.x());
           double v1y = CGAL::to_double(point1.y() - point0.y());
           double v1z =
               point1.is3D() ? CGAL::to_double(point1.z() - point0.z()) : 0;

           double v2x = CGAL::to_double(point2.x() - point0.x());
           double v2y = CGAL::to_double(point2.y() - point0.y());
           double v2z =
               point2.is3D() ? CGAL::to_double(point2.z() - point0.z()) : 0;

           double nx = (v1y * v2z) - (v1z * v2y);
           double ny = (v1z * v2x) - (v1x * v2z);
           double nz = (v1x * v2y) - (v1y * v2x);

           return std::make_unique<SFCGAL::Point>(nx, ny, nz);
         }
       }
       return std::nullopt;
     }},

    // Constructors
    {"make_sphere", "Constructors", "Create a 3D sphere primitive", false,
     "Parameters:\n  x=X_COORD: X coordinate of center (default: 0.0)\n  "
     "y=Y_COORD: Y coordinate of center (default: 0.0)\n  z=Z_COORD: Z "
     "coordinate of center (default: 0.0)\n  radius=VALUE: Sphere radius "
     "(default: 1.0)\n  num_vertical=N: Number of vertical divisions (default: "
     "16)\n  num_horizontal=N: Number of horizontal divisions (default: "
     "32)\n\nExample:\n  sfcgalop make_sphere "
     "\"x=0,y=0,z=0,radius=2.5,num_vertical=20,num_horizontal=40\"",
     "params", "G",
     [](const std::string &args, const SFCGAL::Geometry *,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params       = parse_params(args);
       double x            = params.count("x") ? params["x"] : 0.0;
       double y            = params.count("y") ? params["y"] : 0.0;
       double z            = params.count("z") ? params["z"] : 0.0;
       double radius       = params.count("radius") ? params["radius"] : 1.0;
       auto   num_vertical = static_cast<unsigned int>(
           params.count("num_vertical") ? params["num_vertical"] : 16);
       auto num_horizontal = static_cast<unsigned int>(
           params.count("num_horizontal") ? params["num_horizontal"] : 32);
       return Constructors::make_sphere(x, y, z, radius, num_vertical,
                                        num_horizontal);
     }},

    {"make_cube", "Constructors", "Create a 3D cube primitive", false,
     "Parameters:\n  size=VALUE: Edge length of the cube (default: "
     "1.0)\n\nExample:\n  sfcgalop make_cube \"size=2.0\"",
     "params", "G",
     [](const std::string &args, const SFCGAL::Geometry *,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double size   = params.count("size") ? params["size"] : 1.0;
       return Constructors::make_cube(size);
     }},

    {"make_box", "Constructors", "Create a 3D box primitive", false,
     "Parameters:\n  x_extent=VALUE: Length in X direction (default: 1.0)\n  "
     "y_extent=VALUE: Length in Y direction (default: 1.0)\n  z_extent=VALUE: "
     "Length in Z direction (default: 1.0)\n\nExample:\n  sfcgalop make_box "
     "\"x_extent=2,y_extent=3,z_extent=1\"",
     "params", "G",
     [](const std::string &args, const SFCGAL::Geometry *,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params   = parse_params(args);
       double x_extent = params.count("x_extent") ? params["x_extent"] : 1.0;
       double y_extent = params.count("y_extent") ? params["y_extent"] : 1.0;
       double z_extent = params.count("z_extent") ? params["z_extent"] : 1.0;
       return Constructors::make_box(x_extent, y_extent, z_extent);
     }},

    {"make_cylinder", "Constructors", "Create a 3D cylinder primitive", false,
     "Parameters:\n  base_x=VALUE: X coordinate of base center (default: "
     "0.0)\n  base_y=VALUE: Y coordinate of base center (default: 0.0)\n  "
     "base_z=VALUE: Z coordinate of base center (default: 0.0)\n  "
     "axis_x=VALUE: X component of cylinder axis (default: 0.0)\n  "
     "axis_y=VALUE: Y component of cylinder axis (default: 0.0)\n  "
     "axis_z=VALUE: Z component of cylinder axis (default: 1.0)\n  "
     "radius=VALUE: Cylinder radius (default: 1.0)\n  height=VALUE: Cylinder "
     "height (default: 1.0)\n  num_radial=N: Number of radial divisions "
     "(default: 32)\n\nExample:\n  sfcgalop make_cylinder "
     "\"radius=1.5,height=3,num_radial=16\"",
     "params", "G",
     [](const std::string &args, const SFCGAL::Geometry *,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params     = parse_params(args);
       double base_x     = params.count("base_x") ? params["base_x"] : 0.0;
       double base_y     = params.count("base_y") ? params["base_y"] : 0.0;
       double base_z     = params.count("base_z") ? params["base_z"] : 0.0;
       double axis_x     = params.count("axis_x") ? params["axis_x"] : 0.0;
       double axis_y     = params.count("axis_y") ? params["axis_y"] : 0.0;
       double axis_z     = params.count("axis_z") ? params["axis_z"] : 1.0;
       double radius     = params.count("radius") ? params["radius"] : 1.0;
       double height     = params.count("height") ? params["height"] : 1.0;
       auto   num_radial = static_cast<unsigned int>(
           params.count("num_radial") ? params["num_radial"] : 32);
       return Constructors::make_cylinder(base_x, base_y, base_z, axis_x,
                                          axis_y, axis_z, radius, height,
                                          num_radial);
     }},

    {"make_cone", "Constructors",
     "Create a 3D cone primitive (supports truncated cones)", false,
     "Parameters:\n  base_x=VALUE: X coordinate of base center (default: "
     "0.0)\n  base_y=VALUE: Y coordinate of base center (default: 0.0)\n  "
     "base_z=VALUE: Z coordinate of base center (default: 0.0)\n  "
     "axis_x=VALUE: X component of cone axis (default: 0.0)\n  axis_y=VALUE: Y "
     "component of cone axis (default: 0.0)\n  axis_z=VALUE: Z component of "
     "cone axis (default: 1.0)\n  bottom_radius=VALUE: Cone bottom radius "
     "(default: "
     "1.0)\n  top_radius=VALUE: Cone top radius - 0.0 for regular cone "
     "(default: 0.0)\n  "
     "height=VALUE: Cone height (default: 1.0)\n  num_radial=N: Number "
     "of radial divisions (default: 32)\n\nExamples:\n  sfcgalop make_cone "
     "\"bottom_radius=2,height=4,num_radial=24\"\n  sfcgalop make_cone "
     "\"bottom_radius=3,top_radius=1,height=5\" # Truncated cone",
     "params", "G",
     [](const std::string &args, const SFCGAL::Geometry *,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params = parse_params(args);
       double base_x = params.count("base_x") ? params["base_x"] : 0.0;
       double base_y = params.count("base_y") ? params["base_y"] : 0.0;
       double base_z = params.count("base_z") ? params["base_z"] : 0.0;
       double axis_x = params.count("axis_x") ? params["axis_x"] : 0.0;
       double axis_y = params.count("axis_y") ? params["axis_y"] : 0.0;
       double axis_z = params.count("axis_z") ? params["axis_z"] : 1.0;
       double bottom_radius =
           params.count("bottom_radius") ? params["bottom_radius"] : 1.0;
       double top_radius =
           params.count("top_radius") ? params["top_radius"] : 0.0;
       double height     = params.count("height") ? params["height"] : 1.0;
       auto   num_radial = static_cast<unsigned int>(
           params.count("num_radial") ? params["num_radial"] : 32);
       return Constructors::make_cone(base_x, base_y, base_z, axis_x, axis_y,
                                      axis_z, bottom_radius, top_radius, height,
                                      num_radial);
     }},

    {"make_torus", "Constructors", "Create a 3D torus primitive", false,
     "Parameters:\n  center_x=VALUE: X coordinate of center (default: 0.0)\n  "
     "center_y=VALUE: Y coordinate of center (default: 0.0)\n  center_z=VALUE: "
     "Z coordinate of center (default: 0.0)\n  axis_x=VALUE: X component of "
     "torus axis (default: 0.0)\n  axis_y=VALUE: Y component of torus axis "
     "(default: 0.0)\n  axis_z=VALUE: Z component of torus axis (default: "
     "1.0)\n  major_radius=VALUE: Major radius of torus (default: 2.0)\n  "
     "minor_radius=VALUE: Minor radius of torus (default: 0.5)\n  num_major=N: "
     "Number of major divisions (default: 32)\n  num_minor=N: Number of minor "
     "divisions (default: 16)\n\nExample:\n  sfcgalop make_torus "
     "\"major_radius=3,minor_radius=1,num_major=24,num_minor=12\"",
     "params", "G",
     [](const std::string &args, const SFCGAL::Geometry *,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       auto   params   = parse_params(args);
       double center_x = params.count("center_x") ? params["center_x"] : 0.0;
       double center_y = params.count("center_y") ? params["center_y"] : 0.0;
       double center_z = params.count("center_z") ? params["center_z"] : 0.0;
       double axis_x   = params.count("axis_x") ? params["axis_x"] : 0.0;
       double axis_y   = params.count("axis_y") ? params["axis_y"] : 0.0;
       double axis_z   = params.count("axis_z") ? params["axis_z"] : 1.0;
       double major_radius =
           params.count("major_radius") ? params["major_radius"] : 2.0;
       double minor_radius =
           params.count("minor_radius") ? params["minor_radius"] : 0.5;
       auto num_major = static_cast<unsigned int>(
           params.count("num_major") ? params["num_major"] : 32);
       auto num_minor = static_cast<unsigned int>(
           params.count("num_minor") ? params["num_minor"] : 16);
       return Constructors::make_torus(center_x, center_y, center_z, axis_x,
                                       axis_y, axis_z, major_radius,
                                       minor_radius, num_major, num_minor);
     }},

    {"to_solid", "Conversions", "Convert a PolyhedralSurface to a Solid", false,
     "", "A", "",
     [](const std::string &, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       // Clone the input geometry to pass ownership to make_solid
       auto geom_copy = geom_a->clone();
       return Constructors::make_solid(std::move(geom_copy));
     }},

    // Roof Generation
    {"generate_gable_roof", "Roof Generation",
     "Generate a gable roof using automatic medial axis ridge", false,
     "Parameters:\n"
     "  slope_angle=ANGLE         - Roof slope angle in degrees (required)\n"
     "  building_height=HEIGHT    - Building height in units (optional, "
     "default=0.0)\n"
     "  add_vertical_faces=BOOL   - Add vertical triangular faces at ridge "
     "endpoints (optional, default=false)\n"
     "                              Accepts: true/false, t/f, 1/0, "
     "TRUE/FALSE "
     "(case-insensitive)\n\n"
     "Examples:\n"
     "  sfcgalop -a \"POLYGON((0 0,10 0,10 6,0 6,0 0))\" "
     "generate_gable_roof "
     "\"slope_angle=30.0\"\n"
     "  sfcgalop -a \"POLYGON((0 0,10 0,10 6,0 6,0 0))\" "
     "generate_gable_roof "
     "\"slope_angle=30.0,building_height=5.0,add_vertical_faces=true\"",
     "A, params", "G",
     [](const std::string &params_str, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       const auto *polygon = dynamic_cast<const SFCGAL::Polygon *>(geom_a);
       if (!polygon) {
         std::cerr << "Error: generate_gable_roof requires a POLYGON input\n";
         return std::nullopt;
       }

       auto params = parse_params(params_str);

       auto slope_angle_it = params.find("slope_angle");
       if (slope_angle_it == params.end()) {
         std::cerr << "Error: slope_angle parameter is required\n";
         return std::nullopt;
       }
       double slope_angle = slope_angle_it->second;

       double building_height    = 0.0;
       auto   building_height_it = params.find("building_height");
       if (building_height_it != params.end()) {
         building_height = building_height_it->second;
       }

       bool add_vertical_faces =
           parse_boolean_param(params, "add_vertical_faces", params_str, false);

       SFCGAL::algorithm::RoofParameters roof_params;
       roof_params.type             = SFCGAL::algorithm::RoofType::GABLE;
       roof_params.slopeAngle       = slope_angle;
       roof_params.addVerticalFaces = add_vertical_faces;
       roof_params.buildingHeight   = building_height;
       return SFCGAL::algorithm::generateRoof(*polygon, SFCGAL::LineString(),
                                              roof_params);
     }},

    {"generate_flat_roof", "Roof Generation",
     "Generate a flat roof by extruding a polygon", false,
     "Parameters:\n"
     "  height=VALUE              - Roof height in units (required)\n"
     "  building_height=HEIGHT    - Building height in units (optional, "
     "default=0.0)\n\n"
     "Examples:\n"
     "  sfcgalop -a \"POLYGON((0 0,10 0,10 6,0 6,0 0))\" "
     "generate_flat_roof "
     "\"height=2.0\"\n"
     "  sfcgalop -a \"POLYGON((0 0,10 0,10 6,0 6,0 0))\" "
     "generate_flat_roof "
     "\"height=2.0,building_height=5.0\"",
     "A, params", "G",
     [](const std::string &params_str, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       const auto *polygon = dynamic_cast<const SFCGAL::Polygon *>(geom_a);
       if (!polygon) {
         std::cerr << "Error: generate_flat_roof requires a POLYGON input\n";
         return std::nullopt;
       }

       auto params = parse_params(params_str);

       auto height_it = params.find("height");
       if (height_it == params.end()) {
         std::cerr << "Error: height parameter is required\n";
         return std::nullopt;
       }
       double roof_height = height_it->second;

       double building_height    = 0.0;
       auto   building_height_it = params.find("building_height");
       if (building_height_it != params.end()) {
         building_height = building_height_it->second;
       }

       SFCGAL::algorithm::RoofParameters roof_params;
       roof_params.type           = SFCGAL::algorithm::RoofType::FLAT;
       roof_params.roofHeight     = roof_height;
       roof_params.buildingHeight = building_height;
       return SFCGAL::algorithm::generateRoof(*polygon, SFCGAL::LineString(),
                                              roof_params);
     }},

    {"generate_hipped_roof", "Roof Generation",
     "Generate a hipped roof using straight skeleton extrusion", false,
     "Parameters:\n"
     "  height=VALUE              - Roof height in units (required)\n"
     "  building_height=HEIGHT    - Building height in units (optional, "
     "default=0.0)\n"
     "  angle1=VALUE    - Roof angle for the first edge\n"
     "  angleN=VALUE    - Roof angle for the Nth edge\n\n"
     "Examples:\n"
     "  sfcgalop -a \"POLYGON((0 0,10 0,10 6,0 6,0 0))\" "
     "generate_hipped_roof "
     "\"height=3.0\"\n"
     "  sfcgalop -a \"POLYGON((0 0,10 0,10 6,0 6,0 0))\" "
     "generate_hipped_roof "
     "\"height=3.0,building_height=5.0,angle1=30.0,angle2=45.0,angle3=30.0,"
     "angle4=45.0\"",
     "A, params", "G",
     [](const std::string &params_str, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       const auto *polygon = dynamic_cast<const SFCGAL::Polygon *>(geom_a);
       if (!polygon) {
         std::cerr << "Error: generate_hipped_roof requires a POLYGON input\n";
         return std::nullopt;
       }

       auto params = parse_params(params_str);

       auto height_it = params.find("height");
       if (height_it == params.end()) {
         std::cerr << "Error: height parameter is required\n";
         return std::nullopt;
       }
       double roof_height = height_it->second;

       double building_height    = 0.0;
       auto   building_height_it = params.find("building_height");
       if (building_height_it != params.end()) {
         building_height = building_height_it->second;
       }

       std::vector<SFCGAL::Kernel::FT> angles;
       int                             angleIdx = 0;
       while (true) {
         angleIdx++;
         auto angle_it = params.find("angle" + std::to_string(angleIdx));
         if (angle_it == params.end()) {
           break;
         }
         angles.emplace_back(angle_it->second);
       };

       SFCGAL::algorithm::RoofParameters roof_params;
       roof_params.type           = SFCGAL::algorithm::RoofType::HIPPED;
       roof_params.roofHeight     = roof_height;
       roof_params.buildingHeight = building_height;
       roof_params.angles.emplace_back(angles);
       return SFCGAL::algorithm::generateRoof(*polygon, SFCGAL::LineString(),
                                              roof_params);
     }

    },

    {"generate_skillion_roof", "Roof Generation",
     "Generate a skillion (shed) roof with single slope from ridge "
     "line",
     false,
     "Parameters:\n"
     "  slope_angle=ANGLE         - Roof slope angle in degrees "
     "(required)\n"
     "  add_vertical_faces=BOOL   - Add vertical triangular faces at "
     "roof ends "
     "(optional, default=false)\\n"
     "                              Accepts: true/false, t/f, 1/0, "
     "TRUE/FALSE "
     "(case-insensitive)\n"
     "  building_height=HEIGHT    - Building height in units "
     "(optional, "
     "default=0.0)\n"
     "  ridge_edge=INDEX          - Edge index to use as ridge line "
     "(optional, "
     "default=0)\n\n"
     "Examples:\n"
     "  sfcgalop -a \"POLYGON((0 0,10 0,10 6,0 6,0 0))\" "
     "generate_skillion_roof \"slope_angle=15.0\"\n"
     "  sfcgalop -a \"POLYGON((0 0,10 0,10 6,0 6,0 0))\" "
     "generate_skillion_roof "
     "\"slope_angle=15.0,add_vertical_faces=true\"\n"
     "  sfcgalop -a \"POLYGON((0 0,10 0,10 6,0 6,0 0))\" "
     "generate_skillion_roof "
     "\"slope_angle=15.0,building_height=3.0,ridge_edge=1\"",
     "A, params", "G",
     [](const std::string &params_str, const SFCGAL::Geometry *geom_a,
        const SFCGAL::Geometry *) -> std::optional<OperationResult> {
       const auto *polygon = dynamic_cast<const SFCGAL::Polygon *>(geom_a);
       if (!polygon) {
         std::cerr << "Error: generate_skillion_roof requires a "
                      "POLYGON input\n";
         return std::nullopt;
       }

       auto params = parse_params(params_str);

       auto slope_angle_it = params.find("slope_angle");
       if (slope_angle_it == params.end()) {
         std::cerr << "Error: slope_angle parameter is required for skillion "
                      "roof generation\n";
         std::cerr << "Usage: generate_skillion_roof "
                      "\"slope_angle=30.0,add_vertical_faces=1\"\n";
         std::cerr << "Valid range: 0 < slope_angle < 90 degrees\n";
         return std::nullopt;
       }

       double slope_angle = slope_angle_it->second;
       if (slope_angle <= 0.0 || slope_angle >= 90.0) {
         std::cerr << "Error: slope_angle must be between 0 and 90 "
                      "degrees, got: "
                   << slope_angle << "\n";
         std::cerr << "Common values: 15° (low slope), 30° (moderate), 45° "
                      "(steep)\n";
         return std::nullopt;
       }

       bool add_vertical_faces =
           parse_boolean_param(params, "add_vertical_faces", params_str, false);

       double building_height    = 0.0;
       auto   building_height_it = params.find("building_height");
       if (building_height_it != params.end()) {
         building_height = building_height_it->second;
         if (building_height < 0.0) {
           std::cerr << "Error: building_height must be non-negative, got: "
                     << building_height << "\n";
           return std::nullopt;
         }
       }

       int  ridge_edge    = 0;
       auto ridge_edge_it = params.find("ridge_edge");
       if (ridge_edge_it != params.end()) {
         ridge_edge = static_cast<int>(ridge_edge_it->second);
       }

       // Create a ridge line along specified edge of the polygon
       const auto &ring = polygon->exteriorRing();
       if (ring.numPoints() < 4) {
         std::cerr << "Error: Polygon must have at least 3 vertices for roof "
                      "generation\n";
         std::cerr << "Provided polygon has " << (ring.numPoints() - 1)
                   << " vertices\n";
         return std::nullopt;
       }

       // Validate ridge_edge index
       int num_edges = static_cast<int>(ring.numPoints()) - 1;
       if (ridge_edge < 0 || ridge_edge >= num_edges) {
         std::cerr << "Error: ridge_edge must be between 0 and "
                   << (num_edges - 1) << " for polygon with " << num_edges
                   << " edges, got: " << ridge_edge << "\n";
         std::cerr << "Edge numbering: 0=first edge, " << (num_edges - 1)
                   << "=last edge\n";
         return std::nullopt;
       }

       // Use the specified edge as the ridge line (high edge)
       auto               point1 = ring.pointN(ridge_edge);
       auto               point2 = ring.pointN(ridge_edge + 1);
       SFCGAL::LineString ridgeLine(point1, point2);

       SFCGAL::algorithm::RoofParameters roof_params;
       roof_params.type             = SFCGAL::algorithm::RoofType::SKILLION;
       roof_params.slopeAngle       = slope_angle;
       roof_params.addVerticalFaces = add_vertical_faces;
       roof_params.buildingHeight   = building_height;
       return SFCGAL::algorithm::generateRoof(*polygon, ridgeLine, roof_params);
     }}}; // end operations array

} // namespace

namespace Operations {

/**
 * @brief Execute a named geometry operation from the registry.
 *
 * Looks up the operation identified by op_name and, if found, invokes its
 * callable with op_arg and the optional geometry pointers geom_a and geom_b.
 * Exceptions thrown by the operation are caught and result in a nullopt return.
 *
 * @param op_name Name of the operation to execute (must match an entry in the
 * operations registry).
 * @param op_arg Optional argument string passed to the operation (parsed by the
 * operation as needed).
 * @param geom_a First geometry input (may be null if the operation does not
 * require it).
 * @param geom_b Second geometry input (may be null; many operations that
 * require a second geometry will return nullopt when this is absent).
 * @return std::optional<OperationResult> The operation's result if the
 * operation exists and completes successfully; std::nullopt if the operation
 * name is not found or if an exception occurs during execution.
 */
auto
execute_operation(const std::string &op_name, const std::string &op_arg,
                  const SFCGAL::Geometry *geom_a,
                  const SFCGAL::Geometry *geom_b)
    -> std::optional<OperationResult>
{

  auto operation_it =
      std::find_if(operations.begin(), operations.end(),
                   [&op_name](const Operation &operation) -> bool {
                     return operation.name == op_name;
                   });

  if (operation_it != operations.end()) {
    // Check for null geometry A before calling operation, but allow constructor
    // operations
    bool is_constructor = operation_it->name.find("make_") == 0;
    if (geom_a == nullptr && !is_constructor) {
      return std::nullopt;
    }

    try {
      return operation_it->func(op_arg, geom_a, geom_b);
    } catch (const std::exception &e) {
      std::cerr << "Operation error: " << e.what() << "\n";
      return std::nullopt;
    }
  }

  return std::nullopt;
}

/**
 * @brief Prints a categorized list of all registered geometric operations.
 *
 * Outputs each operation name, short description, and a notice when an
 * operation requires two geometries. Operations are grouped by their category
 * and printed to stdout. Uses simple terminal color/formatting for readability.
 */
void
print_available_operations()
{
  std::map<std::string, std::vector<const Operation *>> by_category;

  for (const auto &operation : operations) {
    by_category[operation.category].push_back(&operation);
  }

  // Use simple formatting for now - can be enhanced with TextUI::Table later
  for (const auto &[category, operation_list] : by_category) {
    std::cout << "\n\033[1;36m" << category << ":\033[0m\n"; // Cyan bold
    for (const auto *operation : operation_list) {
      std::cout << "  \033[1m" << std::setw(25) << std::left << operation->name
                << "\033[0m"
                << " " << operation->description;
      if (operation->requires_b) {
        std::cout << " \033[33m(requires 2 geometries)\033[0m"; // Yellow
      }
      std::cout << "\n";
    }
  }
}

/**
 * @brief Print detailed help for a named operation.
 *
 * Looks up an operation by its name and writes its metadata (name, category,
 * description) to stdout. If the operation requires a second geometry this is
 * indicated. If the name is unknown an error message is written to stderr.
 *
 * @param name Null-terminated operation name to look up; if null the function
 * returns false.
 * @return true if the operation was found and printed; false if the name was
 * null or no matching operation exists.
 */
auto
print_operation_help(const char *name) -> bool
{
  if (name == nullptr) {
    return false;
  }

  auto operation_it = std::find_if(operations.begin(), operations.end(),
                                   [name](const Operation &operation) -> bool {
                                     return operation.name == name;
                                   });

  if (operation_it != operations.end()) {
    std::cout << "\nOperation: " << operation_it->name << "\n"
              << "Category: " << operation_it->category << "\n"
              << "Description: " << operation_it->description << "\n";
    if (operation_it->requires_b) {
      std::cout << "Requires two geometries\n";
    }
    if (!operation_it->param_help.empty()) {
      std::cout << "\n" << operation_it->param_help << "\n";
    }
    return true;
  }

  std::cerr << "Unknown operation: " << name << "\n";
  return false;
}

/**
 * @brief Returns metadata for all registered operations.
 *
 * Produces a list of tuples describing every operation in the registry. Each
 * tuple contains:
 *  - name (std::string): operation identifier
 *  - category (std::string): grouping/category name
 *  - description (std::string): human-readable description
 *  - input (std::string): input specification (A, A,B, A,params)
 *  - output (std::string): output type (G, D, B, T)
 *
 * @return std::vector<std::tuple<std::string, std::string, std::string,
 * std::string, std::string>> Vector of operation metadata tuples in the order
 * they appear in the registry.
 */
auto
get_all_operations_info() -> std::vector<
    std::tuple<std::string, std::string, std::string, std::string, std::string>>
{
  std::vector<std::tuple<std::string, std::string, std::string, std::string,
                         std::string>>
      result;
  result.reserve(operations.size());

  for (const auto &operation : operations) {
    result.emplace_back(operation.name, operation.category,
                        operation.description, operation.input,
                        operation.output);
  }

  return result;
}

/**
 * @brief Check if an operation requires a second geometry parameter.
 *
 * Looks up the operation by name in the operations registry and returns
 * whether it requires a second geometry parameter based on the requires_b
 * field in the operation definition.
 *
 * @param operation_name Name of the operation to check.
 * @return true if the operation requires a second geometry; false if it
 * doesn't require one or if the operation name is not found.
 */
auto
operation_requires_second_geometry(const std::string &operation_name) -> bool
{
  auto operation_it =
      std::find_if(operations.begin(), operations.end(),
                   [&operation_name](const Operation &operation) -> bool {
                     return operation.name == operation_name;
                   });

  if (operation_it != operations.end()) {
    return operation_it->requires_b;
  }

  return false; // Operation not found, assume it doesn't require second
                // geometry
}

} // namespace Operations
