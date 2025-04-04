#pragma once

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>

#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <vector>

namespace SFCGAL {
namespace detail {

namespace PS = CGAL::Polyline_simplification_2;

using K   = CGAL::Projection_traits_xy_3<SFCGAL::Kernel>;
using Vb  = PS::Vertex_base_2<K>;
using Fb  = CGAL::Constrained_triangulation_face_base_2<K>;
using TDS = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using CDT =
    CGAL::Constrained_Delaunay_triangulation_2<K, TDS,
                                               CGAL::Exact_predicates_tag>;
using CT                  = CGAL::Constrained_triangulation_plus_2<CDT>;
using CGALPoint           = CT::Point;
using Constraint_id       = CT::Constraint_id;
using Constraint_iterator = CT::Constraint_iterator;
using Points_in_constraint_iterator = CT::Points_in_constraint_iterator;
using Stop                          = PS::Stop_above_cost_threshold;
using Cost                          = PS::Squared_distance_cost;

/**
 * @brief Represents a segment with endpoints to interpolate Z/M values
 */
struct Segment {
  Point p1;
  Point p2;

  Segment(const Point &p1_, const Point &p2_) : p1(p1_), p2(p2_) {}

  // Distance from point to segment (in XY plane)
  double
  distanceToPoint(double x, double y) const;

  // Interpolation parameter for a point projected onto the segment
  double
  interpolationParameter(double x, double y) const;
};

/**
 * @brief Collection of segments from a geometry for interpolation
 */
class SegmentStore {
private:
  std::vector<Segment> segments;
  bool                 hasZCoord = false;
  bool                 hasMCoord = false;

public:
  SegmentStore() {}

  void
  addSegment(const Segment &segment);

  // Check if store contains segments with Z coordinates
  bool
  hasZ() const;

  // Check if store contains segments with M values
  bool
  hasM() const;

  // Find the nearest segment to a point
  Segment
  findNearestSegment(double x, double y) const;

  // Interpolate Z value for a point
  double
  interpolateZ(double x, double y) const;

  // Interpolate M value for a point
  double
  interpolateM(double x, double y) const;
};

/**
 * @brief Extract segments from a LineString for interpolation
 */
void
extractSegments(const LineString &lineString, SegmentStore &store);

/**
 * @brief Extract segments from a Polygon for interpolation
 */
void
extractSegments(const Polygon &polygon, SegmentStore &store);

/**
 * @brief Extract segments from all geometry types
 */
void
extractSegments(const Geometry &geometry, SegmentStore &store);

/**
 * @brief Create a point with interpolated Z and M values
 */
Point
createPointWithInterpolatedZM(double x, double y, double z,
                              const SegmentStore &store,
                              CoordinateType      dimension);

/**
 * @brief Simplifies a LineString
 */
auto
simplifyLineString(const LineString &lineString, double threshold,
                   bool preserveTopology, const SegmentStore &store)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Simplifies a MultiLineString
 */
auto
simplifyMultiLineString(const MultiLineString &multiLine, double threshold,
                        bool preserveTopology) -> std::unique_ptr<Geometry>;

/**
 * @brief Simplifies a Polygon
 */
auto
simplifyPolygon(const Polygon &polygon, double threshold, bool preserveTopology)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Simplifies a MultiPolygon
 */
auto
simplifyMultiPolygon(const MultiPolygon &multiPolygon, double threshold,
                     bool preserveTopology) -> std::unique_ptr<Geometry>;

/**
 * @brief Simplifies a GeometryCollection with topology preservation
 */
auto
simplifyGeometryCollectionTopology(const GeometryCollection &collection,
                                   double                    threshold)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Simplifies a PolyhedralSurface
 */
auto
simplifyPolyhedralSurface(const PolyhedralSurface &polySurface,
                          double threshold, bool preserveTopology)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Simplifies a GeometryCollection
 */
auto
simplifyGeometryCollection(const GeometryCollection &collection,
                           double threshold, bool preserveTopology)
    -> std::unique_ptr<Geometry>;

} // namespace detail
} // namespace SFCGAL
