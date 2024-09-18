#include "lexicographicOrder.h"
#include <SFCGAL/algorithm/isValid.h>
#include <algorithm>
#include <tuple>

namespace SFCGAL {
namespace transform {

// Utility function to compare two points considering X, Y, Z, and M values
bool
comparePoints(const Point &a, const Point &b)
{
  if (a.x() != b.x())
    return a.x() < b.x();
  if (a.y() != b.y())
    return a.y() < b.y();
  if (a.z() != b.z())
    return a.z() < b.z();
  return a.m() < b.m(); // Compare M if all others are equal
}

// Utility function to find the lowest point of a geometry
Point
findLowestPoint(const Geometry &geom)
{
  switch (geom.geometryTypeId()) {
  case TYPE_POINT:
    return static_cast<const Point &>(geom);
  case TYPE_LINESTRING: {
    const LineString &ls = static_cast<const LineString &>(geom);
    return *std::min_element(ls.begin(), ls.end(), comparePoints);
  }
  case TYPE_POLYGON: {
    const Polygon &poly = static_cast<const Polygon &>(geom);
    return findLowestPoint(poly.exteriorRing());
  }
  case TYPE_TRIANGLE: {
    const Triangle &tri = static_cast<const Triangle &>(geom);
    return *std::min_element(&tri.vertex(0), &tri.vertex(3), comparePoints);
  }
  default:
    // For other types, traverse recursively
    Point lowest = findLowestPoint(geom.geometryN(0));
    for (size_t i = 1; i < geom.numGeometries(); ++i) {
      Point current = findLowestPoint(geom.geometryN(i));
      if (comparePoints(current, lowest))
        lowest = current;
    }
    return lowest;
  }
}

// Implementation for Point
void
lexicographicOrderPoint(Point &, bool reorderRing)
{
  // No modification needed
}

// Implementation for MultiPoint
void
lexicographicOrderMultiPoint(MultiPoint &mp, bool reorderRing)
{
  std::vector<Point> points;
  for (size_t i = 0; i < mp.numGeometries(); ++i) {
    points.push_back(mp.pointN(i));
  }
  if (reorderRing)
    std::sort(points.begin(), points.end(), comparePoints);
  mp = MultiPoint();
  for (const auto &point : points) {
    mp.addGeometry(point);
  }
}

// Implementation for LineString
void
lexicographicOrderLineString(LineString &ls, bool reorderRing)
{
  if (reorderRing) {
    // To reorder the LineString, find the minimal point and rotate the
    // LineString accordingly
    auto minIt = std::min_element(ls.begin(), ls.end(), comparePoints);
    if (minIt != ls.begin()) {
      std::reverse(ls.begin(), ls.end());
    }
  }
}

// Implementation for MultiLineString
void
lexicographicOrderMultiLineString(MultiLineString &mls, bool reorderRing)
{
  std::vector<LineString> lines;
  for (size_t i = 0; i < mls.numGeometries(); ++i) {
    LineString ls = mls.lineStringN(i);
    lexicographicOrderLineString(ls, reorderRing);
    lines.push_back(ls);
  }

  // Sort the LineStrings in MultiLineString based on the lowest point of each
  // LineString
  std::sort(lines.begin(), lines.end(),
            [](const LineString &a, const LineString &b) {
              return comparePoints(findLowestPoint(a), findLowestPoint(b));
            });

  mls = MultiLineString();
  for (const auto &line : lines) {
    mls.addGeometry(line);
  }
}

// Utility function to reorder a ring
LineString
reorderRing(const LineString &ring)
{
  if (ring.numPoints() <= 3) {
    return ring;
  }

  auto   minIt = std::min_element(ring.begin(), ring.end() - 1, comparePoints);
  size_t minIndex = std::distance(ring.begin(), minIt);

  std::vector<Point> newRing;
  newRing.reserve(ring.numPoints());

  newRing.insert(newRing.end(), minIt, ring.end() - 1);
  newRing.insert(newRing.end(), ring.begin(), minIt);

  newRing.push_back(newRing.front());

  return LineString(newRing);
}

// Implementation for Polygon
void
lexicographicOrderPolygon(Polygon &poly, bool doReorderRing)
{
  LineString newExterior = poly.exteriorRing();
  if (doReorderRing) {
    newExterior = reorderRing(poly.exteriorRing());
  }

  std::vector<LineString> interiors;
  for (size_t i = 0; i < poly.numInteriorRings(); ++i) {
    if (doReorderRing) {
      LineString newInterior = reorderRing(poly.interiorRingN(i));
      interiors.push_back(newInterior);
    } else {
      interiors.push_back(poly.interiorRingN(i));
    }
  }

  std::sort(interiors.begin(), interiors.end(),
            [](const LineString &a, const LineString &b) {
              return comparePoints(
                  *std::min_element(a.begin(), a.end(), comparePoints),
                  *std::min_element(b.begin(), b.end(), comparePoints));
            });

  Polygon newPoly(newExterior);
  for (const auto &ring : interiors) {
    newPoly.addInteriorRing(ring);
  }

  poly = newPoly;
}

// Implementation for MultiPolygon
void
lexicographicOrderMultiPolygon(MultiPolygon &mp, bool reorderRing)
{
  std::vector<Polygon> polygons;
  for (size_t i = 0; i < mp.numGeometries(); ++i) {
    Polygon poly = mp.polygonN(i);
    lexicographicOrderPolygon(poly, reorderRing);
    polygons.push_back(poly);
  }
  std::sort(polygons.begin(), polygons.end(),
            [](const Polygon &a, const Polygon &b) {
              return comparePoints(findLowestPoint(a), findLowestPoint(b));
            });
  mp = MultiPolygon();
  for (const auto &poly : polygons) {
    mp.addGeometry(poly);
  }
}

// Implementation for PolyhedralSurface
void
lexicographicOrderPolyhedralSurface(PolyhedralSurface &ps, bool reorderRing)
{
  std::vector<Polygon> polygons;
  for (size_t i = 0; i < ps.numPolygons(); ++i) {
    Polygon poly = ps.polygonN(i);
    lexicographicOrderPolygon(poly, reorderRing);
    polygons.push_back(poly);
  }

  std::sort(polygons.begin(), polygons.end(),
            [](const Polygon &a, const Polygon &b) {
              return comparePoints(findLowestPoint(a), findLowestPoint(b));
            });

  ps = PolyhedralSurface(polygons);
}

// Implementation for Solid
void
lexicographicOrderSolid(Solid &solid, bool reorderRing)
{
  lexicographicOrderPolyhedralSurface(solid.exteriorShell(), reorderRing);
  std::vector<PolyhedralSurface> interiorShells;
  for (size_t i = 0; i < solid.numInteriorShells(); ++i) {
    PolyhedralSurface shell = solid.interiorShellN(i);
    lexicographicOrderPolyhedralSurface(shell, reorderRing);
    interiorShells.push_back(shell);
  }
  std::sort(interiorShells.begin(), interiorShells.end(),
            [](const PolyhedralSurface &a, const PolyhedralSurface &b) {
              return comparePoints(findLowestPoint(a), findLowestPoint(b));
            });

  Solid newSolid(solid.exteriorShell());
  for (const auto &shell : interiorShells) {
    newSolid.addInteriorShell(shell);
  }
  solid = newSolid;
}

// Implementation for GeometryCollection
void
lexicographicOrderGeometryCollection(GeometryCollection &gc, bool reorderRing)
{
  // Appliquer lexicographicOrder à chaque géométrie
  for (size_t i = 0; i < gc.numGeometries(); ++i) {
    lexicographicOrder(gc.geometryN(i), reorderRing);
  }

  // Créer un vecteur de pointeurs vers les géométries
  std::vector<Geometry *> geom_ptrs;
  for (size_t i = 0; i < gc.numGeometries(); ++i) {
    geom_ptrs.push_back(&gc.geometryN(i));
  }

  // Trier les pointeurs vers les géométries par type WKB, puis par leur contenu
  std::sort(geom_ptrs.begin(), geom_ptrs.end(),
            [](const Geometry *a, const Geometry *b) {
              if (a->geometryTypeId() != b->geometryTypeId()) {
                return a->geometryTypeId() < b->geometryTypeId();
              }
              // Si les types sont les mêmes, comparer le contenu
              return a->asText() < b->asText();
            });

  // Reconstruire la GeometryCollection avec les géométries triées
  GeometryCollection newGc;
  for (const auto &geom_ptr : geom_ptrs) {
    newGc.addGeometry(geom_ptr->clone());
  }

  gc = newGc;
}

// Main function
void
lexicographicOrder(Geometry &geometry, bool reorderRing)
{
  switch (geometry.geometryTypeId()) {
  case TYPE_POINT:
    lexicographicOrderPoint(static_cast<Point &>(geometry), reorderRing);
    break;
  case TYPE_MULTIPOINT:
    lexicographicOrderMultiPoint(static_cast<MultiPoint &>(geometry),
                                 reorderRing);
    break;
  case TYPE_LINESTRING:
    lexicographicOrderLineString(static_cast<LineString &>(geometry),
                                 reorderRing);
    break;
  case TYPE_MULTILINESTRING:
    lexicographicOrderMultiLineString(static_cast<MultiLineString &>(geometry),
                                      reorderRing);
    break;
  case TYPE_POLYGON:
    lexicographicOrderPolygon(static_cast<Polygon &>(geometry), reorderRing);
    break;
  case TYPE_MULTIPOLYGON:
    lexicographicOrderMultiPolygon(static_cast<MultiPolygon &>(geometry),
                                   reorderRing);
    break;
  case TYPE_POLYHEDRALSURFACE:
    lexicographicOrderPolyhedralSurface(
        static_cast<PolyhedralSurface &>(geometry), reorderRing);
    break;
  case TYPE_SOLID:
    lexicographicOrderSolid(static_cast<Solid &>(geometry), reorderRing);
    break;
  case TYPE_GEOMETRYCOLLECTION:
    lexicographicOrderGeometryCollection(
        static_cast<GeometryCollection &>(geometry), reorderRing);
    break;

  default:
    // Unsupported geometry, do nothing
    break;
  }
}

} // namespace transform
} // namespace SFCGAL
