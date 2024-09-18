#include "SFCGAL/algorithm/compare.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/detail/GeometrySet.h"
#include "SFCGAL/detail/tools/Registry.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace SFCGAL::algorithm::compare {

namespace {

bool
nearEqual(const Kernel::FT &a, const Kernel::FT &b, double epsilon = 1e-9)
{
  return CGAL::abs(CGAL::to_double(a) - CGAL::to_double(b)) < epsilon;
}

template <typename T, typename... Args>
bool
fuzzyEqual(T epsilon, const Args &...args) noexcept
{
  static_assert((sizeof...(args) % 2 == 0 && sizeof...(args) != 0),
                "The number of arguments must be greater than 0 and even");
  constexpr size_t numArgs  = sizeof...(args);
  bool             result   = true;
  T                values[] = {static_cast<T>(CGAL::to_double(args))...};
  for (size_t i = 0; i < numArgs / 2; ++i) {
    result = result && nearEqual(values[i], values[i + numArgs / 2], epsilon);
  }
  return result;
}

template <typename T, typename... Args>
bool
fuzzyDistanceEqual(T epsilon, const Args &...args) noexcept
{
  static_assert((sizeof...(args) % 2 == 0 && sizeof...(args) >= 4),
                "The number of arguments must be greater than 4 and even");
  constexpr size_t numArgs        = sizeof...(args);
  const T          squaredEpsilon = epsilon * epsilon;
  T                sum            = 0;
  T                values[]       = {static_cast<T>(CGAL::to_double(args))...};
  for (size_t i = 0; i < numArgs / 2; ++i) {
    const T diff = values[i] - values[i + numArgs / 2];
    sum += diff * diff;
  }
  return sum < squaredEpsilon;
}

bool
fuzzyComparePoints(const Point &p1, const Point &p2, double epsilon,
                   FuzzyCompareMethod method)
{
  if (p1.isEmpty() != p2.isEmpty())
    return false;
  if (p1.isEmpty())
    return true;

  if (p1.is3D() != p2.is3D())
    return false;
  if (p1.isMeasured() != p2.isMeasured())
    return false;

  if (method == FuzzyCompareMethod::Equal) {
    if (p1.is3D()) {
      if (p1.isMeasured()) {
        return fuzzyEqual(epsilon, p1.x(), p1.y(), p1.z(), p1.m(), p2.x(),
                          p2.y(), p2.z(), p2.m());
      } else {
        return fuzzyEqual(epsilon, p1.x(), p1.y(), p1.z(), p2.x(), p2.y(),
                          p2.z());
      }
    } else {
      if (p1.isMeasured()) {
        return fuzzyEqual(epsilon, p1.x(), p1.y(), p1.m(), p2.x(), p2.y(),
                          p2.m());
      } else {
        return fuzzyEqual(epsilon, p1.x(), p1.y(), p2.x(), p2.y());
      }
    }
  } else {
    if (p1.is3D()) {
      if (p1.isMeasured()) {
        return fuzzyDistanceEqual(epsilon, p1.x(), p1.y(), p1.z(), p1.m(),
                                  p2.x(), p2.y(), p2.z(), p2.m());
      } else {
        return fuzzyDistanceEqual(epsilon, p1.x(), p1.y(), p1.z(), p2.x(),
                                  p2.y(), p2.z());
      }
    } else {
      if (p1.isMeasured()) {
        return fuzzyDistanceEqual(epsilon, p1.x(), p1.y(), p1.m(), p2.x(),
                                  p2.y(), p2.m());
      } else {
        return fuzzyDistanceEqual(epsilon, p1.x(), p1.y(), p2.x(), p2.y());
      }
    }
  }
}

bool
fuzzyCompareLineStrings(const LineString &ls1, const LineString &ls2,
                        double epsilon, FuzzyCompareMethod method)
{
  if (ls1.numPoints() != ls2.numPoints())
    return false;
  for (size_t i = 0; i < ls1.numPoints(); ++i) {
    if (!fuzzyComparePoints(ls1.pointN(i), ls2.pointN(i), epsilon, method)) {
      return false;
    }
  }
  return true;
}

bool
fuzzyComparePolygons(const Polygon &poly1, const Polygon &poly2, double epsilon,
                     FuzzyCompareMethod method)
{
  if (!fuzzyCompareLineStrings(poly1.exteriorRing(), poly2.exteriorRing(),
                               epsilon, method))
    return false;
  if (poly1.numInteriorRings() != poly2.numInteriorRings())
    return false;
  for (size_t i = 0; i < poly1.numInteriorRings(); ++i) {
    if (!fuzzyCompareLineStrings(poly1.interiorRingN(i), poly2.interiorRingN(i),
                                 epsilon, method)) {
      return false;
    }
  }
  return true;
}

bool
fuzzyCompareTriangles(const Triangle &t1, const Triangle &t2, double epsilon,
                      FuzzyCompareMethod method)
{
  return fuzzyComparePoints(t1.vertex(0), t2.vertex(0), epsilon, method) &&
         fuzzyComparePoints(t1.vertex(1), t2.vertex(1), epsilon, method) &&
         fuzzyComparePoints(t1.vertex(2), t2.vertex(2), epsilon, method);
}

// Helper function to find the lowest point in a Geometry
Point
findLowestPoint(const Geometry *geom)
{
  switch (geom->geometryTypeId()) {
  case TYPE_POINT:
    return *static_cast<const Point *>(geom);
  case TYPE_LINESTRING: {
    const auto &ls = static_cast<const LineString *>(geom);
    return *std::min_element(
        ls->begin(), ls->end(), [](const Point &p1, const Point &p2) {
          return p1.x() < p2.x() ||
                 (p1.x() == p2.x() &&
                  (p1.y() < p2.y() || (p1.y() == p2.y() && p1.z() < p2.z())));
        });
  }
  case TYPE_POLYGON: {
    const auto &poly = static_cast<const Polygon *>(geom);
    return findLowestPoint(&(poly->exteriorRing()));
  }
  case TYPE_TRIANGLE: {
    const auto &tri = static_cast<const Triangle *>(geom);
    return *std::min_element(
        &tri->vertex(0), &tri->vertex(0) + 3,
        [](const Point &p1, const Point &p2) {
          return p1.x() < p2.x() ||
                 (p1.x() == p2.x() &&
                  (p1.y() < p2.y() || (p1.y() == p2.y() && p1.z() < p2.z())));
        });
  }
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
  case TYPE_SOLID: {
    const auto &gc = static_cast<const GeometryCollection *>(geom);
    Point       lowest;
    for (size_t i = 0; i < gc->numGeometries(); ++i) {
      Point current = findLowestPoint(&(gc->geometryN(i)));
      if (i == 0 || current.x() < lowest.x() ||
          (current.x() == lowest.x() &&
           (current.y() < lowest.y() ||
            (current.y() == lowest.y() && current.z() < lowest.z())))) {
        lowest = current;
      }
    }
    return lowest;
  }
  default:
    BOOST_THROW_EXCEPTION(Exception(
        (boost::format("Unsupported geometry type %s") % geom->geometryType())
            .str()));
  }
}

// Helper function to sort Polygons based on their lowest point
void
sortPolygon(Polygon &poly)
{
  std::vector<LineString> rings;
  rings.reserve(poly.numInteriorRings());
  for (size_t i = 0; i < poly.numInteriorRings(); ++i) {
    rings.push_back(poly.interiorRingN(i));
  }

  std::sort(rings.begin(), rings.end(),
            [](const LineString &a, const LineString &b) {
              return findLowestPoint(&a) < findLowestPoint(&b);
            });

  poly.removeInteriorRings();
  for (const auto &ring : rings) {
    poly.addInteriorRing(ring);
  }
}

// Helper function to sort GeometryCollections
template <typename T>
void
sortGeometryCollection(T &gc)
{
  std::vector<std::pair<Point, size_t>> sortedIndices;
  sortedIndices.reserve(gc.numGeometries());
  for (size_t i = 0; i < gc.numGeometries(); ++i) {
    sortedIndices.emplace_back(findLowestPoint(&gc.geometryN(i)), i);
  }

  std::sort(sortedIndices.begin(), sortedIndices.end(),
            [](const auto &a, const auto &b) { return a.first < b.first; });

  T sortedGc;
  for (const auto &[_, index] : sortedIndices) {
    sortedGc.addGeometry(gc.geometryN(index).clone());
  }

  gc = std::move(sortedGc);
}

template <typename CompareFunc>
bool
compareGeometries(const Geometry &g1, const Geometry &g2,
                  CompareFunc compareFunc, bool sorted = false)
{
  if (g1.geometryTypeId() != g2.geometryTypeId())
    return false;
  if (g1.isEmpty() && g2.isEmpty())
    return true;
  if (g1.isEmpty() != g2.isEmpty())
    return false;

  switch (g1.geometryTypeId()) {
  case TYPE_POINT:
    return compareFunc(g1.as<Point>(), g2.as<Point>());
  case TYPE_LINESTRING: {
    const auto &ls1 = g1.as<LineString>();
    const auto &ls2 = g2.as<LineString>();
    if (ls1.numPoints() != ls2.numPoints())
      return false;
    for (size_t i = 0; i < ls1.numPoints(); ++i) {
      if (!compareFunc(ls1.pointN(i), ls2.pointN(i)))
        return false;
    }
    return true;
  }
  case TYPE_POLYGON: {
    const auto &poly1 = g1.as<Polygon>();
    const auto &poly2 = g2.as<Polygon>();
    if (!compareGeometries(poly1.exteriorRing(), poly2.exteriorRing(),
                           compareFunc, sorted))
      return false;
    if (poly1.numInteriorRings() != poly2.numInteriorRings())
      return false;
    for (size_t i = 0; i < poly1.numInteriorRings(); ++i) {
      if (!compareGeometries(poly1.interiorRingN(i), poly2.interiorRingN(i),
                             compareFunc, sorted))
        return false;
    }
    return true;
  }
  case TYPE_TRIANGLE:
    return compareFunc(g1.as<Triangle>(), g2.as<Triangle>());
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_MULTISOLID:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE: {
    const auto &gc1 = g1.as<GeometryCollection>();
    const auto &gc2 = g2.as<GeometryCollection>();
    if (gc1.numGeometries() != gc2.numGeometries())
      return false;

    if (sorted) {
      std::vector<const Geometry *> geoms1(gc1.numGeometries());
      std::vector<const Geometry *> geoms2(gc2.numGeometries());

      for (size_t i = 0; i < gc1.numGeometries(); ++i) {
        geoms1[i] = &gc1.geometryN(i);
        geoms2[i] = &gc2.geometryN(i);
      }

      auto geomCompare = [](const Geometry *a, const Geometry *b) {
        return findLowestPoint(a) < findLowestPoint(b);
      };

      std::sort(geoms1.begin(), geoms1.end(), geomCompare);
      std::sort(geoms2.begin(), geoms2.end(), geomCompare);

      return std::equal(
          geoms1.begin(), geoms1.end(), geoms2.begin(),
          [&compareFunc, sorted](const Geometry *a, const Geometry *b) {
            return compareGeometries(*a, *b, compareFunc, sorted);
          });
    } else {
      for (size_t i = 0; i < gc1.numGeometries(); ++i) {
        if (!compareGeometries(gc1.geometryN(i), gc2.geometryN(i), compareFunc,
                               sorted))
          return false;
      }
      return true;
    }
  }
  case TYPE_SOLID: {
    const auto &s1 = g1.as<Solid>();
    const auto &s2 = g2.as<Solid>();

    if (!s1.exteriorShell().isEmpty() && !s2.exteriorShell().isEmpty()) {
      if (!compareGeometries(s1.exteriorShell(), s2.exteriorShell(),
                             compareFunc, sorted))
        return false;
    } else if (s1.exteriorShell().isEmpty() != s2.exteriorShell().isEmpty()) {
      return false;
    }

    if (s1.numInteriorShells() != s2.numInteriorShells())
      return false;

    if (sorted) {
      std::vector<const PolyhedralSurface *> shells1(s1.numInteriorShells());
      std::vector<const PolyhedralSurface *> shells2(s2.numInteriorShells());

      for (size_t i = 0; i < s1.numInteriorShells(); ++i) {
        shells1[i] = &s1.interiorShellN(i);
        shells2[i] = &s2.interiorShellN(i);
      }

      auto shellCompare = [](const PolyhedralSurface *a,
                             const PolyhedralSurface *b) {
        return findLowestPoint(a) < findLowestPoint(b);
      };

      std::sort(shells1.begin(), shells1.end(), shellCompare);
      std::sort(shells2.begin(), shells2.end(), shellCompare);

      return std::equal(shells1.begin(), shells1.end(), shells2.begin(),
                        [&compareFunc, sorted](const PolyhedralSurface *a,
                                               const PolyhedralSurface *b) {
                          return compareGeometries(*a, *b, compareFunc, sorted);
                        });
    } else {
      for (size_t i = 0; i < s1.numInteriorShells(); ++i) {
        if (!compareGeometries(s1.interiorShellN(i), s2.interiorShellN(i),
                               compareFunc, sorted))
          return false;
      }
      return true;
    }
  }
  default:
    BOOST_THROW_EXCEPTION(Exception(
        (boost::format("Unsupported geometry type %s") % g1.geometryType())
            .str()));
  }
}

} // anonymous namespace

bool
strictCompare(const Geometry &ga, const Geometry &gb)
{
  return compareGeometries(
      ga, gb, [](const auto &a, const auto &b) { return a == b; }, false);
}
std::unique_ptr<Geometry>
sortGeometry(const Geometry &g)
{
  switch (g.geometryTypeId()) {
  case TYPE_POINT:
  case TYPE_LINESTRING:
    return std::unique_ptr<Geometry>(g.clone());
  case TYPE_POLYGON: {
    auto sortedPoly =
        std::unique_ptr<Polygon>(static_cast<Polygon *>(g.clone()));
    std::vector<std::unique_ptr<LineString>> rings;
    for (size_t i = 0; i < sortedPoly->numInteriorRings(); ++i) {
      rings.push_back(std::unique_ptr<LineString>(
          static_cast<LineString *>(sortedPoly->interiorRingN(i).clone())));
    }
    // Tri basé sur les critères appropriés pour les anneaux intérieurs
    std::sort(rings.begin(), rings.end(),
              [](const std::unique_ptr<LineString> &a,
                 const std::unique_ptr<LineString> &b) {
                return findLowestPoint(a.get()) < findLowestPoint(b.get());
              });
    sortedPoly->removeInteriorRings();
    for (const auto &ring : rings) {
      sortedPoly->addRing(*ring);
    }
    return sortedPoly;
  }
  case TYPE_TRIANGLE: {
    auto sortedTri =
        std::unique_ptr<Triangle>(static_cast<Triangle *>(g.clone()));
    std::vector<Point> points = {sortedTri->vertex(0), sortedTri->vertex(1),
                                 sortedTri->vertex(2)};
    std::sort(points.begin(), points.end(), [](const Point &a, const Point &b) {
      return a.x() < b.x() ||
             (a.x() == b.x() &&
              (a.y() < b.y() || (a.y() == b.y() && a.z() < b.z())));
    });
    *sortedTri = Triangle(points[0], points[1], points[2]);
    return sortedTri;
  }
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_MULTIPOLYGON:
  case TYPE_GEOMETRYCOLLECTION:
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE: {
    auto sortedGC = std::unique_ptr<GeometryCollection>(
        static_cast<GeometryCollection *>(g.clone()));
    std::vector<std::unique_ptr<Geometry>> geoms;
    for (size_t i = 0; i < sortedGC->numGeometries(); ++i) {
      std::unique_ptr<Geometry> sortedGeom =
          sortGeometry(sortedGC->geometryN(i)); // Appel de sortGeometry
      geoms.push_back(std::move(sortedGeom));
    }
    // Tri des géométries au sein de la collection
    std::sort(geoms.begin(), geoms.end(),
              [](const std::unique_ptr<Geometry> &a,
                 const std::unique_ptr<Geometry> &b) {
                return findLowestPoint(a.get()) < findLowestPoint(b.get());
              });
    GeometryCollection newGC;
    for (const auto &geom : geoms) {
      newGC.addGeometry(geom->clone());
    }
    return std::unique_ptr<Geometry>(newGC.clone());
  }
  case TYPE_SOLID: {
    auto sortedSolid = std::unique_ptr<Solid>(static_cast<Solid *>(g.clone()));
    std::vector<std::unique_ptr<PolyhedralSurface>> shells;
    for (size_t i = 0; i < sortedSolid->numInteriorShells(); ++i) {
      std::unique_ptr<Geometry> sortedShell =
          sortGeometry(sortedSolid->interiorShellN(i)); // Appel de sortGeometry
      shells.push_back(std::unique_ptr<PolyhedralSurface>(
          static_cast<PolyhedralSurface *>(sortedShell.release())));
    }
    // Tri des "shells" basés sur un critère géométrique (lowest point)
    std::sort(shells.begin(), shells.end(),
              [](const std::unique_ptr<PolyhedralSurface> &a,
                 const std::unique_ptr<PolyhedralSurface> &b) {
                return findLowestPoint(a.get()) < findLowestPoint(b.get());
              });
    sortedSolid->removeInteriorShells();
    for (const auto &shell : shells) {
      sortedSolid->addInteriorShell(*shell);
    }
    return sortedSolid;
  }
  default:
    BOOST_THROW_EXCEPTION(Exception(
        (boost::format("Unsupported geometry type %s") % g.geometryType())
            .str()));
  }
}

bool
sortedCompare(const Geometry &ga, const Geometry &gb, bool useFuzzy,
              double epsilon, FuzzyCompareMethod method)
{
  // Lambda pour trier les géométries en fonction de critères plus complexes
  std::function<std::unique_ptr<Geometry>(const Geometry &)> sortGeometry =
      [&](const Geometry &g) -> std::unique_ptr<Geometry> {
    switch (g.geometryTypeId()) {
    case TYPE_POINT:
    case TYPE_LINESTRING:
      return std::unique_ptr<Geometry>(g.clone());
    case TYPE_POLYGON: {
      auto sortedPoly =
          std::unique_ptr<Polygon>(static_cast<Polygon *>(g.clone()));
      std::vector<std::unique_ptr<LineString>> rings;
      for (size_t i = 0; i < sortedPoly->numInteriorRings(); ++i) {
        rings.push_back(std::unique_ptr<LineString>(
            static_cast<LineString *>(sortedPoly->interiorRingN(i).clone())));
      }
      // Tri des anneaux en fonction de la distance cumulative des points pour
      // éviter les erreurs de tri
      std::sort(rings.begin(), rings.end(),
                [](const std::unique_ptr<LineString> &a,
                   const std::unique_ptr<LineString> &b) {
                  return findLowestPoint(a.get()) < findLowestPoint(b.get());
                });
      sortedPoly->removeInteriorRings();
      for (const auto &ring : rings) {
        sortedPoly->addRing(*ring);
      }
      return sortedPoly;
    }
    case TYPE_TRIANGLE: {
      auto sortedTri =
          std::unique_ptr<Triangle>(static_cast<Triangle *>(g.clone()));
      std::vector<Point> points = {sortedTri->vertex(0), sortedTri->vertex(1),
                                   sortedTri->vertex(2)};
      std::sort(points.begin(), points.end(),
                [](const Point &a, const Point &b) {
                  return a.x() < b.x() ||
                         (a.x() == b.x() &&
                          (a.y() < b.y() || (a.y() == b.y() && a.z() < b.z())));
                });
      *sortedTri = Triangle(points[0], points[1], points[2]);
      return sortedTri;
    }
    case TYPE_MULTIPOINT:
    case TYPE_MULTILINESTRING:
    case TYPE_MULTIPOLYGON:
    case TYPE_GEOMETRYCOLLECTION:
    case TYPE_TRIANGULATEDSURFACE:
    case TYPE_POLYHEDRALSURFACE: {
      auto sortedGC = std::unique_ptr<GeometryCollection>(
          static_cast<GeometryCollection *>(g.clone()));
      std::vector<std::unique_ptr<Geometry>> geoms;
      for (size_t i = 0; i < sortedGC->numGeometries(); ++i) {
        geoms.push_back(sortGeometry(sortedGC->geometryN(i)));
      }
      std::sort(geoms.begin(), geoms.end(),
                [](const std::unique_ptr<Geometry> &a,
                   const std::unique_ptr<Geometry> &b) {
                  return findLowestPoint(a.get()) < findLowestPoint(b.get());
                });
      GeometryCollection newGC;
      for (const auto &geom : geoms) {
        newGC.addGeometry(geom->clone());
      }
      return std::unique_ptr<Geometry>(newGC.clone());
    }
    case TYPE_SOLID: {
      auto sortedSolid =
          std::unique_ptr<Solid>(static_cast<Solid *>(g.clone()));
      std::vector<std::unique_ptr<PolyhedralSurface>> shells;
      for (size_t i = 0; i < sortedSolid->numInteriorShells(); ++i) {
        shells.push_back(
            std::unique_ptr<PolyhedralSurface>(static_cast<PolyhedralSurface *>(
                sortGeometry(sortedSolid->interiorShellN(i)).release())));
      }
      // Tri des shells en fonction de la géométrie
      std::sort(shells.begin(), shells.end(),
                [](const std::unique_ptr<PolyhedralSurface> &a,
                   const std::unique_ptr<PolyhedralSurface> &b) {
                  return findLowestPoint(a.get()) < findLowestPoint(b.get());
                });
      sortedSolid->removeInteriorShells();
      for (const auto &shell : shells) {
        sortedSolid->addInteriorShell(*shell);
      }
      return sortedSolid;
    }
    default:
      BOOST_THROW_EXCEPTION(Exception(
          (boost::format("Unsupported geometry type %s") % g.geometryType())
              .str()));
    }
  };

  // Tri des géométries avant la comparaison
  std::unique_ptr<Geometry> sortedGa = sortGeometry(ga);
  std::unique_ptr<Geometry> sortedGb = sortGeometry(gb);

  // Comparaison avec ou sans fuzzy
  auto compareFunc = [useFuzzy, epsilon, method](const auto &a, const auto &b) {
    if (useFuzzy) {
      if constexpr (std::is_same_v<std::decay_t<decltype(a)>, Point>) {
        return fuzzyComparePoints(a, b, epsilon, method);
      } else if constexpr (std::is_same_v<std::decay_t<decltype(a)>,
                                          LineString>) {
        return fuzzyCompareLineStrings(a, b, epsilon, method);
      } else if constexpr (std::is_same_v<std::decay_t<decltype(a)>, Polygon>) {
        return fuzzyComparePolygons(a, b, epsilon, method);
      } else if constexpr (std::is_same_v<std::decay_t<decltype(a)>,
                                          Triangle>) {
        return fuzzyCompareTriangles(a, b, epsilon, method);
      } else {
        return a == b;
      }
    } else {
      return a == b;
    }
  };

  return compareGeometries(*sortedGa, *sortedGb, compareFunc, true);
}

bool
fuzzyCompare(const Geometry &ga, const Geometry &gb, double epsilon,
             FuzzyCompareMethod method)
{
  return compareGeometries(
      ga, gb,
      [epsilon, method](const auto &a, const auto &b) {
        if constexpr (std::is_same_v<std::decay_t<decltype(a)>, Point>) {
          return fuzzyComparePoints(a, b, epsilon, method);
        } else if constexpr (std::is_same_v<std::decay_t<decltype(a)>,
                                            LineString>) {
          return fuzzyCompareLineStrings(a, b, epsilon, method);
        } else if constexpr (std::is_same_v<std::decay_t<decltype(a)>,
                                            Polygon>) {
          return fuzzyComparePolygons(a, b, epsilon, method);
        } else if constexpr (std::is_same_v<std::decay_t<decltype(a)>,
                                            Triangle>) {
          return fuzzyCompareTriangles(a, b, epsilon, method);
        } else {
          return a == b;
        }
      },
      false);
}

} // namespace SFCGAL::algorithm::compare
