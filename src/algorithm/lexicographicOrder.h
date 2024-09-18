#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryCollection.h"
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

namespace SFCGAL {
namespace transform {

/**
 * Réorganise une géométrie selon un ordre lexicographique.
 * @param geometry La géométrie à réorganiser.
 * @param invertRings Option pour inverser le sens des anneaux (par défaut:
 * false).
 */
void
lexicographicOrder(Geometry &geometry, bool invertRings = false);

} // namespace transform
} // namespace SFCGAL
