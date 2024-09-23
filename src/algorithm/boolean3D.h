#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <SFCGAL/Geometry.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/TriangulatedSurface.h>

namespace SFCGAL {
namespace algorithm {

class Boolean3D {
private:
  using Nef_polyhedron_3 = CGAL::Nef_polyhedron_3<Kernel>;
  using Polyhedron_3     = CGAL::Polyhedron_3<Kernel>;

  static Nef_polyhedron_3
  toNef(const Geometry &g)
  {
    if (g.is<TriangulatedSurface>() || g.is<PolyhedralSurface>()) {
      std::unique_ptr<Polyhedron_3> poly(
          g.as<PolyhedralSurface>().toPolyhedron_3<Kernel, Polyhedron_3>());
      return Nef_polyhedron_3(*poly);
    }
    throw std::invalid_argument(
        "Geometry must be TriangulatedSurface or PolyhedralSurface");
  }

  static std::unique_ptr<Geometry>
  fromNef(const Nef_polyhedron_3 &nef)
  {
    Polyhedron_3 poly;
    nef.convert_to_polyhedron(poly);
    return std::make_unique<PolyhedralSurface>(poly);
  }

public:
  enum class Operation {
    UNION,
    DIFFERENCE,
    INTERSECTION,
    SYMMETRIC_DIFFERENCE
  };

  static std::unique_ptr<Geometry>
  operation(const Geometry &gA, const Geometry &gB, Operation op)
  {
    Nef_polyhedron_3 nefA = toNef(gA);
    Nef_polyhedron_3 nefB = toNef(gB);

    switch (op) {
    case Operation::UNION:
      return fromNef(nefA + nefB);
    case Operation::DIFFERENCE:
      return fromNef(nefA - nefB);
    case Operation::INTERSECTION:
      return fromNef(nefA * nefB);
    case Operation::SYMMETRIC_DIFFERENCE:
      return fromNef(nefA ^ nefB);
    }

    throw std::invalid_argument("Invalid operation");
  }

  static std::unique_ptr<Geometry>
  union3D(const Geometry &gA, const Geometry &gB)
  {
    return operation(gA, gB, Operation::UNION);
  }

  static std::unique_ptr<Geometry>
  difference3D(const Geometry &gA, const Geometry &gB)
  {
    return operation(gA, gB, Operation::DIFFERENCE);
  }

  static std::unique_ptr<Geometry>
  intersection3D(const Geometry &gA, const Geometry &gB)
  {
    return operation(gA, gB, Operation::INTERSECTION);
  }

  static std::unique_ptr<Geometry>
  symmetricDifference3D(const Geometry &gA, const Geometry &gB)
  {
    return operation(gA, gB, Operation::SYMMETRIC_DIFFERENCE);
  }
};

} // namespace algorithm
} // namespace SFCGAL
