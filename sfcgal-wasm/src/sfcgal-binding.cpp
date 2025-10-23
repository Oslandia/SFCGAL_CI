#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <string>
#include <memory>
#include <sstream>

// SFCGAL headers
#include <SFCGAL/Geometry.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/TriangulatedSurface.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/io/wkt.h>
#include <SFCGAL/version.h>

// SFCGAL algorithms
#include <SFCGAL/algorithm/area.h>
#include <SFCGAL/algorithm/centroid.h>
#include <SFCGAL/algorithm/distance.h>
#include <SFCGAL/algorithm/distance3d.h>
#include <SFCGAL/algorithm/length.h>
#include <SFCGAL/algorithm/intersection.h>
#include <SFCGAL/algorithm/union.h>
#include <SFCGAL/algorithm/difference.h>
#include <SFCGAL/algorithm/offset.h>
#include <SFCGAL/algorithm/extrude.h>
#include <SFCGAL/algorithm/convexHull.h>
#include <SFCGAL/algorithm/isValid.h>
#include <SFCGAL/algorithm/translate.h>
#include <SFCGAL/algorithm/rotate.h>
#include <SFCGAL/algorithm/volume.h>
#include <SFCGAL/detail/transform/AffineTransform3.h>

using namespace emscripten;

class SFCGALWrapper {
public:
    SFCGALWrapper() {}

    void initialize() {
        // SFCGAL initialization if needed
    }

    std::string version() const {
        return SFCGAL_VERSION;
    }

    // Validation
    bool isValid(const std::string& wkt) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));
            if (!geom || geom->isEmpty()) return false;
            return SFCGAL::algorithm::isValid(*geom);
        } catch (...) {
            return false;
        }
    }

    // Get geometry information for debugging
    val getGeometryInfo(const std::string& wkt) const {
        val info = val::object();
        try {
            std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));
            if (!geom) {
                info.set("error", std::string("Failed to parse WKT"));
                return info;
            }

            info.set("isEmpty", geom->isEmpty());
            info.set("is3D", geom->is3D());
            info.set("dimension", static_cast<int>(geom->dimension()));
            info.set("coordinateDimension", static_cast<int>(geom->coordinateDimension()));

            std::string geomType;
            switch (geom->geometryTypeId()) {
                case SFCGAL::TYPE_POINT: geomType = "POINT"; break;
                case SFCGAL::TYPE_LINESTRING: geomType = "LINESTRING"; break;
                case SFCGAL::TYPE_POLYGON: geomType = "POLYGON"; break;
                case SFCGAL::TYPE_MULTIPOINT: geomType = "MULTIPOINT"; break;
                case SFCGAL::TYPE_MULTILINESTRING: geomType = "MULTILINESTRING"; break;
                case SFCGAL::TYPE_MULTIPOLYGON: geomType = "MULTIPOLYGON"; break;
                case SFCGAL::TYPE_GEOMETRYCOLLECTION: geomType = "GEOMETRYCOLLECTION"; break;
                case SFCGAL::TYPE_POLYHEDRALSURFACE: geomType = "POLYHEDRALSURFACE"; break;
                case SFCGAL::TYPE_TRIANGULATEDSURFACE: geomType = "TRIANGULATEDSURFACE"; break;
                case SFCGAL::TYPE_SOLID: geomType = "SOLID"; break;
                default: geomType = "UNKNOWN";
            }
            info.set("geometryType", geomType);
            info.set("isValid", SFCGAL::algorithm::isValid(*geom));

            // Count elements for POLYHEDRALSURFACE
            if (geom->geometryTypeId() == SFCGAL::TYPE_POLYHEDRALSURFACE) {
                const SFCGAL::PolyhedralSurface& surf = geom->as<SFCGAL::PolyhedralSurface>();
                info.set("numPolygons", static_cast<int>(surf.numPolygons()));
            }

        } catch (const std::exception& e) {
            info.set("error", std::string(e.what()));
        } catch (...) {
            info.set("error", std::string("Unknown error"));
        }
        return info;
    }

    // Convert POLYHEDRALSURFACE to SOLID
    std::string toSolid(const std::string& wkt) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));
            if (!geom) {
                return "ERROR: Failed to parse WKT";
            }

            // If already a SOLID, return as-is
            if (geom->geometryTypeId() == SFCGAL::TYPE_SOLID) {
                return geom->asText(10);
            }

            // If POLYHEDRALSURFACE, wrap it in a SOLID
            if (geom->geometryTypeId() == SFCGAL::TYPE_POLYHEDRALSURFACE) {
                SFCGAL::PolyhedralSurface* surf = static_cast<SFCGAL::PolyhedralSurface*>(geom->clone());
                SFCGAL::Solid solid(surf);
                return solid.asText(10);
            }

            return "ERROR: Can only convert POLYHEDRALSURFACE to SOLID";
        } catch (const std::exception& e) {
            std::string msg = "ERROR: toSolid failed: ";
            msg += e.what();
            return msg;
        } catch (...) {
            return "ERROR: Unknown error in toSolid";
        }
    }

    // Metrics
    double area(const std::string& wkt) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));
            if (!geom) return 0.0;
            if (!SFCGAL::algorithm::isValid(*geom)) {
                std::cerr << "Warning: Invalid geometry for area calculation" << std::endl;
                return 0.0;
            }
            return CGAL::to_double(SFCGAL::algorithm::area(*geom));
        } catch (...) {
            return 0.0;
        }
    }

    double area3D(const std::string& wkt) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));
            if (!geom) return 0.0;
            return CGAL::to_double(SFCGAL::algorithm::area3D(*geom));
        } catch (...) {
            return 0.0;
        }
    }

    double volume(const std::string& wkt) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));
            if (!geom) return 0.0;
            return CGAL::to_double(SFCGAL::algorithm::volume(*geom));
        } catch (...) {
            return 0.0;
        }
    }

    double length(const std::string& wkt) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));
            if (!geom) return 0.0;
            return CGAL::to_double(SFCGAL::algorithm::length(*geom));
        } catch (...) {
            return 0.0;
        }
    }

    double length3D(const std::string& wkt) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));
            if (!geom) return 0.0;
            return CGAL::to_double(SFCGAL::algorithm::length3D(*geom));
        } catch (...) {
            return 0.0;
        }
    }

    double perimeter(const std::string& wkt) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));
            if (!geom) return 0.0;

            // For polygons, calculate perimeter of exterior ring
            if (geom->geometryTypeId() == SFCGAL::TYPE_POLYGON) {
                const SFCGAL::Polygon& poly = geom->as<SFCGAL::Polygon>();
                return CGAL::to_double(SFCGAL::algorithm::length(poly.exteriorRing()));
            }

            // For other geometries, use length
            return CGAL::to_double(SFCGAL::algorithm::length(*geom));
        } catch (...) {
            return 0.0;
        }
    }

    double perimeter3D(const std::string& wkt) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));
            if (!geom) return 0.0;

            // For polygons, calculate 3D perimeter of exterior ring
            if (geom->geometryTypeId() == SFCGAL::TYPE_POLYGON) {
                const SFCGAL::Polygon& poly = geom->as<SFCGAL::Polygon>();
                return CGAL::to_double(SFCGAL::algorithm::length3D(poly.exteriorRing()));
            }

            // For other geometries, use length3D
            return CGAL::to_double(SFCGAL::algorithm::length3D(*geom));
        } catch (...) {
            return 0.0;
        }
    }

    double distance(const std::string& wkt1, const std::string& wkt2) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom1(SFCGAL::io::readWkt(wkt1));
            std::unique_ptr<SFCGAL::Geometry> geom2(SFCGAL::io::readWkt(wkt2));
            if (!geom1 || !geom2) return 0.0;
            return CGAL::to_double(SFCGAL::algorithm::distance(*geom1, *geom2));
        } catch (...) {
            return 0.0;
        }
    }

    double distance3D(const std::string& wkt1, const std::string& wkt2) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom1(SFCGAL::io::readWkt(wkt1));
            std::unique_ptr<SFCGAL::Geometry> geom2(SFCGAL::io::readWkt(wkt2));
            if (!geom1 || !geom2) return 0.0;
            return CGAL::to_double(SFCGAL::algorithm::distance3D(*geom1, *geom2));
        } catch (...) {
            return 0.0;
        }
    }

    // Geometric operations
    std::string centroid(const std::string& wkt) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));
            if (!geom) return "POINT EMPTY";

            std::unique_ptr<SFCGAL::Geometry> result(SFCGAL::algorithm::centroid(*geom));
            return result ? result->asText(10) : "POINT EMPTY";

        } catch (...) {
            return "POINT EMPTY";
        }
    }

    std::string intersection(const std::string& wkt1, const std::string& wkt2) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom1(SFCGAL::io::readWkt(wkt1));
            std::unique_ptr<SFCGAL::Geometry> geom2(SFCGAL::io::readWkt(wkt2));
            if (!geom1 || !geom2) return "GEOMETRYCOLLECTION EMPTY";

            std::unique_ptr<SFCGAL::Geometry> result(SFCGAL::algorithm::intersection(*geom1, *geom2));
            return result ? result->asText(10) : "GEOMETRYCOLLECTION EMPTY";
        } catch (...) {
            return "GEOMETRYCOLLECTION EMPTY";
        }
    }

    std::string intersection3D(const std::string& wkt1, const std::string& wkt2) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom1(SFCGAL::io::readWkt(wkt1));
            std::unique_ptr<SFCGAL::Geometry> geom2(SFCGAL::io::readWkt(wkt2));
            if (!geom1 || !geom2) {
                return "ERROR: Failed to parse WKT for intersection3D";
            }

            if (!SFCGAL::algorithm::isValid(*geom1)) {
                return "ERROR: Invalid geometry 1 for intersection3D. Check geometry validity with isValid().";
            }
            if (!SFCGAL::algorithm::isValid(*geom2)) {
                return "ERROR: Invalid geometry 2 for intersection3D. Check geometry validity with isValid().";
            }

            std::unique_ptr<SFCGAL::Geometry> result(SFCGAL::algorithm::intersection3D(*geom1, *geom2));
            return result ? result->asText(10) : "GEOMETRYCOLLECTION EMPTY";
        } catch (const std::exception& e) {
            std::string msg = "ERROR: intersection3D failed: ";
            msg += e.what();
            return msg;
        } catch (...) {
            return "ERROR: Unknown error in intersection3D";
        }
    }

    std::string union_(const std::string& wkt1, const std::string& wkt2) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom1(SFCGAL::io::readWkt(wkt1));
            std::unique_ptr<SFCGAL::Geometry> geom2(SFCGAL::io::readWkt(wkt2));
            if (!geom1 || !geom2) return "GEOMETRYCOLLECTION EMPTY";

            std::unique_ptr<SFCGAL::Geometry> result(SFCGAL::algorithm::union_(*geom1, *geom2));
            return result ? result->asText(10) : "GEOMETRYCOLLECTION EMPTY";
        } catch (...) {
            return "GEOMETRYCOLLECTION EMPTY";
        }
    }

    std::string union3D(const std::string& wkt1, const std::string& wkt2) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom1(SFCGAL::io::readWkt(wkt1));
            std::unique_ptr<SFCGAL::Geometry> geom2(SFCGAL::io::readWkt(wkt2));
            if (!geom1 || !geom2) {
                return "ERROR: Failed to parse WKT for union3D";
            }

            // Check if empty
            if (geom1->isEmpty() || geom2->isEmpty()) {
                return "ERROR: One or both geometries are empty";
            }

            // Check validity
            if (!SFCGAL::algorithm::isValid(*geom1)) {
                return "ERROR: Invalid geometry 1 for union3D. Check geometry validity with isValid().";
            }
            if (!SFCGAL::algorithm::isValid(*geom2)) {
                return "ERROR: Invalid geometry 2 for union3D. Check geometry validity with isValid().";
            }

            // Perform union3D
            std::unique_ptr<SFCGAL::Geometry> result(SFCGAL::algorithm::union3D(*geom1, *geom2));
            return result ? result->asText(10) : "GEOMETRYCOLLECTION EMPTY";
        } catch (const std::bad_alloc& e) {
            return "ERROR: Memory allocation failed. The geometries may be too complex. Try simplifying them or converting POLYHEDRALSURFACE to SOLID with toSolid().";
        } catch (const std::exception& e) {
            std::string msg = "ERROR: union3D failed: ";
            msg += e.what();
            return msg;
        } catch (...) {
            return "ERROR: Unknown error in union3D";
        }
    }

    std::string difference(const std::string& wkt1, const std::string& wkt2) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom1(SFCGAL::io::readWkt(wkt1));
            std::unique_ptr<SFCGAL::Geometry> geom2(SFCGAL::io::readWkt(wkt2));
            if (!geom1 || !geom2) return "GEOMETRYCOLLECTION EMPTY";

            std::unique_ptr<SFCGAL::Geometry> result(SFCGAL::algorithm::difference(*geom1, *geom2));
            return result ? result->asText(10) : "GEOMETRYCOLLECTION EMPTY";
        } catch (...) {
            return "GEOMETRYCOLLECTION EMPTY";
        }
    }

    std::string difference3D(const std::string& wkt1, const std::string& wkt2) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom1(SFCGAL::io::readWkt(wkt1));
            std::unique_ptr<SFCGAL::Geometry> geom2(SFCGAL::io::readWkt(wkt2));
            if (!geom1 || !geom2) {
                return "ERROR: Failed to parse WKT for difference3D";
            }

            if (!SFCGAL::algorithm::isValid(*geom1)) {
                return "ERROR: Invalid geometry 1 for difference3D. Check geometry validity with isValid().";
            }
            if (!SFCGAL::algorithm::isValid(*geom2)) {
                return "ERROR: Invalid geometry 2 for difference3D. Check geometry validity with isValid().";
            }

            std::unique_ptr<SFCGAL::Geometry> result(SFCGAL::algorithm::difference3D(*geom1, *geom2));
            return result ? result->asText(10) : "GEOMETRYCOLLECTION EMPTY";
        } catch (const std::exception& e) {
            std::string msg = "ERROR: difference3D failed: ";
            msg += e.what();
            return msg;
        } catch (...) {
            return "ERROR: Unknown error in difference3D";
        }
    }

    std::string convexHull(const std::string& wkt) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));
            if (!geom) return "GEOMETRYCOLLECTION EMPTY";

            std::unique_ptr<SFCGAL::Geometry> result(SFCGAL::algorithm::convexHull(*geom));
            return result ? result->asText(10) : "GEOMETRYCOLLECTION EMPTY";
        } catch (...) {
            return "GEOMETRYCOLLECTION EMPTY";
        }
    }

    std::string buffer(const std::string& wkt, double distance) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));
            if (!geom) return "GEOMETRYCOLLECTION EMPTY";

            std::unique_ptr<SFCGAL::Geometry> result(SFCGAL::algorithm::offset(*geom, distance));
            return result ? result->asText(10) : "GEOMETRYCOLLECTION EMPTY";
        } catch (...) {
            return "GEOMETRYCOLLECTION EMPTY";
        }
    }

    // 3D operations
    std::string extrude(const std::string& wkt, double height) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));
            if (!geom) return "GEOMETRYCOLLECTION EMPTY";

            std::unique_ptr<SFCGAL::Geometry> result(SFCGAL::algorithm::extrude(*geom, 0, 0, height));
            // Use higher precision (10 decimals) for accurate 3D coordinates
            return result ? result->asText(10) : "GEOMETRYCOLLECTION EMPTY";
        } catch (...) {
            return "GEOMETRYCOLLECTION EMPTY";
        }
    }

    val extrudeDetailed(const std::string& wkt, double height, const val& options) {
        val result = val::object();

        try {
            std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));
            if (!geom) {
                result.set("wkt", std::string("GEOMETRYCOLLECTION EMPTY"));
                return result;
            }

            std::unique_ptr<SFCGAL::Geometry> extruded(SFCGAL::algorithm::extrude(*geom, 0, 0, height));

            result.set("wkt", extruded ? extruded->asText(10) : std::string("GEOMETRYCOLLECTION EMPTY"));
            result.set("baseArea", CGAL::to_double(SFCGAL::algorithm::area(*geom)));
            result.set("volume", CGAL::to_double(SFCGAL::algorithm::volume(*extruded)));
            result.set("perimeter", this->length(wkt));

            val stats = val::object();
            if (extruded && extruded->geometryTypeId() == SFCGAL::TYPE_SOLID) {
                const SFCGAL::Solid& solid = extruded->as<SFCGAL::Solid>();
                stats.set("numFaces", solid.numShells() > 0 ? solid.shellN(0).numPolygons() : 0);
            } else if (extruded && extruded->geometryTypeId() == SFCGAL::TYPE_POLYHEDRALSURFACE) {
                const SFCGAL::PolyhedralSurface& surf = extruded->as<SFCGAL::PolyhedralSurface>();
                stats.set("numFaces", surf.numPolygons());
            }
            result.set("stats", stats);

        } catch (...) {
            result.set("wkt", std::string("GEOMETRYCOLLECTION EMPTY"));
        }

        return result;
    }

    // Transformations
    std::string translate(const std::string& wkt, double dx, double dy, double dz) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));
            if (!geom) return wkt;

            SFCGAL::algorithm::translate(*geom, dx, dy, dz);
            return geom->asText(10);
        } catch (...) {
            return wkt;
        }
    }

    std::string rotate(const std::string& wkt, double angle) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));
            if (!geom) return wkt;

            SFCGAL::algorithm::rotate(*geom, angle);
            return geom->asText(10);
        } catch (...) {
            return wkt;
        }
    }

    std::string scale(const std::string& wkt, double sx, double sy, double sz) const {
        try {
            std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));
            if (!geom) return wkt;

            SFCGAL::transform::AffineTransform3 transform(
                CGAL::Aff_transformation_3<SFCGAL::Kernel>(
                    sx, 0, 0, 0,
                    0, sy, 0, 0,
                    0, 0, sz, 1
                )
            );
            geom->accept(transform);
            return geom->asText(10);
        } catch (...) {
            return wkt;
        }
    }
};

EMSCRIPTEN_BINDINGS(sfcgal_module) {
    class_<SFCGALWrapper>("SFCGAL")
        .constructor<>()
        .function("initialize", &SFCGALWrapper::initialize)
        .function("version", &SFCGALWrapper::version)
        .function("isValid", &SFCGALWrapper::isValid)
        .function("getGeometryInfo", &SFCGALWrapper::getGeometryInfo)
        .function("toSolid", &SFCGALWrapper::toSolid)
        .function("area", &SFCGALWrapper::area)
        .function("area3D", &SFCGALWrapper::area3D)
        .function("length", &SFCGALWrapper::length)
        .function("length3D", &SFCGALWrapper::length3D)
        .function("perimeter", &SFCGALWrapper::perimeter)
        .function("perimeter3D", &SFCGALWrapper::perimeter3D)
        .function("distance", &SFCGALWrapper::distance)
        .function("distance3D", &SFCGALWrapper::distance3D)
        .function("centroid", &SFCGALWrapper::centroid)
        .function("intersection", &SFCGALWrapper::intersection)
        .function("intersection3D", &SFCGALWrapper::intersection3D)
        .function("union", &SFCGALWrapper::union_)
        .function("union3D", &SFCGALWrapper::union3D)
        .function("difference", &SFCGALWrapper::difference)
        .function("difference3D", &SFCGALWrapper::difference3D)
        .function("convexHull", &SFCGALWrapper::convexHull)
        .function("buffer", &SFCGALWrapper::buffer)
        .function("extrude", &SFCGALWrapper::extrude)
        .function("extrudeDetailed", &SFCGALWrapper::extrudeDetailed)
        .function("translate", &SFCGALWrapper::translate)
        .function("rotate", &SFCGALWrapper::rotate)
        .function("scale", &SFCGALWrapper::scale)
        .function("volume", &SFCGALWrapper::volume);
}
