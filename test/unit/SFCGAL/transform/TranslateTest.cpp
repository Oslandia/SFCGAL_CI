#include <boost/test/unit_test.hpp>

#include <SFCGAL/algorithm/translate.h>
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

using namespace SFCGAL;
using namespace SFCGAL::algorithm;

BOOST_AUTO_TEST_SUITE(SFCGAL_transform_TranslateTest)

BOOST_AUTO_TEST_CASE(testTranslatePoint2D)
{
    Point point(0.0, 0.0);
    Kernel::Vector_2 translation(1.0, 2.0);
    
    translate(point, translation);
    
    BOOST_CHECK_CLOSE(point.x(), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(point.y(), 2.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(testTranslatePoint3D)
{
    Point point(0.0, 0.0, 0.0);
    Kernel::Vector_3 translation(1.0, 2.0, 3.0);
    
    translate(point, translation);
    
    BOOST_CHECK_CLOSE(point.x(), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(point.y(), 2.0, 1e-10);
    BOOST_CHECK_CLOSE(point.z(), 3.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(testTranslateLineString2D)
{
    LineString lineString;
    lineString.addPoint(Point(0.0, 0.0));
    lineString.addPoint(Point(1.0, 1.0));
    
    translate(lineString, 2.0, 3.0, 0.0);
    
    BOOST_CHECK_CLOSE(lineString.pointN(0).x(), 2.0, 1e-10);
    BOOST_CHECK_CLOSE(lineString.pointN(0).y(), 3.0, 1e-10);
    BOOST_CHECK_CLOSE(lineString.pointN(1).x(), 3.0, 1e-10);
    BOOST_CHECK_CLOSE(lineString.pointN(1).y(), 4.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(testTranslatePolygon3D)
{
    LineString exteriorRing;
    exteriorRing.addPoint(Point(0.0, 0.0, 0.0));
    exteriorRing.addPoint(Point(1.0, 0.0, 0.0));
    exteriorRing.addPoint(Point(1.0, 1.0, 0.0));
    exteriorRing.addPoint(Point(0.0, 1.0, 0.0));
    exteriorRing.addPoint(Point(0.0, 0.0, 0.0));

    Polygon polygon(exteriorRing);
    
    translate(polygon, 1.0, 2.0, 3.0);
    
    BOOST_CHECK_CLOSE(polygon.exteriorRing().pointN(0).x(), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(polygon.exteriorRing().pointN(0).y(), 2.0, 1e-10);
    BOOST_CHECK_CLOSE(polygon.exteriorRing().pointN(0).z(), 3.0, 1e-10);
    
    BOOST_CHECK_CLOSE(polygon.exteriorRing().pointN(2).x(), 2.0, 1e-10);
    BOOST_CHECK_CLOSE(polygon.exteriorRing().pointN(2).y(), 3.0, 1e-10);
    BOOST_CHECK_CLOSE(polygon.exteriorRing().pointN(2).z(), 3.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(testTranslateGeometryCollection)
{
    GeometryCollection collection;
    collection.addGeometry(Point(0.0, 0.0));
    collection.addGeometry(LineString(Point(0.0, 0.0), Point(1.0, 1.0)));
    
    translate(collection, 1.0, 2.0, 3.0);
    
    const Point& translatedPoint = collection.geometryN(0).as<Point>();
    BOOST_CHECK_CLOSE(translatedPoint.x(), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(translatedPoint.y(), 2.0, 1e-10);
    BOOST_CHECK_CLOSE(translatedPoint.z(), 3.0, 1e-10);
    
    const LineString& translatedLineString = collection.geometryN(1).as<LineString>();
    BOOST_CHECK_CLOSE(translatedLineString.pointN(0).x(), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(translatedLineString.pointN(0).y(), 2.0, 1e-10);
    BOOST_CHECK_CLOSE(translatedLineString.pointN(0).z(), 3.0, 1e-10);
    BOOST_CHECK_CLOSE(translatedLineString.pointN(1).x(), 2.0, 1e-10);
    BOOST_CHECK_CLOSE(translatedLineString.pointN(1).y(), 3.0, 1e-10);
    BOOST_CHECK_CLOSE(translatedLineString.pointN(1).z(), 3.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(testTranslateNonFiniteValues)
{
    Point point(0.0, 0.0, 0.0);
    BOOST_CHECK_THROW(translate(point, std::numeric_limits<double>::infinity(), 0.0, 0.0), NonFiniteValueException);
    BOOST_CHECK_THROW(translate(point, 0.0, std::numeric_limits<double>::quiet_NaN(), 0.0), NonFiniteValueException);
}

BOOST_AUTO_TEST_SUITE_END()
