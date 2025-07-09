// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <memory>
#include <string>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/detail/io/Serialization.h"
#include "SFCGAL/io/ewkt.h"
#include "SFCGAL/io/wkt.h"

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_io_SerializationTest)

BOOST_AUTO_TEST_CASE(textTest)
{
  Coordinate const pt1(2.3, 8.9);
  Coordinate       rpt1;
  // 2/3
  CGAL::Gmpq const q1(2, 3);
  CGAL::Gmpq       rq1;

  std::ostringstream            ostr;
  boost::archive::text_oarchive arc(ostr);
  arc << pt1;
  arc << q1;

  std::string const str = ostr.str();

  std::istringstream            istr(str);
  boost::archive::text_iarchive iarc(istr);
  iarc >> rpt1;
  iarc >> rq1;

  BOOST_CHECK(pt1 == rpt1);
  BOOST_CHECK(q1 == rq1);
}

BOOST_AUTO_TEST_CASE(binaryTest)
{
  Coordinate const pt(2.3, 4.5);
  Coordinate       rpt;

  std::ostringstream   ostr;
  io::BinarySerializer arc(ostr);
  arc << pt;

  std::string const str = ostr.str();

  std::istringstream     istr(str);
  io::BinaryUnserializer iarc(istr);
  iarc >> rpt;

  BOOST_CHECK(pt == rpt);
}

BOOST_AUTO_TEST_CASE(geometryTest)
{
  std::unique_ptr<Geometry> g1 = io::readWkt("POINT ( 3.4 4.5 5.6 )");
  std::unique_ptr<Geometry> g2 =
      io::readWkt("LINESTRING ( 3.4 4.5 5.6, 5 6 8 )");
  std::unique_ptr<Geometry> g3 =
      io::readWkt("TRIANGLE (( 0 0 0, 3.4 5.6 6.7,2 3 4, 0 0 0))");
  std::unique_ptr<Geometry> g4 =
      io::readWkt("POLYGON (( 0 0 0, 1 1 1, 3.4 5.6 6.7,2 3 4, 0 0 0))");
  std::unique_ptr<Geometry> g5 =
      io::readWkt("TIN ((( 0 0 0, 3.4 5.6 6.7,2 3 4, 0 0 0)), ((0 0 0,0 1 0,1 "
                  "1 0,0 0 0)))");
  std::unique_ptr<Geometry> g6 = io::readWkt(
      "POLYHEDRALSURFACE ((( 0 0 0, 3.4 5.6 6.7,2 3 4, 0 0 0)), ((0 "
      "0 0,0 1 0,1 1 0,1 0 0,0 0 0)))");
  std::unique_ptr<Geometry> g7 =
      io::readWkt("SOLID (((( 0 0 0, 3.4 5.6 6.7,2 3 4, 0 0 0)), ((0 0 0,0 1 "
                  "0,1 1 0,1 0 0,0 0 0))))");
  std::unique_ptr<Geometry> g8 = io::readWkt("MULTIPOINT (( 3.4 4.5 5.6 ))");
  std::unique_ptr<Geometry> g9 =
      io::readWkt("MULTILINESTRING (( 3.4 4.5 5.6, 5 6 8 ))");
  std::unique_ptr<Geometry> g10 =
      io::readWkt("MULTIPOLYGON ((( 0 0 0, 1 1 1, 3.4 5.6 6.7,2 3 4, 0 0 0)))");
  std::unique_ptr<Geometry> g11 =
      io::readWkt("MULTISOLID ((((( 0 0 0, 3.4 5.6 6.7,2 3 4, 0 0 0)), ((0 0 "
                  "0,0 1 0,1 1 0,1 0 0,0 0 0)))))");

  BOOST_CHECK(io::readBinaryGeometry(io::writeBinaryGeometry(*g1))->asText() ==
              g1->asText());
  BOOST_CHECK(io::readBinaryGeometry(io::writeBinaryGeometry(*g2))->asText() ==
              g2->asText());
  BOOST_CHECK(io::readBinaryGeometry(io::writeBinaryGeometry(*g3))->asText() ==
              g3->asText());
  BOOST_CHECK(io::readBinaryGeometry(io::writeBinaryGeometry(*g4))->asText() ==
              g4->asText());
  BOOST_CHECK(io::readBinaryGeometry(io::writeBinaryGeometry(*g5))->asText() ==
              g5->asText());
  BOOST_CHECK(io::readBinaryGeometry(io::writeBinaryGeometry(*g6))->asText() ==
              g6->asText());
  BOOST_CHECK(io::readBinaryGeometry(io::writeBinaryGeometry(*g7))->asText() ==
              g7->asText());
  BOOST_CHECK(io::readBinaryGeometry(io::writeBinaryGeometry(*g8))->asText() ==
              g8->asText());
  BOOST_CHECK(io::readBinaryGeometry(io::writeBinaryGeometry(*g9))->asText() ==
              g9->asText());
  BOOST_CHECK(io::readBinaryGeometry(io::writeBinaryGeometry(*g10))->asText() ==
              g10->asText());
  BOOST_CHECK(io::readBinaryGeometry(io::writeBinaryGeometry(*g11))->asText() ==
              g11->asText());
}

BOOST_AUTO_TEST_CASE(caseTest)
{
  std::unique_ptr<Geometry> g1 = io::readWkt("PoInT( 3.4 4.5 5.6 )");
  BOOST_CHECK(io::readBinaryGeometry(io::writeBinaryGeometry(*g1))->asText() ==
              g1->asText());
}

BOOST_AUTO_TEST_CASE(extBinaryTest)
{
  std::unique_ptr<Geometry> ng1 = io::readWkt(
      "POLYGON ((9.70065196277574 -2.37991360574961,3.74221071600914 "
      "5.33515858836472,-4.30443822173402 -2.37991360574961,3.74221071600914 "
      "-12.1891073728912,9.70065196277574 -2.37991360574961))");
  SFCGAL::Geometry *mg1 = ng1.release();
  SFCGAL::Geometry *ng2 = nullptr;

  std::string const    str;
  std::ostringstream   ostr(str);
  io::BinarySerializer arc(ostr);
  arc << mg1;

  std::istringstream     istr(ostr.str());
  io::BinaryUnserializer iarc(istr);
  iarc >> ng2;

  BOOST_CHECK(mg1->asText() == ng2->asText());
  delete mg1;
  delete ng2;
}

BOOST_AUTO_TEST_CASE(preparedGeometryTest)
{
  std::unique_ptr<PreparedGeometry> g1 = io::readEwkt("POINT ( 3.4 4.5 5.6 )");
  std::unique_ptr<PreparedGeometry> g2 =
      io::readEwkt("srid=0;POINT ( 3.4 4.5 5.6 )");
  std::unique_ptr<PreparedGeometry> g3 =
      io::readEwkt("srid=4326;POINT ( 3.4 4.5 5.6 )");

  BOOST_CHECK(io::readBinaryPrepared(io::writeBinaryPrepared(*g1))->asEWKT() ==
              g1->asEWKT());
  BOOST_CHECK(io::readBinaryPrepared(io::writeBinaryPrepared(*g2))->asEWKT() ==
              g2->asEWKT());
  BOOST_CHECK(io::readBinaryPrepared(io::writeBinaryPrepared(*g3))->asEWKT() ==
              g3->asEWKT());
}

BOOST_AUTO_TEST_SUITE_END()
