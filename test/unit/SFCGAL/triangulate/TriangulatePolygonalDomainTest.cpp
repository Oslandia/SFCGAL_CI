/**
 *   SFCGAL
 *
 *   Copyright (C) 2012-2013 Oslandia <infos@oslandia.com>
 *   Copyright (C) 2012-2013 IGN (http://www.ign.fr)
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Library General Public
 *   License as published by the Free Software Foundation; either
 *   version 2 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Library General Public License for more details.

 *   You should have received a copy of the GNU Library General Public
 *   License along with this library; if not, see
 <http://www.gnu.org/licenses/>.
 */
#include <boost/test/unit_test.hpp>

#include <SFCGAL/Exception.h>
#include <SFCGAL/TriangulatedSurface.h>
#include <SFCGAL/io/wkt.h>
#include <SFCGAL/triangulate/triangulatePolygonalDomain.h>

using namespace boost::unit_test;
using namespace SFCGAL;
using namespace SFCGAL::triangulate;

BOOST_AUTO_TEST_SUITE(SFCGAL_triangulate_TriangulatePolygonalDomain)

BOOST_AUTO_TEST_CASE(testPolygon)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))"));
  TriangulatedSurface triangulation = constrainedPolygonalDomain(
      g->as<Polygon>(), MultiLineString(), MultiPoint());
  BOOST_CHECK_EQUAL(triangulation.numTriangles(), 2U);
  std::string expectedWKT("TIN(((0 1,1 0,1 1,0 1)),((0 1,0 0,1 0,0 1)))");
  BOOST_CHECK_EQUAL(triangulation.asText(0), expectedWKT);
}

BOOST_AUTO_TEST_CASE(testPolygon_line_constrained)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))"));
  std::unique_ptr<Geometry> l(
      io::readWkt("MULTILINESTRING((0.2 0.2,0.2 0.8,0.8 0.8,0.8 0.2,0.20.2))"));
  TriangulatedSurface triangulation = constrainedPolygonalDomain(
      g->as<Polygon>(), l->as<MultiLineString>(), MultiPoint());
  BOOST_CHECK_EQUAL(triangulation.numTriangles(), 8U);
  std::string expectedWKT(
      "TIN(((0.8 0.2,0.2 0.2,1.0 0.0,0.8 0.2)),((0.2 0.2,0.0 0.0,1.0 0.0,0.2 "
      "0.2)),((1.0 1.0,0.8 0.8,0.8 0.2,1.0 1.0)),((0.0 1.0,0.0 0.0,0.2 0.2,0.0 "
      "1.0)),((0.0 1.0,0.2 0.8,1.0 1.0,0.0 1.0)),((0.0 1.0,0.2 0.2,0.2 0.8,0.0 "
      "1.0)),((0.2 0.8,0.8 0.8,1.0 1.0,0.2 0.8)),((1.0 1.0,0.8 0.2,1.0 0.0,1.0 "
      "1.0)))");
  BOOST_CHECK_EQUAL(triangulation.asText(1), expectedWKT);
}

BOOST_AUTO_TEST_SUITE_END()
