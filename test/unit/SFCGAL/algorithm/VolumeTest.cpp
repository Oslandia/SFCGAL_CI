// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/volume.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/io/wkt.h"
#include <CGAL/Nef_polyhedron_3.h>

#include <boost/test/unit_test.hpp>

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_VolumeTest)

BOOST_AUTO_TEST_CASE(cubeVolume)
{
  const std::unique_ptr<Geometry> s =
      io::readWkt("SOLID ((((0 0 0,0 0 1,0 1 1,0 1 0,0 0 0)),\
                                                     ((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),\
                                                     ((0 0 0,1 0 0,1 0 1,0 0 1,0 0 0)),\
                                                     ((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0)),\
                                                     ((0 0 1,1 0 1,1 1 1,0 1 1,0 0 1)),\
                                                     ((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0))))");
  BOOST_CHECK_EQUAL(algorithm::volume(*s), 1);
}

BOOST_AUTO_TEST_CASE(cubeWithHoleVolume)
{
  const std::unique_ptr<Geometry> s =
      io::readWkt("SOLID ((((0 0 0,0 0 1,0 1 1,0 1 0,0 0 0)),\
                ((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),\
                ((0 0 0,1 0 0,1 0 1,0 0 1,0 0 0)),\
                ((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0)),\
                ((0 0 1,1 0 1,1 1 1,0 1 1,0 0 1)),\
                ((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0))),\
               (((.2 .2 .2,.2 .8 .2,.2 .8 .8,.2 .2 .8,.2 .2 .2)),\
                ((.2 .2 .2,.8 .2 .2,.8 .8 .2,.2 .8 .2,.2 .2 .2)),\
                ((.2 .2 .2,.2 .2 .8,.8 .2 .8,.8 .2 .2,.2 .2 .2)),\
                ((.8 .2 .2,.8 .2 .8,.8 .8 .8,.8 .8 .2,.8 .2 .2)),\
                ((.2 .2 .8,.2 .8 .8,.8 .8 .8,.8 .2 .8,.2 .2 .8)),\
                ((.2 .8 .2,.8 .8 .2,.8 .8 .8,.2 .8 .8,.2 .8 .2))))");
  const Kernel::FT c(.6);
  const Kernel::FT ref = 1 - c * c * c;
  BOOST_CHECK(algorithm::volume(s->as<Solid>(), algorithm::NoValidityCheck()) -
                  ref <
              0.001);
}

BOOST_AUTO_TEST_CASE(invertedCubeVolume)
{
  std::unique_ptr<Geometry> const s =
      io::readWkt("SOLID ((((0 0 0,0 1 0,0 1 1,0 0 1,0 0 0)),\
                                                     ((0 0 0,1 0 0,1 1 0,0 1 0,0 0 0)),\
                                                     ((0 0 0,0 0 1,1 0 1,1 0 0,0 0 0)),\
                                                     ((1 0 0,1 0 1,1 1 1,1 1 0,1 0 0)),\
                                                     ((0 0 1,0 1 1,1 1 1,1 0 1,0 0 1)),\
                                                     ((0 1 0,1 1 0,1 1 1,0 1 1,0 1 0))))");
  BOOST_CHECK_EQUAL(algorithm::volume(*s), -1);
}

BOOST_AUTO_TEST_CASE(polyhedronVolume)
{
  const std::string block0(
      "POLYHEDRALSURFACE Z (((0 0 0, 0 1 0, 1 0 0, 0 0 0 )), "
      "((0 0 0, 1 0 0, 0 0 1, 0 0 0 )), "
      "((0 0 0, 0 0 1, 0 1 0, 0 0 0 )), "
      "((1 0 0, 0 1 0, 0 0 1, 1 0 0 )) )");

  std::unique_ptr<Geometry> geometry0(io::readWkt(block0));
  Solid const               solid(geometry0->as<PolyhedralSurface>());
  auto                      vol{algorithm::volume(solid)};
  BOOST_CHECK_EQUAL(vol * 6, 1.0);
  CGAL::Nef_polyhedron_3<SFCGAL::Kernel> const n;
}

BOOST_AUTO_TEST_SUITE_END()
