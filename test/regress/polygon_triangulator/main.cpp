// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <filesystem>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include "SFCGAL/Exception.h"

#include "SFCGAL/GeometryCollection.h"
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

#include "SFCGAL/algorithm/area.h"
#include "SFCGAL/algorithm/force2D.h"
#include "SFCGAL/algorithm/orientation.h"
#include "SFCGAL/io/vtk.h"
#include "SFCGAL/io/wkt.h"
#include "SFCGAL/triangulate/triangulatePolygon.h"

#include <boost/chrono.hpp>

using namespace SFCGAL;

namespace po = boost::program_options;

/*
 * Triangulate each polygon in an input file containing lines in the following
 * format : <id> "|" ( <wkt polygon> | <wkt multipolygon> )
 */
auto
main(int argc, char *argv[]) -> int
{
  /*
   * declare options
   */
  po::options_description desc("polygon triangulator options : ");
  desc.add_options()("help", "produce help message")(
      "progress", "display progress")("verbose", "verbose mode")(
      "force2d", "force 2d polygon")("line", po::value<int>(), "line to test")(
      "filename", po::value<std::string>(),
      "input filename (id|wkt_[multi]polygon on each line)");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help") != 0U) {
    std::cout << desc << '\n';
    return 0;
  }

  bool const verbose  = vm.count("verbose") != 0;
  bool const progress = vm.count("progress") != 0;
  bool const force2d  = vm.count("force2d") != 0;

  std::string filename;

  if (vm.count("filename") != 0U) {
    filename = vm["filename"].as<std::string>();
  } else {
    std::cerr << "missing input file" << '\n';
    std::cout << desc << '\n';
    return 1;
  }

  int oneLine = -1;

  if (vm.count("line") != 0U) {
    oneLine = vm["line"].as<int>();
  }

  /*
   * open file
   */
  std::ifstream ifs(filename.c_str());

  if (!ifs.good()) {
    std::cerr << "fail to open : " << filename << '\n';
    return 1;
  }

  const std::string tri_filename(filename + ".tri.wkt");

  std::ofstream tri_ofs(tri_filename.c_str());

  if (!tri_ofs.good()) {
    std::cerr << "fail to write : " << tri_filename << '\n';
    return 1;
  }

  const std::string error_filename(filename + ".error.wkt");

  std::ofstream ofs_error(error_filename.c_str());

  if (!ofs_error.good()) {
    std::cerr << "fail to write : " << error_filename << '\n';
    return 1;
  }

  // boost::timer timer ;
  boost::chrono::system_clock::time_point const start =
      boost::chrono::system_clock::now();

  std::vector<std::string> invalidGeom;
  std::vector<std::string> inapropriateGeom;
  /*
   * process file
   */
  int         lineNumber = 0;
  int         numFailed  = 0;
  int         numSuccess = 0;
  std::string line;

  while (std::getline(ifs, line)) {
    lineNumber++;

    if (-1 != oneLine && oneLine != lineNumber) {
      continue;
    }

    boost::algorithm::trim(line);

    if (line.empty()) {
      continue;
    }

    if (verbose) {
      std::cout << "#" << line << '\n';
      std::cout.flush();
    }

    if (progress && lineNumber % 1000 == 0) {
      std::cout.width(12);
      boost::chrono::duration<double> const elapsed =
          boost::chrono::system_clock::now() - start;
      std::cout << std::left << lineNumber << "(" << elapsed << " s)" << '\n';
    }

    std::vector<std::string> tokens;
    boost::algorithm::split(tokens, line, boost::is_any_of("|"));

    std::string const &wkt = tokens.back();
    std::string        id;

    if (tokens.size() > 1) {
      id = tokens.front();
    } else {
      std::ostringstream oss;
      oss << lineNumber;
      id = oss.str();
    }

    // std::cout << "process " << id << std::endl;

    bool                failed = true;
    TriangulatedSurface triangulatedSurface;

    try {
      std::unique_ptr<Geometry> g(io::readWkt(wkt));

      // io::vtk( *g, (boost::format("/tmp/polygon_%s.vtk") % id ).str() );

      if (force2d) {
        algorithm::force2D(*g);
      }

      triangulate::triangulatePolygon3D(*g, triangulatedSurface);

      // check area
      double const areaPolygons  = algorithm::area3D(*g);
      double const areaTriangles = algorithm::area3D(triangulatedSurface);

      double const ratio = fabs(areaPolygons - areaTriangles) /
                           std::max(areaPolygons, areaTriangles);

      if (ratio > 0.1) {
        std::cerr << filename << ":" << lineNumber << " error:" << id << "|"
                  << "area(polygon) != area(tin) ( " << areaPolygons
                  << " !=" << areaTriangles << ")"
                  << "|" << g->asText() << "|" << triangulatedSurface.asText()
                  << '\n';
      }

      numSuccess++;
      failed = false;
    } catch (InappropriateGeometryException &e) {
      inapropriateGeom.push_back(id);
      numFailed++;
    } catch (GeometryInvalidityException &e) {
      invalidGeom.push_back(id);
      numFailed++;
    } catch (std::exception &e) {
      BOOST_ASSERT_MSG(false,
                       (boost::format("%s:%d: unhandled std::exception: %s") %
                        filename % lineNumber % e.what())
                           .str()
                           .c_str());
    } catch (...) {
      BOOST_ASSERT_MSG(false, (boost::format("%s:%d: unknown exception") %
                               filename % lineNumber)
                                  .str()
                                  .c_str());
    }

    // output triangulated surface
    tri_ofs << id << "|" << failed << "|" << triangulatedSurface.asText(5)
            << '\n';
  } // end for each line

  ofs_error.close();
  tri_ofs.close();

  boost::chrono::duration<double> const elapsed =
      boost::chrono::system_clock::now() - start;

  for (auto &i : invalidGeom) {
    std::cout << "    " << i << " is invalid\n";
  }

  for (auto &i : inapropriateGeom) {
    std::cout << "    " << i << " is inapropriate for triangulation\n";
  }

  std::cout << filename << " complete (" << elapsed << " s)---" << '\n';
  std::cout << numFailed << " failed /" << (numFailed + numSuccess) << '\n';

  if (numFailed == 0) {
    // delete empty error file
    std::filesystem::remove(error_filename);
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
