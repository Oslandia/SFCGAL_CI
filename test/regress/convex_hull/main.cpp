// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

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

#include "SFCGAL/algorithm/convexHull.h"
#include "SFCGAL/io/wkt.h"

#include <boost/chrono.hpp>

#include <boost/filesystem.hpp>

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
  po::options_description desc("convex hull options : ");
  desc.add_options()("help", "produce help message")(
      "progress", "display progress")("verbose", "verbose mode")(
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

  std::string filename;

  if (vm.count("filename") != 0U) {
    filename = vm["filename"].as<std::string>();
  } else {
    std::cerr << "missing input file" << '\n';
    std::cout << desc << '\n';
    return 1;
  }

  /*
   * open file
   */
  std::ifstream ifs(filename.c_str());

  if (!ifs.good()) {
    std::cerr << "fail to open : " << filename << '\n';
    return 1;
  }

  std::string const tri_filename(filename + ".convex.wkt");
  std::ofstream     ofs_result(tri_filename.c_str());

  if (!ofs_result.good()) {
    std::cerr << "fail to write : " << tri_filename << '\n';
    return 1;
  }

  std::string const error_filename(filename + ".error.wkt");
  std::ofstream     ofs_error(error_filename.c_str());

  if (!ofs_error.good()) {
    std::cerr << "fail to write : " << error_filename << '\n';
    return 1;
  }

  // boost::timer timer ;
  boost::chrono::system_clock::time_point const start =
      boost::chrono::system_clock::now();

  /*
   * process file
   */
  int         lineNumber = 0;
  int         numFailed  = 0;
  int         numSuccess = 0;
  std::string line;

  while (std::getline(ifs, line)) {
    lineNumber++;

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

    bool failed = true;

    std::unique_ptr<Geometry> hull;
    std::unique_ptr<Geometry> hull3D;

    try {
      std::unique_ptr<Geometry> g;
      g      = io::readWkt(wkt);
      hull   = algorithm::convexHull(*g);
      hull3D = algorithm::convexHull3D(*g);
      failed = false;
    } catch (Exception &e) {
      std::cerr << "[Exception]" << id << "|" << e.what() << "|" << wkt << '\n';
    } catch (std::exception &e) {
      std::cerr << "[std::exception]" << id << "|" << e.what() << "|" << wkt
                << '\n';
    } catch (...) {
      std::cerr << "[...]" << id << "|" << wkt << '\n';
    }

    if (failed) {
      numFailed++;
      ofs_error << line << '\n';
    } else {
      numSuccess++;
    }

    // output triangulated surface
    ofs_result << id << "|" << failed << "|" << hull->asText(5) << '\n';
  } // end for each line

  ofs_error.close();
  ofs_result.close();

  boost::chrono::duration<double> const elapsed =
      boost::chrono::system_clock::now() - start;
  std::cout << filename << " complete (" << elapsed << " s)---" << '\n';
  std::cout << numFailed << " failed /" << (numFailed + numSuccess) << '\n';

  if (numFailed == 0) {
    // delete empty error file
    boost::filesystem::remove(error_filename);
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
