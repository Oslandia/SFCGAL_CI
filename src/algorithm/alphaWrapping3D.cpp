// Copyright (c) 2024-2024, SFCGAL Contributors and Oslandia
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/alphaWrapping3D.h"
#include "SFCGAL/detail/GetPointsVisitor.h"

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/alpha_wrap_3.h>

namespace SFCGAL::algorithm {

using InexactKernel   = CGAL::Exact_predicates_inexact_constructions_kernel;
using Inexact_Point_3 = InexactKernel::Point_3;
using Mesh            = CGAL::Surface_mesh<Inexact_Point_3>;
using ExactMesh       = CGAL::Surface_mesh<Kernel::Point_3>;
using EK_to_IK        = CGAL::Cartesian_converter<Kernel, InexactKernel>;
using IK_to_EK        = CGAL::Cartesian_converter<InexactKernel, Kernel>;

auto
alphaWrapping3D(const Geometry &geom, size_t relativeAlpha,
                size_t relativeOffset) -> std::unique_ptr<PolyhedralSurface>
{
  if (geom.isEmpty()) {
    return std::make_unique<PolyhedralSurface>();
  }

  // Collect points from geometry
  SFCGAL::detail::GetPointsVisitor getPointVisitor;
  const_cast<Geometry &>(geom).accept(getPointVisitor);

  // Need at least 4 points for 3D alpha wrapping
  if (getPointVisitor.points.size() < 4) {
    return std::make_unique<PolyhedralSurface>();
  }

  // Create points vector
  EK_to_IK                     toInexact;
  std::vector<Inexact_Point_3> points;
  points.reserve(getPointVisitor.points.size());
  for (const auto &point : getPointVisitor.points) {
    points.push_back(toInexact(point->toPoint_3()));
  }

  // compute alpha and offset
  CGAL::Bbox_3 bbox        = CGAL::bbox_3(points.begin(), points.end());
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));
  const double alpha       = diag_length / static_cast<double>(relativeAlpha);

  Mesh wrapMesh;
  if (relativeOffset == 0) {
    CGAL::alpha_wrap_3(points, alpha, wrapMesh);
  } else {
    const double offset = diag_length / static_cast<double>(relativeOffset);
    CGAL::alpha_wrap_3(points, alpha, offset, wrapMesh);
  }

  return std::make_unique<PolyhedralSurface>(PolyhedralSurface(wrapMesh));
}
} // namespace SFCGAL::algorithm
