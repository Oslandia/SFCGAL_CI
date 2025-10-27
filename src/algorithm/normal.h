// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_NORMAL_H_
#define SFCGAL_ALGORITHM_NORMAL_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Polygon.h"

namespace SFCGAL::algorithm {

/**
 * Returns the 3D normal to 3 consecutive points.
 * @param point1 the first point
 * @param point2 the second point
 * @param point3 the third point
 * @return the 3D normal vector
 */
template <typename Kernel>
CGAL::Vector_3<Kernel>
normal3D(const CGAL::Point_3<Kernel> &point1,
         const CGAL::Point_3<Kernel> &point2,
         const CGAL::Point_3<Kernel> &point3)
{
  // bc ^ ba
  return CGAL::cross_product(point3 - point2, point1 - point2);
}

/**
 * Returns the 3D normal to a ring (supposed to be planar and closed).
 * @param lineString the line string representing the ring
 * @param exact whether to use exact computation (default: true)
 * @return the 3D normal vector
 * @warning exact allows to avoid double rounding at the end of the computation
 */
template <typename Kernel>
CGAL::Vector_3<Kernel>
normal3D(const LineString &lineString, bool exact = true)
{
  // Newell's formula
  typename Kernel::FT nx, ny, nz;
  nx = ny = nz = 0.0;

  for (size_t i = 0; i < lineString.numPoints(); ++i) {
    const Point &pointI = lineString.pointN(i);
    const Point &pointJ = lineString.pointN((i + 1) % lineString.numPoints());
    typename Kernel::FT zi = pointI.z();
    typename Kernel::FT zj = pointJ.z();
    nx += (pointI.y() - pointJ.y()) * (zi + zj);
    ny += (zi - zj) * (pointI.x() + pointJ.x());
    nz += (pointI.x() - pointJ.x()) * (pointI.y() + pointJ.y());
  }

  if (exact) {
    return CGAL::Vector_3<Kernel>(nx, ny, nz);
  } else {
    return CGAL::Vector_3<Kernel>(CGAL::to_double(nx), CGAL::to_double(ny),
                                  CGAL::to_double(nz));
  }
}

/**
 * Returns the 3D normal to a polygon (supposed to be planar).
 * @param polygon the polygon to compute the normal for
 * @param exact whether to use exact computation (default: true)
 * @return the 3D normal vector
 * @warning exact allows to avoid double rounding at the end of the computation
 */
template <typename Kernel>
CGAL::Vector_3<Kernel>
normal3D(const Polygon &polygon, bool exact = true)
{
  return normal3D<Kernel>(polygon.exteriorRing(), exact);
}

} // namespace SFCGAL::algorithm

#endif
