// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/Triangle.h>
// GOTCHA do not include the header, this is a trick to avoid ambiguous def in
// CGAL

namespace SFCGAL {
namespace algorithm {
CGAL::Object
intersection(const CGAL::Triangle_3<Kernel> &a,
             const CGAL::Triangle_3<Kernel> &b)
{
  return CGAL::intersection(a, b);
}
} // namespace algorithm
} // namespace SFCGAL
