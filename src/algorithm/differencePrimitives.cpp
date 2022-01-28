// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/Triangle.h>
// GOTCHA do not include the header, this is a trick to avoid ambiguous def in
// CGAL

namespace SFCGAL {
namespace algorithm {
auto
intersection(const CGAL::Triangle_3<Kernel> &a,
             const CGAL::Triangle_3<Kernel> &b) -> CGAL::Object
{
  return CGAL::intersection(a, b);
}
} // namespace algorithm
} // namespace SFCGAL
