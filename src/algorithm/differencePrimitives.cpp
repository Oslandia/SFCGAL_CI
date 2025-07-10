// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/Triangle.h"
// GOTCHA do not include the header, this is a trick to avoid ambiguous def in
// CGAL

namespace SFCGAL::algorithm {
auto
intersection(const CGAL::Triangle_3<Kernel> &a,
             const CGAL::Triangle_3<Kernel> &b) -> CGAL::Object
{
  return CGAL::intersection(a, b);
}
} // namespace SFCGAL::algorithm
