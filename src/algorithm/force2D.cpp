// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/algorithm/force2D.h>
#include <SFCGAL/detail/transform/Force2D.h>

namespace SFCGAL {
namespace algorithm {

///
///
///
void
force2D(Geometry &g)
{
  transform::Force2D t;
  g.accept(t);
}

} // namespace algorithm
} // namespace SFCGAL
