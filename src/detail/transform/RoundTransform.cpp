// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/detail/transform/RoundTransform.h>

#include <SFCGAL/Point.h>

namespace SFCGAL {
namespace transform {

///
///
///
RoundTransform::RoundTransform(const long &scale) : _scale(scale) {}

///
///
///
void
RoundTransform::transform(Point &p)
{
  p.coordinate().round(_scale);
}

} // namespace transform
} // namespace SFCGAL
