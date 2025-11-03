// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/transform/RoundTransform.h"

#include "SFCGAL/Point.h"

namespace SFCGAL::transform {

RoundTransform::RoundTransform(const long &scale) : _scale(scale) {}

void
RoundTransform::transform(Point &point)
{
  point.coordinate().round(_scale);
}

} // namespace SFCGAL::transform
