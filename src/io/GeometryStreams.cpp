// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/io/GeometryStreams.h"
#include "SFCGAL/Envelope.h"
#include "SFCGAL/Geometry.h"

namespace SFCGAL {

///
///
///
auto
operator<<(std::ostream &ostr, const Envelope &env) -> std::ostream &
{
  return env.print(ostr);
}

///
///
///
auto
operator<<(std::ostream &ostr, const Geometry &g) -> std::ostream &
{
  ostr << g.asText();
  return ostr;
}

} // namespace SFCGAL
