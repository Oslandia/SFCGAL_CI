// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/Envelope.h>
#include <SFCGAL/Geometry.h>
#include <SFCGAL/io/GeometryStreams.h>

namespace SFCGAL {

///
///
///
std::ostream &
operator<<(std::ostream &ostr, const Envelope &env)
{
  return env.print(ostr);
}

///
///
///
std::ostream &
operator<<(std::ostream &ostr, const Geometry &g)
{
  ostr << g.asText();
  return ostr;
}

} // namespace SFCGAL
