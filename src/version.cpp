// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/version.h>

namespace SFCGAL {

const char _sfcgal_version[]      = SFCGAL_VERSION;
const char _sfcgal_full_version[] = SFCGAL_FULL_VERSION;

auto
Version() -> const char *
{
  return _sfcgal_version;
}

auto
Full_Version() -> const char *
{
  return _sfcgal_full_version;
}

} // namespace SFCGAL
