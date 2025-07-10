// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/version.h"

namespace SFCGAL {

/// @private
const char _sfcgal_version[] = SFCGAL_VERSION;
/// @private
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
