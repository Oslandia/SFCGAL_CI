// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/version.hpp"

namespace SFCGAL {

// NOLINTBEGIN(modernize-avoid-c-arrays)
/// @private
static const char _sfcgal_version[] = SFCGAL_VERSION;

/// @private
static const char _sfcgal_full_version[] = SFCGAL_FULL_VERSION;
// NOLINTEND(modernize-avoid-c-arrays)

/// API C++
/// @private
auto
Version() -> const char *
{
  return _sfcgal_version;
}

/// API C++
/// @private
auto
Full_Version() -> const char *
{
  return _sfcgal_full_version;
}

} // namespace SFCGAL

// --- API C ---

extern "C" {

auto
SFCGAL_Version() -> const char *
{
  return SFCGAL::Version();
}

auto
SFCGAL_Full_Version() -> const char *
{
  return SFCGAL::Full_Version();
}
}
