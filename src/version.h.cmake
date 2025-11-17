// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_VERSION_H_
#define SFCGAL_VERSION_H_

#include "SFCGAL/export.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SFCGAL_VERSION_MAJOR @SFCGAL_VERSION_MAJOR@
#define SFCGAL_VERSION_MINOR @SFCGAL_VERSION_MINOR@
#define SFCGAL_VERSION_PATCH @SFCGAL_VERSION_PATCH@

#define SFCGAL_VERSION "@SFCGAL_VERSION@"

// CGAL version macros
#define SFCGAL_CGAL_VERSION_MAJOR @SFCGAL_CGAL_VERSION_MAJOR@
#define SFCGAL_CGAL_VERSION_MINOR @SFCGAL_CGAL_VERSION_MINOR@
#define SFCGAL_CGAL_VERSION_PATCH @SFCGAL_CGAL_VERSION_PATCH@

#define SFCGAL_FULL_VERSION \
    "SFCGAL @SFCGAL_VERSION@, CGAL @CGAL_VERSION@, BOOST @Boost_VERSION_STRING@"

// --- API C ---

SFCGAL_API const char *SFCGAL_Version();
SFCGAL_API const char *SFCGAL_Full_Version();

#ifdef __cplusplus
}
#endif

#endif
