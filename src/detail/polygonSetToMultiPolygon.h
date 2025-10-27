// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_DETAIL_POLYGONSETTOMULTIPOLYGON_H_
#define SFCGAL_DETAIL_POLYGONSETTOMULTIPOLYGON_H_

#include "SFCGAL/config.h"

#include "SFCGAL/MultiPolygon.h"
#include <CGAL/Polygon_set_2.h>

namespace SFCGAL::detail {

/**
 * @brief convert a CGAL::Polygon_set_2 to a MultiPolygon
 * @param polygonSet The CGAL polygon set to convert
 * @return A unique pointer to the resulting MultiPolygon
 * @todo unittest
 */
SFCGAL_API auto
polygonSetToMultiPolygon(const CGAL::Polygon_set_2<Kernel> &polygonSet)
    -> std::unique_ptr<MultiPolygon>;

} // namespace SFCGAL::detail

#endif
