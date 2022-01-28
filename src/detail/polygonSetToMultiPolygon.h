// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef _SFCGAL_DETAIL_POLYGONSETTOMULTIPOLYGON_H_
#define _SFCGAL_DETAIL_POLYGONSETTOMULTIPOLYGON_H_

#include <SFCGAL/config.h>

#include <CGAL/Polygon_set_2.h>
#include <SFCGAL/MultiPolygon.h>

namespace SFCGAL {
namespace detail {

/**
 * @brief convert a CGAL::Polygon_set_2 to a MultiPolygon
 * @todo unittest
 */
SFCGAL_API std::unique_ptr<MultiPolygon>
polygonSetToMultiPolygon(const CGAL::Polygon_set_2<Kernel> &polygonSet);

} // namespace detail
} // namespace SFCGAL

#endif
