// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, SFCGAL team.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_NAMESPACE_H
#define SFCGAL_NAMESPACE_H

// first we configure Doxygen groups (names and hierarchy)
/**
 * @defgroup public_api Public C++ API
 * @{
 * @defgroup algorithm Main algorithms
 *
 * @defgroup detail Implementation details
 * @{
 * @defgroup algorithm_detail Algorithm implementation details
 * @defgroup io_detail IO implementation details
 * @}
 *
 * @defgroup generator Geometry generator
 *
 * @defgroup graph Geometry graph
 *
 * @defgroup incomplete Incomplete or buggy functions
 *
 * @defgroup io Input/output tools
 *
 * @defgroup transform Geometry transform
 *
 * @}
 *
 * @defgroup capi Public C API
 */

// second, for each namespace we assign a Doxygen group and a bit of
// documentation. By doing this we do not need to set the group on each C/C++
// items
/**
 * @brief Default SFCGAL namespace
 * @ingroup public_api
 */
namespace SFCGAL {

/**
 * @brief Main algorithm namespace
 * @ingroup algorithm
 */
namespace algorithm {
}

/**
 * @brief Implementation details namespace
 * @ingroup detail
 */
namespace detail {

/**
 * @brief Algorithm detail namespace
 * @ingroup algorithm_detail
 */
namespace algorithm {
}

/**
 * @brief Input/output detail namespace
 * @ingroup io_detail
 */
namespace io {
}
} // namespace detail

/**
 * @brief graph namespace
 * @ingroup graph
 * @todo should be moved from detail to SFCGAL folder
 */
namespace graph {
/**
 * @brief Graph algorithm namespace
 * @ingroup algorithm
 */
namespace algorithm {
}
} // namespace graph

/**
 * @brief Generator namespace
 * @ingroup generator
 * @todo should be moved from detail to SFCGAL folder
 */
namespace generator {
}

/**
 * @brief Main Input/output namespace
 * @ingroup io
 */
namespace io {

/**
 * @brief OBJ file format tools
 */
namespace OBJ {
}

/**
 * @brief VKT file format tools
 */
namespace VKT {
}
} // namespace io

/**
 * @brief transform namespace
 * @ingroup transform
 * @todo should be moved from detail to SFCGAL folder
 */
namespace transform {
}

/**
 * @brief Triangulate dedicated namespace
 * @ingroup algorithm
 * @todo should be moved into algorithm namespace/folder
 */
namespace triangulate {
/**
 * @brief Triangulate detail namespace
 * @ingroup detail
 */
namespace detail {
}
} // namespace triangulate

} // namespace SFCGAL

#endif
