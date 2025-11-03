// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_ISVALID_H_
#define SFCGAL_ALGORITHM_ISVALID_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Validity.h"
#include "SFCGAL/algorithm/force2D.h"
#include "SFCGAL/algorithm/force3D.h"

namespace SFCGAL {

/**
 * Functions used to assert for geometry validity
 * @note exception message is apparently limited in length, thus print the
 * reason for invalidity before its text representation (that can be very long)
 */

/**
 * @brief Assert geometry validity without context
 * @param geometry The geometry to validate
 * @return void
 */
void SFCGAL_API
SFCGAL_ASSERT_GEOMETRY_VALIDITY(const Geometry &geometry);

/**
 * @brief Assert 2D geometry validity
 * @param geometry The geometry to validate in 2D
 * @return void
 */
void SFCGAL_API
SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(const Geometry &geometry);

/**
 * @brief Assert 3D geometry validity
 * @param geometry The geometry to validate in 3D
 * @return void
 */
void SFCGAL_API
SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(const Geometry &geometry);

/**
 * @brief Assert geometry validity on a plane (not implemented)
 * @param geometry The geometry to validate
 * @return void
 */
void SFCGAL_API
SFCGAL_ASSERT_GEOMETRY_VALIDITY_ON_PLANE(const Geometry &geometry);

namespace algorithm {

/**
 * @brief Check validity of generic Geometry
 * @param geometry The Geometry to validate
 * @param toleranceAbs Absolute tolerance for validation
 * @return Validity status
 */
SFCGAL_API auto
isValid(const Geometry &geometry, const double &toleranceAbs = 1e-9)
    -> Validity;

/**
 * @brief Propagate validity flag to geometry
 * @param geometry The geometry to mark
 * @param valid The validity flag to set
 */
SFCGAL_API void
propagateValidityFlag(Geometry &geometry, bool valid);

/**
 * Tag used for variants of algorithm that do not do validity check
 */
struct NoValidityCheck {};

} // namespace algorithm
} // namespace SFCGAL

#endif
