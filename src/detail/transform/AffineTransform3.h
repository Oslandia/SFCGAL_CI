// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TRANSFORM_AFFINETRANSFORM3_H_
#define SFCGAL_TRANSFORM_AFFINETRANSFORM3_H_

#include "SFCGAL/Kernel.h"
#include "SFCGAL/Transform.h"
#include "SFCGAL/config.h"

#include <CGAL/Aff_transformation_3.h>

namespace SFCGAL::transform {

/**
 * Wrapper for CGAL::Aff_transform_3
 * @todo unittest
 */
class SFCGAL_API AffineTransform3 : public Transform {
public:
  /**
   * @brief Constructor with a transform
   * @param transform The CGAL 3D affine transformation
   */
  AffineTransform3(CGAL::Aff_transformation_3<Kernel> transform);

  /**
   * @brief Transform a point using the 3D affine transformation
   * @param point The point to transform
   */
  void
  transform(Point &point) override;

  /**
   * @brief Transform a linestring using the 3D affine transformation
   * @param linestring The linestring to transform
   */
  virtual void
  transform(LineString &linestring);
  /**
   * @brief Transform a triangle using the 3D affine transformation
   * @param triangle The triangle to transform
   */
  virtual void
  transform(Triangle &triangle);
  /**
   * @brief Transform a polygon using the 3D affine transformation
   * @param polygon The polygon to transform
   */
  virtual void
  transform(Polygon &polygon);

  /**
   * @brief Transform a polyhedral surface using the 3D affine transformation
   * @param surface The polyhedral surface to transform
   */
  virtual void
  transform(PolyhedralSurface &surface);

  /**
   * @brief Transform a triangulated surface using the 3D affine transformation
   * @param surface The triangulated surface to transform
   */
  virtual void
  transform(TriangulatedSurface &surface);

  /**
   * @brief Transform a solid using the 3D affine transformation
   * @param solid The solid to transform
   */
  virtual void
  transform(Solid &solid);

private:
  CGAL::Aff_transformation_3<Kernel> _transform;
};

} // namespace SFCGAL::transform

#endif
