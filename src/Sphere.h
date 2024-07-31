// Copyright (c) 2024-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef _SFCGAL_SPHERE_H_
#define _SFCGAL_SPHERE_H_

#include <optional>
#include <vector>

#include <CGAL/Polyhedron_3.h>

#include "SFCGAL/Kernel.h"
#include "SFCGAL/export.h"

namespace SFCGAL {

/**
 * @brief Represents a sphere in 3D space
 *
 * This class provides methods to generate a polyhedron and a point cloud
 * representation of a sphere. It uses SFCGAL's Kernel for exact computations.
 */
class SFCGAL_API Sphere {
public:
  /**
   * @brief Constructs a Sphere object
   * @param radius The radius of the sphere
   * @param center The center point of the sphere
   * @param num_vertical The number of vertical divisions
   * @param num_horizontal The number of horizontal divisions
   */
  Sphere(const Kernel::FT      &radius = 1.0,
         const Kernel::Point_3 &center = Kernel::Point_3(0, 0, 0),
         int num_vertical = 16, int num_horizontal = 32);

  /**
   * @brief Copy constructor
   */
  Sphere(const Sphere &other) = default;

  /**
   * @brief Assignment operator
   */
  Sphere &
  operator=(Sphere other);

  /**
   * @brief Destructor
   */
  ~Sphere() = default;

  /**
   * @brief Sets the radius of the sphere
   * @param radius The new radius
   */
  inline void
  setRadius(const Kernel::FT &radius)
  {
    m_radius = radius;
    invalidateCache();
  }

  /**
   * @brief Sets the center of the sphere
   * @param center The new center point
   */
  inline void
  setCenter(const Kernel::Point_3 &center)
  {
    m_center = center;
    invalidateCache();
  }

  /**
   * @brief Sets the number of vertical divisions
   * @param num The new number of vertical divisions
   */
  inline void
  setNumVertical(int num)
  {
    m_num_vertical = num;
    invalidateCache();
  }

  /**
   * @brief Sets the number of horizontal divisions
   * @param num The new number of horizontal divisions
   */
  inline void
  setNumHorizontal(int num)
  {
    m_num_horizontal = num;
    invalidateCache();
  }

  /**
   * @brief Gets the radius of the sphere
   * @return The radius
   */
  inline const Kernel::FT &
  radius() const
  {
    return m_radius;
  }

  /**
   * @brief Gets the center of the sphere
   * @return The center point
   */
  inline const Kernel::Point_3 &
  center() const
  {
    return m_center;
  }

  /**
   * @brief Gets the number of vertical divisions
   * @return The number of vertical divisions
   */
  inline int
  numVertical() const
  {
    return m_num_vertical;
  }

  /**
   * @brief Gets the number of horizontal divisions
   * @return The number of horizontal divisions
   */
  inline int
  numHorizontal() const
  {
    return m_num_horizontal;
  }

  /**
   * @brief Generates a polyhedron representation of the sphere
   * @return A CGAL::Polyhedron_3 object representing the sphere
   */
  CGAL::Polyhedron_3<Kernel>
  generatePolyhedron();

  /**
   * @brief Generates a point cloud representation of the sphere
   * @return A vector of Point_3 objects representing points on the sphere's
   * surface
   */
  std::vector<Kernel::Point_3>
  generatePoints();

private:
  Kernel::FT                                  m_radius;
  Kernel::Point_3                             m_center;
  int                                         m_num_vertical;
  int                                         m_num_horizontal;
  std::optional<CGAL::Polyhedron_3<Kernel>>   m_polyhedron;
  std::optional<std::vector<Kernel::Point_3>> m_points;

  void
  invalidateCache();
  CGAL::Polyhedron_3<Kernel>
  generateSpherePolyhedron();
  std::vector<Kernel::Point_3>
  generateSpherePoints();
};

} // namespace SFCGAL

#endif // _SFCGAL_SPHERE_H_
