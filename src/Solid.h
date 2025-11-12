// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_SOLID_H_
#define SFCGAL_SOLID_H_

#include <boost/assert.hpp>
#include <memory>
#include <vector>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/unique_ptr.hpp>
#include <boost/serialization/vector.hpp>

#include "SFCGAL/DereferenceIterator.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"

namespace SFCGAL {

/**
 * A Solid modeled with an exteriorShell and interiorShells materialized by
 * PolyhedralSurface.
 * @note A shell is supposed to be closed.
 * @warning GM_Solid, from ISO 19107 is defined in CityGML, but not in SFA.
 * Without Solid concept,
 * @note Volume concept is missing.
 */
class SFCGAL_API Solid : public GeometryImpl<Solid, Geometry> {
public:
  /// @brief Iterator type for solid shells
  using iterator = DereferenceIterator<
      std::vector<std::unique_ptr<PolyhedralSurface>>::iterator>;
  /// @brief Const iterator type for solid shells
  using const_iterator = DereferenceIterator<
      std::vector<std::unique_ptr<PolyhedralSurface>>::const_iterator>;

  /**
   * Empty Solid constructor
   */
  Solid();
  /**
   * Constructor with an exterior shell
   * @param exteriorShell The exterior shell of the solid
   */
  Solid(const PolyhedralSurface &exteriorShell);
  /**
   * Constructor with an exterior shell (takes ownership)
   * @param exteriorShell The exterior shell of the solid
   */
  Solid(PolyhedralSurface *exteriorShell);
  /**
   * Constructor with a vector of shells (PolyhedralSurface)
   * @param shells Vector of polyhedral surfaces forming the solid
   */
  Solid(const std::vector<PolyhedralSurface> &shells);
  /**
   * Copy constructor
   * @param other The solid to copy from
   */
  Solid(const Solid &other);
  /**
   * assign operator
   * @param other The solid to assign from
   * @return Reference to this solid
   */
  auto
  operator=(Solid other) -> Solid &;
  /**
   * destructor
   */
  ~Solid() override;

  //-- SFCGAL::Geometry
  /// @brief Get the geometry type as string
  /// @return "Solid"
  [[nodiscard]] auto
  geometryType() const -> std::string override;
  //-- SFCGAL::Geometry
  /// @brief Get the geometry type identifier
  /// @return TYPE_SOLID
  [[nodiscard]] auto
  geometryTypeId() const -> GeometryType override;
  //-- SFCGAL::Geometry
  /// @brief Get the dimension of the solid
  /// @return 3 (solids are 3-dimensional)
  [[nodiscard]] auto
  dimension() const -> int override;
  //-- SFCGAL::Geometry
  /// @brief Get the coordinate dimension
  /// @return Number of coordinates per point
  [[nodiscard]] auto
  coordinateDimension() const -> int override;
  //-- SFCGAL::Geometry
  /// @brief Check if the solid is empty
  /// @return true if empty, false otherwise
  [[nodiscard]] auto
  isEmpty() const -> bool override;
  //-- SFCGAL::Geometry
  /// @brief Check if the solid has 3D coordinates
  /// @return true if 3D, false otherwise
  [[nodiscard]] auto
  is3D() const -> bool override;
  //-- SFCGAL::Geometry
  /// @brief Check if the solid has measured coordinates
  /// @return true if measured, false otherwise
  [[nodiscard]] auto
  isMeasured() const -> bool override;

  /// @brief Drop Z coordinate from all surfaces
  /// @return true if Z was dropped, false otherwise
  auto
  dropZ() -> bool override;

  /// @brief Drop M coordinate from all surfaces
  /// @return true if M was dropped, false otherwise
  auto
  dropM() -> bool override;

  /// @brief Swap X and Y coordinates of all surfaces
  auto
  swapXY() -> void override;

  /**
   * Returns the exterior shell
   * @return Const reference to the exterior shell
   */
  [[nodiscard]] auto
  exteriorShell() const -> const PolyhedralSurface &
  {
    return *_shells[0];
  }
  /**
   * Returns the exterior shell
   * @return Reference to the exterior shell
   */
  auto
  exteriorShell() -> PolyhedralSurface &
  {
    return *_shells[0];
  }

  /**
   * Returns the number of interior shells
   * @return Number of interior shells
   */
  [[nodiscard]] auto
  numInteriorShells() const -> size_t
  {
    return _shells.size() - 1;
  }
  /**
   * Returns the n-th interior shell
   * @param n The index of the interior shell to get
   * @return Const reference to the nth interior shell
   */
  [[nodiscard]] auto
  interiorShellN(size_t const &n) const -> const PolyhedralSurface &
  {
    return *_shells[n + 1];
  }
  /**
   * Returns the n-th interior shell
   * @param n The index of the interior shell to get
   * @return Reference to the nth interior shell
   */
  auto
  interiorShellN(size_t const &n) -> PolyhedralSurface &
  {
    return *_shells[n + 1];
  }
  /**
   * adds an interior shell to the Solid
   * @param shell The polyhedral surface to add as interior shell
   */
  void
  addInteriorShell(const PolyhedralSurface &shell)
  {
    addInteriorShell(shell.clone());
  }
  /**
   * adds an interior shell to the Solid
   * @param shell The polyhedral surface to add as interior shell
   * @deprecated The unique_ptr version should be used instead
   */
  void
  addInteriorShell(PolyhedralSurface *shell)
  {
    addInteriorShell(std::unique_ptr<PolyhedralSurface>(shell));
  }
  /**
   * @brief Adds an interior shell to the Solid.
   *
   * @param shell A unique pointer to the PolyhedralSurface object representing
   * the new shell to add. Ownership of this shell is moved into the Solid.
   */
  void
  addInteriorShell(std::unique_ptr<PolyhedralSurface> shell)
  {
    BOOST_ASSERT(shell != nullptr);
    _shells.push_back(std::move(shell));
  }

  /**
   * @brief Sets the Solid exterior shell.
   *
   * @param shell A unique pointer to the new shell. Ownership is transferred to
   * this class.
   */
  void
  setExteriorShell(std::unique_ptr<PolyhedralSurface> shell)
  {
    _shells[0] = std::move(shell);
  }

  /**
   * @brief Sets the Solid exterior shell
   * @param shell PolyhedralSurface to set as exterior shell
   */
  void
  setExteriorShell(const PolyhedralSurface &shell)
  {
    setExteriorShell(shell.clone());
  }

  /**
   * @brief Sets the Solid exterior shell
   * @param shell Pointer to PolyhedralSurface to set as exterior shell
   * @note The ownership of the shell is taken. The caller is not responsible
   * anymore of its deallocation.
   * @deprecated The unique_ptr version should be used instead
   */
  void
  setExteriorShell(PolyhedralSurface *shell)
  {
    setExteriorShell(std::unique_ptr<PolyhedralSurface>(shell));
  }

  /**
   * @brief Returns the number of shells
   * @return Number of shells in the solid
   */
  [[nodiscard]] auto
  numShells() const -> size_t
  {
    return _shells.size();
  }
  /**
   * @brief Returns the n-th shell, 0 is exteriorShell
   * @param n Index of the shell to get
   * @return Const reference to the nth shell
   * @warning not standard, avoid conditionnal to access rings
   */
  [[nodiscard]] auto
  shellN(const size_t &n) const -> const PolyhedralSurface &
  {
    BOOST_ASSERT(n < numShells());
    return *_shells[n];
  }
  /**
   * @brief Returns the n-th shell, 0 is exteriorShell
   * @param n Index of the shell to get
   * @return Reference to the nth shell
   * @warning not standard, avoid conditionnal to access rings
   */
  auto
  shellN(const size_t &n) -> PolyhedralSurface &
  {
    BOOST_ASSERT(n < numShells());
    return *_shells[n];
  }

  //-- iterators

  /**
   * @brief Get iterator to beginning of shells
   * @return Iterator to first shell
   */
  auto
  begin() -> iterator
  {
    return dereference_iterator(_shells.begin());
  }
  /**
   * @brief Get const iterator to beginning of shells
   * @return Const iterator to first shell
   */
  [[nodiscard]] auto
  begin() const -> const_iterator
  {
    return dereference_iterator(_shells.begin());
  }

  /**
   * @brief Get iterator to end of shells
   * @return Iterator to past-the-end
   */
  auto
  end() -> iterator
  {
    return dereference_iterator(_shells.end());
  }
  /**
   * @brief Get const iterator to end of shells
   * @return Const iterator to past-the-end
   */
  [[nodiscard]] auto
  end() const -> const_iterator
  {
    return dereference_iterator(_shells.end());
  }

  //-- visitors

  //-- SFCGAL::Geometry
  /// @brief Accept a geometry visitor
  /// @param visitor Visitor to accept
  void
  accept(GeometryVisitor &visitor) override;
  //-- SFCGAL::Geometry
  /// @brief Accept a const geometry visitor
  /// @param visitor Const visitor to accept
  void
  accept(ConstGeometryVisitor &visitor) const override;

  /**
   * @brief Convert Solid to CGAL::Surface_mesh (only exterior shell)
   * @return CGAL Surface_mesh representation of the solid's exterior shell
   */
  auto
  toSurfaceMesh() const -> Surface_mesh_3;

  /**
   * @brief Serializer
   * @param ar Archive for serialization
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int /*version*/)
  {
    ar &boost::serialization::base_object<Geometry>(*this);
    ar & _shells;
  }

private:
  std::vector<std::unique_ptr<PolyhedralSurface>> _shells;

  void
  swap(Solid &other) noexcept
  {
    _shells.swap(other._shells);
  }
};

} // namespace SFCGAL

#endif
