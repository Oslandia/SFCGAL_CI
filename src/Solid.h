// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_SOLID_H_
#define SFCGAL_SOLID_H_

#include <boost/assert.hpp>
#include <vector>

#include <boost/ptr_container/ptr_vector.hpp>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

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
class SFCGAL_API Solid : public Geometry {
public:
  typedef boost::ptr_vector<PolyhedralSurface>::iterator       iterator;
  typedef boost::ptr_vector<PolyhedralSurface>::const_iterator const_iterator;

  /**
   * Empty Solid constructor
   */
  Solid();
  /**
   * Constructor with an exterior shell
   */
  Solid(const PolyhedralSurface &exteriorShell);
  /**
   * Constructor with an exterior shell (takes ownership)
   */
  Solid(PolyhedralSurface *exteriorShell);
  /**
   * Constructor with a vector of shells (PolyhedralSurface)
   */
  Solid(const std::vector<PolyhedralSurface> &shells);
  /**
   * Copy constructor
   */
  Solid(const Solid &other);
  /**
   * assign operator
   */
  Solid &
  operator=(Solid other);
  /**
   * destructor
   */
  ~Solid();

  //-- SFCGAL::Geometry
  Solid *
  clone() const override;

  //-- SFCGAL::Geometry
  std::string
  geometryType() const override;
  //-- SFCGAL::Geometry
  GeometryType
  geometryTypeId() const override;
  //-- SFCGAL::Geometry
  int
  dimension() const override;
  //-- SFCGAL::Geometry
  int
  coordinateDimension() const override;
  //-- SFCGAL::Geometry
  bool
  isEmpty() const override;
  //-- SFCGAL::Geometry
  bool
  is3D() const override;
  //-- SFCGAL::Geometry
  bool
  isMeasured() const override;

  /**
   * Returns the exterior shell
   */
  inline const PolyhedralSurface &
  exteriorShell() const
  {
    return _shells[0];
  }
  /**
   * Returns the exterior shell
   */
  inline PolyhedralSurface &
  exteriorShell()
  {
    return _shells[0];
  }

  /**
   * Returns the number of interior shells
   */
  inline size_t
  numInteriorShells() const
  {
    return _shells.size() - 1;
  }
  /**
   * Returns the n-th interior shell
   */
  inline const PolyhedralSurface &
  interiorShellN(size_t const &n) const
  {
    return _shells[n + 1];
  }
  /**
   * Returns the n-th interior shell
   */
  inline PolyhedralSurface &
  interiorShellN(size_t const &n)
  {
    return _shells[n + 1];
  }
  /**
   * adds an interior shell to the Solid
   */
  inline void
  addInteriorShell(const PolyhedralSurface &shell)
  {
    _shells.push_back(shell.clone());
  }
  /**
   * adds an interior shell to the Solid
   */
  inline void
  addInteriorShell(PolyhedralSurface *shell)
  {
    BOOST_ASSERT(shell != NULL);
    _shells.push_back(shell);
  }

  /**
   * Sets the Solid exterior shell
   */
  inline void
  setExteriorShell(const PolyhedralSurface &shell)
  {
    _shells.replace(0, shell.clone());
  }

  /**
   * Sets the Solid exterior shell
   * The ownership of the shell is taken. The caller is not responsible
   * anymore of its deallocation.
   */
  inline void
  setExteriorShell(PolyhedralSurface *shell)
  {
    _shells.replace(0, shell);
  }

  /**
   * Returns the number of shells
   */
  inline size_t
  numShells() const
  {
    return _shells.size();
  }
  /**
   * Returns the n-th shell, 0 is exteriorShell
   * @warning not standard, avoid conditionnal to access rings
   */
  inline const PolyhedralSurface &
  shellN(const size_t &n) const
  {
    BOOST_ASSERT(n < numShells());
    return _shells[n];
  }
  /**
   * Returns the n-th shell, 0 is exteriorShell
   * @warning not standard, avoid conditionnal to access rings
   */
  inline PolyhedralSurface &
  shellN(const size_t &n)
  {
    BOOST_ASSERT(n < numShells());
    return _shells[n];
  }

  //-- iterators

  inline iterator
  begin()
  {
    return _shells.begin();
  }
  inline const_iterator
  begin() const
  {
    return _shells.begin();
  }

  inline iterator
  end()
  {
    return _shells.end();
  }
  inline const_iterator
  end() const
  {
    return _shells.end();
  }

  //-- visitors

  //-- SFCGAL::Geometry
  void
  accept(GeometryVisitor &visitor) override;
  //-- SFCGAL::Geometry
  void
  accept(ConstGeometryVisitor &visitor) const override;

  /**
   * Serializer
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int /*version*/)
  {
    ar &boost::serialization::base_object<Geometry>(*this);
    ar &_shells;
  }

private:
  boost::ptr_vector<PolyhedralSurface> _shells;

  void
  swap(Solid &other)
  {
    _shells.swap(other._shells);
  }
};

} // namespace SFCGAL

#endif
