// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_SIMPLICITY_H_
#define SFCGAL_SIMPLICITY_H_

namespace SFCGAL {

/**
 * @brief the class, convertible to bool, that stores the reason why a geom is
 * not simple
 */
struct Simplicity {
  /**
   * @note the class has private ctor to force the use of functions valid() and
   * invalid(reason) that are clearer in the code than to remember that "Valid
   * constructed with a reason is invalid"
   */
  static const Simplicity
  simple()
  {
    return Simplicity();
  }
  static const Simplicity
  complex(const std::string &reason)
  {
    return Simplicity(reason);
  }
  operator bool() const { return _simple; }
  const std::string &
  reason() const
  {
    return _reason;
  }

private:
  bool        _simple; // not const to allow default copy
  std::string _reason;
  /**
   * @brief default ctor for simple
   */
  Simplicity() : _simple(true) {}
  /**
   * @brief if we construct with a reason, the class is complex
   */
  Simplicity(const std::string &reason) : _simple(false), _reason(reason) {}
};

} // namespace SFCGAL
#endif
