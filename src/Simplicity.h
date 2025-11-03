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
  /**
   * @brief Create a simple (valid) geometry state
   * @return Simplicity object representing simple geometry
   */
  static const Simplicity
  simple()
  {
    return Simplicity();
  }
  /**
   * @brief Create a complex (invalid) geometry state with reason
   * @param reason The reason why the geometry is complex
   * @return Simplicity object representing complex geometry
   */
  static const Simplicity
  complex(const std::string &reason)
  {
    return Simplicity(reason);
  }
  /**
   * @brief Convert to boolean indicating if geometry is simple
   * @return true if simple, false if complex
   */
  operator bool() const { return _simple; }
  /**
   * @brief Get the reason why geometry is complex
   * @return The reason string (empty if simple)
   */
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
