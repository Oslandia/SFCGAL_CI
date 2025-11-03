// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_VALIDITY_H_
#define SFCGAL_VALIDITY_H_

namespace SFCGAL {

/**
 * @brief the class, convertible to bool, that stores the reason why a geom is
 * invalid
 */
struct Validity {
  /**
   * @note the class has private ctor to force the use of functions valid() and
   * invalid(reason) that are clearer in the code than to remember that "Valid
   * constructed with a reason is invalid"
   */
  /**
   * @brief Create a valid geometry status
   * @return Valid geometry status
   */
  static auto
  valid() -> const Validity
  {
    return Validity();
  }
  /**
   * @brief Create an invalid geometry status with reason
   * @param reason The reason why the geometry is invalid
   * @return Invalid geometry status
   */
  static auto
  invalid(const std::string &reason) -> const Validity
  {
    return Validity(reason);
  }
  /**
   * @brief Convert to bool (check if geometry is valid)
   * @return true if geometry is valid, false otherwise
   */
  operator bool() const { return _valid; }
  /**
   * @brief Get the reason why geometry is invalid
   * @return The reason string (empty if valid)
   */
  [[nodiscard]] auto
  reason() const -> const std::string &
  {
    return _reason;
  }

private:
  bool        _valid; // not const to allow default copy
  std::string _reason;
  /**
   * @brief default ctor for valid
   */
  Validity() : _valid(true) {}
  /**
   * @brief if we construct with a reason, the class is invalid
   */
  Validity(const std::string &reason) : _valid(false), _reason(reason) {}
};

} // namespace SFCGAL
#endif
