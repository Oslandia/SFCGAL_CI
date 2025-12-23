// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

/**
 * util.hpp - Utility function declarations
 */

#ifndef SFCGALOP_UTIL_HPP
#define SFCGALOP_UTIL_HPP

#include <string>
#include <string_view>
#include <utility>
#include <vector>

/**
 * Setup error handlers for SFCGAL
 */
void
setup_error_handlers();

/**
 * Parse a parameter string in the format "name=value" or just "value"
 *
 * @param param_str Parameter string
 * @param name Reference to store the parameter name (empty if not found)
 * @param value Reference to store the parameter value
 * @return true if parsing was successful, false otherwise
 */
auto
parse_parameter(const std::string &param_str, std::string &name,
                std::string &value) -> bool;

/**
 * Parse multiple parameters from a string
 *
 * @param params_str Parameters string (e.g., "radius=5,segments=10")
 * @return Vector of name-value pairs
 */
auto
parse_parameters(const std::string &params_str)
    -> std::vector<std::pair<std::string, std::string>>;

namespace SFCGAL::sfcgalop::util {

/**
 * @brief Structure to hold the result of an operation match
 */
struct OperationMatchResult {
  /**
   * @brief Indicates whether the operation was found
   */
  bool found;
  /**
   * @brief The resolved name after alias processing
   */
  std::string resolved_name;
  /**
   * @brief Iterator to the operation in the operations vector
   */
  void *operation_it; // Using void* as a placeholder; will be cast
                      // appropriately in implementation
};

/**
 * @brief Convert operation name to underscore convention for display purposes.
 * This function converts camelCase names to snake_case for consistent display.
 */
auto
to_underscore_convention(std::string_view op_name) -> std::string;

/**
 * @brief Normalize operation name by converting to lowercase and removing
 * underscores This allows matching of different naming conventions (e.g.,
 * alpha_wrapping_3d, alphawrapping3d, alphaWrapping3D)
 */
auto
normalize_operation_name(std::string_view op_name) -> std::string;

} // namespace SFCGAL::sfcgalop::util

#endif /* SFCGALOP_UTIL_HPP */
