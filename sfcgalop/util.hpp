// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

/**
 * util.hpp - Utility function declarations
 */

#ifndef SFCGALOP_UTIL_HPP
#define SFCGALOP_UTIL_HPP

#include <string>
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

#endif /* SFCGALOP_UTIL_HPP */
