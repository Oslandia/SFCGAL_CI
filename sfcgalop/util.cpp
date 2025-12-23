// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

/**
 * util.cpp - Utility function implementations
 */

#include "util.hpp"
#include <algorithm>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

// Error handlers removed - using C++ exception handling

} // anonymous namespace

/**
 * @brief No-op placeholder for SFCGAL error handler setup.
 *
 * Historically this registered C-style error handlers for SFCGAL. Error
 * handling is now managed by the C++ API, so this function intentionally does
 * nothing and is retained only for API compatibility.
 */
void
setup_error_handlers()
{
  // Error handlers now managed by C++ API
}

/**
 * @brief Parse a single parameter of the form "name=value" or an unlabeled
 * "value".
 *
 * Parses param_str and sets output references:
 * - If param_str contains '=', splits on the first '=' into name (left) and
 * value (right).
 * - If no '=', sets name to an empty string and value to param_str.
 *
 * Leading and trailing spaces and tabs are trimmed from both name and value.
 *
 * @param name[out] Receives the parameter name; empty when the input had no
 * '='.
 * @param value[out] Receives the parameter value (trimmed). Must be non-empty
 * for the function to succeed.
 * @return true if a non-empty value was produced and assigned; false if
 * param_str is empty or value is empty after trimming.
 */
namespace {
// Safe trim helper that handles all-whitespace strings
void
safe_trim(std::string &str)
{
  if (str.empty()) {
    return;
  }

  size_t first = str.find_first_not_of(" \t");
  if (first == std::string::npos) {
    str.clear();
    return;
  }

  size_t last = str.find_last_not_of(" \t");
  str         = str.substr(first, last - first + 1);
}
} // namespace

auto
parse_parameter(const std::string &param_str, std::string &name,
                std::string &value) -> bool
{
  if (param_str.empty()) {
    return false;
  }

  size_t eq_pos = param_str.find('=');
  if (eq_pos != std::string::npos) {
    // Format: name=value
    name  = param_str.substr(0, eq_pos);
    value = param_str.substr(eq_pos + 1);

    // Trim whitespace safely
    safe_trim(name);
    safe_trim(value);
  } else {
    // Format: just value
    name.clear();
    value = param_str;

    // Trim whitespace safely
    safe_trim(value);
  }

  return !value.empty();
}

auto
parse_parameters(const std::string &params_str)
    -> std::vector<std::pair<std::string, std::string>>
{
  std::vector<std::pair<std::string, std::string>> result;

  if (params_str.empty()) {
    return result;
  }

  std::stringstream string_stream(params_str);
  std::string       param;

  while (std::getline(string_stream, param, ',')) {
    std::string parameter_name;
    std::string parameter_value;
    if (parse_parameter(param, parameter_name, parameter_value)) {
      result.emplace_back(parameter_name, parameter_value);
    }
  }

  return result;
}

namespace SFCGAL::sfcgalop::util {

auto
to_underscore_convention(std::string_view op_name) -> std::string
{
  std::string result;
  for (size_t i = 0; i < op_name.length(); ++i) {
    char current_char = op_name[i];

    // Add underscore before uppercase letters (except at the beginning)
    if (current_char >= 'A' && current_char <= 'Z' && i > 0) {
      // Check if previous character is lowercase or if previous is digit and
      // current is not 'D' to avoid adding underscore in sequences like "3D"
      if ((i > 0 && op_name[i - 1] >= 'a' && op_name[i - 1] <= 'z') ||
          (i > 0 && op_name[i - 1] >= '0' && op_name[i - 1] <= '9' &&
           current_char != 'D')) {
        result += '_';
      }
    }

    // Convert uppercase to lowercase
    result += static_cast<char>(std::tolower(current_char));
  }

  // Handle common patterns like 3D, 2D
  std::regex pattern(R"((\d)([A-Z]))");
  result = std::regex_replace(result, pattern, "$1_$2");

  return result;
}

/**
 * @brief Normalize operation name by converting to lowercase and removing
 * underscores This allows matching of different naming conventions (e.g.,
 * alpha_wrapping_3d, alphawrapping3d, alphaWrapping3D)
 */
auto
normalize_operation_name(std::string_view op_name) -> std::string
{
  std::string result;
  result.reserve(op_name.length());

  for (char c : op_name) {
    if (c != '_') { // Skip underscores
      result += std::tolower(c);
    }
  }

  return result;
}

} // namespace SFCGAL::sfcgalop::util
