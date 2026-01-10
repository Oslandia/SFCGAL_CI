// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGALOP_IO_HPP
#define SFCGALOP_IO_HPP

#include <SFCGAL/Geometry.h>
#include <cstdint>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <variant>

enum class OutputFormat : std::uint8_t { WKT = 0, WKB = 1, TXT = 2, OBJ = 3, GEOJSON = 4 };

[[nodiscard]] auto
load_geometry(const std::string &source) -> std::unique_ptr<SFCGAL::Geometry>;

using OperationResult =
    std::variant<std::unique_ptr<SFCGAL::Geometry>, bool, double, std::string>;

namespace IO {
/**
 * @brief Writes an OperationResult to an output stream using the requested
 * format and precision.
 *
 * This prints the value contained in `result`:
 * - If the variant holds `std::unique_ptr<SFCGAL::Geometry>` and the pointer is
 * non-null:
 *   - `OutputFormat::WKT` or `OutputFormat::TXT`: prints WKT text using
 * `precision`.
 *   - `OutputFormat::WKB`: prints the geometry's WKB as lowercase hex bytes.
 *   - `OutputFormat::OBJ`: prints the geometry in OBJ format.
 *   - A null geometry produces no output.
 * - If the variant holds `bool`: prints "true" or "false".
 * - If the variant holds `double`: prints the numeric value with the given
 * `precision`.
 * - If the variant holds `std::string`: prints the string as-is.
 *
 * @param result Optional variant containing the operation result (geometry,
 * bool, double, or string).
 * @param format Output format to use for geometry values (WKT, WKB, TXT, or
 * OBJ).
 * @param precision Number of digits of precision for textual numeric output.
 * @param out Output stream to write to (defaults to std::cout).
 */
void
print_result(const std::optional<OperationResult> &result, OutputFormat format,
             int precision, std::ostream &out = std::cout);
} // namespace IO

[[nodiscard]] auto
parse_output_format(const char *format_str, OutputFormat &format) -> bool;

#endif // SFCGALOP_IO_HPP
