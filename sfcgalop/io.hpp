#ifndef SFCGALOP_IO_HPP
#define SFCGALOP_IO_HPP

#include <SFCGAL/Geometry.h>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <variant>

enum class OutputFormat : std::uint8_t { WKT = 0, WKB = 1, TXT = 2 };

auto
load_geometry(const std::string &source) -> std::unique_ptr<SFCGAL::Geometry>;

using OperationResult =
    std::variant<std::unique_ptr<SFCGAL::Geometry>, bool, double, std::string>;

namespace IO {
void
print_result(const std::optional<OperationResult> &result, OutputFormat format,
             int precision);
}

auto
parse_output_format(const char *format_str, OutputFormat &format) -> bool;

#endif // SFCGALOP_IO_HPP
