// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "io.hpp"

#include <SFCGAL/io/OBJ.h>
#include <SFCGAL/io/ewkt.h>
#include <SFCGAL/io/geojson.h>
#include <SFCGAL/io/wkb.h>
#include <SFCGAL/io/wkt.h>

#include <algorithm>
#include <boost/endian/arithmetic.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#ifdef _WIN32
  #include <fcntl.h>
  #include <io.h>
#endif

/**
 * @brief Load a geometry from various sources with WKT, WKB, and HEXWKB
 * support.
 *
 * Reads geometry data from one of three sources and returns a parsed SFCGAL
 * geometry. If `source` is empty, returns nullptr. If `source` equals "stdin",
 * the function reads all of standard input, trims surrounding whitespace, and
 * auto-detects the format (WKT, WKB, or HEXWKB). Otherwise it treats `source`
 * as a filesystem path: if the file can be opened its entire contents (trimmed)
 * are parsed with format auto-detection; if the file cannot be opened, `source`
 * itself is parsed.
 *
 * Format detection:
 * - HEXWKB: Detects "0x" prefix, binary WKB patterns (00/01), or hex digit
 * strings
 * - WKB: Converts detected hex strings to binary and uses SFCGAL::io::readWkb
 * - WKT: Falls back to SFCGAL::io::readWkt for non-hex data
 *
 * All exceptions from file I/O and SFCGAL parsing are caught and return
 * nullptr.
 *
 * @param source Input selector: either "stdin", a filesystem path, or geometry
 * data.
 * @return std::unique_ptr<SFCGAL::Geometry> Parsed geometry, or nullptr on
 * error.
 */
// NOLINTBEGIN(readability-function-cognitive-complexity)
auto
load_geometry(const std::string &source) -> std::unique_ptr<SFCGAL::Geometry>
{
  if (source.empty()) {
    return nullptr;
  }

  // Helper to safely trim whitespace
  auto trim = [](const std::string &str) -> std::string {
    if (str.empty()) {
      return str;
    }

    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) {
      return {};
    }

    size_t last = str.find_last_not_of(" \t\n\r");
    if (last == std::string::npos || first > last) {
      return {};
    }

    return str.substr(first, (last - first + 1));
  };

  // Helper to detect and parse WKB (binary/hex), WKT, or OBJ
  // NOLINTBEGIN(readability-function-cognitive-complexity)
  auto try_parse_geometry =
      [&trim](const std::string &input) -> std::unique_ptr<SFCGAL::Geometry> {
    std::string data = trim(input);
    if (data.empty()) {
      return nullptr;
    }

    try {
      // Check if this looks like an OBJ file
      // OBJ files typically start with comments (#) or vertex declarations (v)
      std::istringstream iss(data);
      std::string        first_line;
      bool               looks_like_obj = false;

      while (std::getline(iss, first_line)) {
        first_line = trim(first_line);
        if (first_line.empty()) {
          continue; // Skip empty lines
        }
        if (first_line[0] == '#') {
          continue; // Skip comments
        }
        if (first_line.size() >= 2 && first_line[0] == 'v' &&
            first_line[1] == ' ') {
          looks_like_obj = true;
          break;
        }
        if (first_line.size() >= 2 && first_line[0] == 'f' &&
            first_line[1] == ' ') {
          looks_like_obj = true;
          break;
        }
        if (first_line.size() >= 2 && first_line[0] == 'l' &&
            first_line[1] == ' ') {
          looks_like_obj = true;
          break;
        }
        if (first_line.size() >= 2 && first_line[0] == 'p' &&
            first_line[1] == ' ') {
          looks_like_obj = true;
          break;
        }
        // If we find a line that doesn't look like OBJ, stop checking
        break;
      }

      if (looks_like_obj) {
        return SFCGAL::io::OBJ::load(data);
      }
      // Check if this looks like GeoJSON
      // GeoJSON typically starts with '{' and contains "type" field
      if (!data.empty() && data[0] == '{') {
        // Check for common GeoJSON type values
        if (data.find("\"type\"") != std::string::npos &&
            (data.find("\"Point\"") != std::string::npos ||
             data.find("\"LineString\"") != std::string::npos ||
             data.find("\"Polygon\"") != std::string::npos ||
             data.find("\"MultiPoint\"") != std::string::npos ||
             data.find("\"MultiLineString\"") != std::string::npos ||
             data.find("\"MultiPolygon\"") != std::string::npos ||
             data.find("\"GeometryCollection\"") != std::string::npos ||
             data.find("\"Feature\"") != std::string::npos ||
             data.find("\"FeatureCollection\"") != std::string::npos)) {
          return SFCGAL::io::readGeoJSON(data);
        }
      }
      // First, detect raw (binary) WKB: first byte 0x00/0x01 or presence of
      // NUL/control bytes.
      if (!data.empty() && (static_cast<unsigned char>(data[0]) == 0x00 ||
                            static_cast<unsigned char>(data[0]) == 0x01)) {
        return SFCGAL::io::readWkb(data);
      }
      if (std::any_of(data.begin(), data.end(), [](unsigned char character) {
            return character < 0x09 || (character < 0x20 && character != '\n' &&
                                        character != '\r' && character != '\t');
          })) {
        return SFCGAL::io::readWkb(data);
      }

      // Check for HEXWKB patterns
      bool is_hex = false;
      if (data.length() >= 2) {
        // Check for "0x" prefix
        if (data.substr(0, 2) == "0x" || data.substr(0, 2) == "0X") {
          data   = data.substr(2);
          is_hex = true;
        }
        // Plain hex digits (even length)
        else if (data.length() % 2 == 0 &&
                 std::all_of(data.begin(), data.end(),
                             [](unsigned char character) {
                               return std::isxdigit(character) != 0;
                             })) {
          is_hex = true;
        }
      }

      if (is_hex) {
        // Fast hex -> binary
        auto hex_to_int = [](unsigned char character) -> int {
          if (character >= '0' && character <= '9') {
            return character - '0';
          }
          character = static_cast<unsigned char>(std::tolower(character));
          if (character >= 'a' && character <= 'f') {
            return 10 + (character - 'a');
          }
          return -1;
        };
        std::string bin;
        bin.reserve(data.size() / 2);
        for (size_t i = 0; i < data.size(); i += 2) {
          int high_nibble = hex_to_int(static_cast<unsigned char>(data[i]));
          int low_nibble  = hex_to_int(static_cast<unsigned char>(data[i + 1]));
          if (high_nibble < 0 || low_nibble < 0) {
            return SFCGAL::io::readWkt(input); // fallback
          }
          bin.push_back(static_cast<char>((high_nibble << 4) | low_nibble));
        }
        return SFCGAL::io::readWkb(bin);
      }
      // Try WKT parsing
      return SFCGAL::io::readWkt(data);
    } catch (const std::exception &) {
      // Return nullptr on any parsing error
      return nullptr;
    }
  };
  // NOLINTEND(readability-function-cognitive-complexity)

  try {
    if (source == "stdin") {
      // Set stdin to binary mode on Windows to preserve WKB bytes
#ifdef _WIN32
      _setmode(_fileno(stdin), _O_BINARY);
#endif
      std::ostringstream string_stream(std::ios::binary);
      string_stream << std::cin.rdbuf();
      return try_parse_geometry(string_stream.str());
    }

    // Open in binary mode to preserve raw WKB bytes on Windows
    std::ifstream file(source, std::ios::binary);
    if (file.good()) {
      std::ostringstream string_stream(std::ios::binary);
      string_stream << file.rdbuf();
      file.close();
      return try_parse_geometry(string_stream.str());
    }
    // File doesn't exist or can't be opened, treat source as geometry data
    return try_parse_geometry(source);
  } catch (const std::exception &) {
    // Return nullptr for any I/O or parsing errors
    return nullptr;
  }
}
// NOLINTEND(readability-function-cognitive-complexity)

namespace IO {

void
print_result(const std::optional<OperationResult> &result, OutputFormat format,
             int precision, std::ostream &out)
{
  if (!result.has_value()) {
    return;
  }

  std::visit(
      [format, precision, &out](auto &&arg) {
        using T = std::decay_t<decltype(arg)>;

        if constexpr (std::is_same_v<T, std::unique_ptr<SFCGAL::Geometry>>) {
          if (arg) {
            switch (format) {
            case OutputFormat::WKT:
              out << arg->asText(precision) << "\n";
              break;
            case OutputFormat::WKB:
              out << arg->asWkb(boost::endian::order::native, true) << "\n";
              break;
            case OutputFormat::TXT:
              // Use WKT with full precision as EWKT equivalent
              out << arg->asText(precision) << "\n";
              break;
            case OutputFormat::OBJ:
              out << SFCGAL::io::OBJ::saveToString(*arg);
              break;
            case OutputFormat::GEOJSON:
              out << SFCGAL::io::writeGeoJSON(*arg) << "\n";
              break;
            }
          }
        } else if constexpr (std::is_same_v<T, bool>) {
          out << (arg ? "true" : "false") << "\n";
        } else if constexpr (std::is_same_v<T, double>) {
          out << std::setprecision(precision) << arg << "\n";
        } else if constexpr (std::is_same_v<T, std::string>) {
          out << arg << "\n";
        }
      },
      result.value());
}

} // namespace IO

/**
 * @brief Parse a case-insensitive output format name into an OutputFormat enum.
 *
 * Recognizes "wkt" -> OutputFormat::WKT, "wkb" -> OutputFormat::WKB,
 * "txt" or "ewkt" -> OutputFormat::TXT, and "obj" -> OutputFormat::OBJ.
 * Comparison is performed in lowercase.
 *
 * @param format_str Null-terminated C string containing the format name; if
 * null the function returns false.
 * @param[out] format Destination enum set on success.
 * @return true if the input matched a known format and `format` was set; false
 * otherwise.
 */
auto
parse_output_format(const char *format_str, OutputFormat &format) -> bool
{
  if (format_str == nullptr) {
    return false;
  }

  std::string fmt(format_str);
  std::transform(fmt.begin(), fmt.end(), fmt.begin(), ::tolower);

  if (fmt == "wkt") {
    format = OutputFormat::WKT;
  } else if (fmt == "wkb") {
    format = OutputFormat::WKB;
  } else if (fmt == "txt" || fmt == "ewkt") {
    format = OutputFormat::TXT;
  } else if (fmt == "obj") {
    format = OutputFormat::OBJ;
  } else if (fmt == "geojson" || fmt == "json") {
    format = OutputFormat::GEOJSON;
  } else {
    return false;
  }
  return true;
}
