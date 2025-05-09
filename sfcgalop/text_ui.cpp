// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "text_ui.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <sstream>

namespace TextUI {

// Table constructor and destructor definitions
Table::Table() noexcept  = default;
Table::~Table() noexcept = default;

// Compute the display width of a UTF-8 string, ignoring ANSI escape sequences.
// This function strips ANSI SGR sequences then counts visible characters.
static auto
visual_width(const std::string &str) -> size_t
{
  // Remove ANSI escape sequences
  static const std::regex ansi_regex("\033\\[[0-9;]*m");
  std::string             clean = std::regex_replace(str, ansi_regex, "");

  // Count visual width (UTF-8 aware)
  size_t width = 0;
  size_t index = 0;
  while (index < clean.length()) {
    auto character = static_cast<unsigned char>(clean[index]);
    if (character < 0x80) {
      // ASCII
      width++;
      index++;
    } else if ((character & 0xE0) == 0xC0) {
      // 2-byte UTF-8
      width++;
      index += 2;
    } else if ((character & 0xF0) == 0xE0) {
      // 3-byte UTF-8 (most Unicode box drawing)
      width++;
      index += 3;
    } else if ((character & 0xF8) == 0xF0) {
      // 4-byte UTF-8
      width++;
      index += 4;
    } else {
      index++;
    }
  }
  return width;
}

void
Table::add_column(const std::string &header, size_t width,
                  Column::Alignment align)
{
  Column col{header, (width != 0U) ? width : visual_width(header), align};
  columns.push_back(col);
}

void
Table::add_row(const std::vector<std::string> &row)
{
  rows.push_back(row);
  for (size_t i = 0; i < std::min(row.size(), columns.size()); ++i) {
    columns[i].width = std::max(columns[i].width, visual_width(row[i]));
  }
}

/**
 * @brief Aligns a text string to fit within a column of the given visual width.
 *
 * Computes visual width using UTF-8- and ANSI-escape-aware measurement (via
 * visual_width) and pads the returned string with ASCII spaces so the visible
 * content occupies exactly `width` columns when shorter. If the visible text
 * width is greater than or equal to `width`, the original string is returned
 * unchanged.
 *
 * Padding behavior:
 * - RIGHT: pads on the left.
 * - CENTER: splits padding between left and right; if an odd number of spaces
 * is required, the extra space is placed on the right.
 * - LEFT: pads on the right.
 *
 * @param text The input text (may contain ANSI color codes or UTF-8
 * characters).
 * @param width The target visible column width.
 * @param align Alignment directive (LEFT, CENTER, RIGHT).
 * @return std::string A new string containing the original text with added
 * ASCII-space padding so its visible width matches `width`, or the original
 * text if no padding is needed.
 */
auto
Table::align_text(const std::string &text, size_t width,
                  Column::Alignment align) -> std::string
{
  size_t text_width = visual_width(text);

  if (text_width >= width) {
    // For now, just return the text as-is if too long
    return text;
  }

  size_t padding = width - text_width;

  switch (align) {
  case Column::RIGHT:
    return std::string(padding, ' ') + text;
  case Column::CENTER: {
    size_t left_pad  = padding / 2;
    size_t right_pad = padding - left_pad;
    return std::string(left_pad, ' ') + text + std::string(right_pad, ' ');
  }
  case Column::LEFT:
  default:
    return text + std::string(padding, ' ');
  }
}

void
Table::get_box_chars(const char *&top_left, const char *&top_right,
                     const char *&bottom_left, const char *&bottom_right,
                     const char *&horizontal, const char *&vertical,
                     const char *&cross, const char *&tee_down,
                     const char *&tee_up, const char *&top_left_s,
                     const char *&top_right_s) const
{
  if (!use_unicode || style == TableStyle::ASCII) {
    top_left     = "+";
    top_right    = "+";
    bottom_left  = "+";
    bottom_right = "+";
    horizontal   = "-";
    vertical     = "|";
    cross        = "+";
    tee_down     = "+";
    tee_up       = "+";
    top_left_s   = "+";
    top_right_s  = "+";
    return;
  }

  switch (style) {
  case TableStyle::DOUBLE:
    top_left     = Box::DOUBLE_TOP_LEFT;
    top_right    = Box::DOUBLE_TOP_RIGHT;
    bottom_left  = Box::DOUBLE_BOTTOM_LEFT;
    bottom_right = Box::DOUBLE_BOTTOM_RIGHT;
    horizontal   = Box::DOUBLE_HORIZONTAL;
    vertical     = Box::DOUBLE_VERTICAL;
    cross        = Box::DOUBLE_CROSS;
    tee_down     = Box::DOUBLE_T_DOWN;
    tee_up       = Box::DOUBLE_T_UP;
    top_left_s   = Box::DOUBLE_T_RIGHT;
    top_right_s  = Box::DOUBLE_T_LEFT;
    break;
  case TableStyle::ROUNDED:
    top_left     = Box::ROUNDED_TOP_LEFT;
    top_right    = Box::ROUNDED_TOP_RIGHT;
    bottom_left  = Box::ROUNDED_BOTTOM_LEFT;
    bottom_right = Box::ROUNDED_BOTTOM_RIGHT;
    horizontal   = Box::HORIZONTAL;
    vertical     = Box::VERTICAL;
    cross        = Box::CROSS;
    tee_down     = Box::T_DOWN;
    tee_up       = Box::T_UP;
    top_left_s   = Box::T_RIGHT;
    top_right_s  = Box::T_LEFT;
    break;
  case TableStyle::HEAVY:
    top_left     = Box::TOP_LEFT;
    top_right    = Box::TOP_RIGHT;
    bottom_left  = Box::BOTTOM_LEFT;
    bottom_right = Box::BOTTOM_RIGHT;
    horizontal   = Box::HEAVY_HORIZONTAL;
    vertical     = Box::HEAVY_VERTICAL;
    cross        = Box::CROSS;
    tee_down     = Box::T_DOWN;
    tee_up       = Box::T_UP;
    top_left_s   = Box::T_RIGHT;
    top_right_s  = Box::T_LEFT;
    break;
  case TableStyle::MINIMAL:
    top_left     = " ";
    top_right    = " ";
    bottom_left  = " ";
    bottom_right = " ";
    horizontal   = " ";
    vertical     = " ";
    cross        = " ";
    tee_down     = " ";
    tee_up       = " ";
    top_left_s   = " ";
    top_right_s  = " ";
    break;
  case TableStyle::SIMPLE:
  default:
    top_left     = Box::TOP_LEFT;
    top_right    = Box::TOP_RIGHT;
    bottom_left  = Box::BOTTOM_LEFT;
    bottom_right = Box::BOTTOM_RIGHT;
    horizontal   = Box::HORIZONTAL;
    vertical     = Box::VERTICAL;
    cross        = Box::CROSS;
    tee_down     = Box::T_DOWN;
    tee_up       = Box::T_UP;
    top_left_s   = Box::T_RIGHT;
    top_right_s  = Box::T_LEFT;
    break;
  }
}

/**
 * @brief Render a horizontal separator line for the table to stdout.
 *
 * Constructs and prints a single separator row composed of the provided
 * left, middle, right and fill strings. For each column, the fill string
 * is repeated across (column.width + 2) character cells (accounting for
 * one-space padding on each side of cell content). Middle separators are
 * inserted between columns. A trailing newline is printed.
 *
 * @param left  String printed at the start of the separator (e.g., corner).
 * @param middle String printed between columns (e.g., column junction).
 * @param right String printed at the end of the separator (e.g., corner).
 * @param fill  String used to fill the horizontal span of each column.
 *
 * @note This function writes directly to std::cout.
 */
void
Table::print_separator(const char *left, const char *middle, const char *right,
                       const char *fill) const
{
  std::cout << left;
  for (size_t i = 0; i < columns.size(); ++i) {
    for (size_t j = 0; j < columns[i].width + 2; ++j) {
      std::cout << fill;
    }
    if (i < columns.size() - 1) {
      std::cout << middle;
    }
  }
  std::cout << right << "\n";
}

/**
 * @brief Render the table title row with a surrounding top border and a
 * connecting bottom border.
 *
 * If the table has no title this function returns immediately. When a title is
 * present it computes the total visual width of the table (including per-column
 * widths, padding, and inter-column separators), prints a top border, a
 * centered title line (optionally colorized when the table's color mode is
 * enabled), and a bottom border that visually connects the title area to the
 * table columns.
 *
 * Side effects:
 * - Writes the title and surrounding border lines to std::cout.
 *
 * Notes:
 * - Uses the table's box-drawing characters (via get_box_chars) and text
 * alignment logic (align_text) to produce the output.
 */
void
Table::print_title() const
{
  if (title.empty()) {
    return;
  }

  const char *top_left     = nullptr;
  const char *top_right    = nullptr;
  const char *bottom_left  = nullptr;
  const char *bottom_right = nullptr;
  const char *horizontal   = nullptr;
  const char *vertical     = nullptr;
  const char *cross        = nullptr;
  const char *tee_down     = nullptr;
  const char *tee_up       = nullptr;
  const char *top_left_s   = nullptr;
  const char *top_right_s  = nullptr;
  get_box_chars(top_left, top_right, bottom_left, bottom_right, horizontal,
                vertical, cross, tee_down, tee_up, top_left_s, top_right_s);

  // Calculate the total table width based on columns
  size_t total_width = 0;
  for (const auto &col : columns) {
    total_width += col.width + 2; // +2 for padding
  }
  if (!columns.empty()) {
    total_width += columns.size() - 1; // Add separators between columns
  }

  // Top border
  std::cout << top_left;
  for (size_t i = 0; i < total_width; ++i) {
    std::cout << horizontal;
  }
  std::cout << top_right << "\n";

  // Title line
  std::cout << vertical;
  std::string formatted_title = title;
  if (use_colors) {
    formatted_title =
        std::string(Colors::BOLD) + Colors::BRIGHT_CYAN + title + Colors::RESET;
  }
  std::cout << align_text(formatted_title, total_width, Column::CENTER);
  std::cout << vertical << "\n";

  // Bottom border of title that connects to columns
  std::cout << top_left_s;
  for (size_t i = 0; i < columns.size(); ++i) {
    for (size_t j = 0; j < columns[i].width + 2; ++j) {
      std::cout << horizontal;
    }
    if (i < columns.size() - 1) {
      std::cout << tee_down;
    }
  }
  std::cout << top_right_s << "\n";
}

/**
 * @brief Render a single table row (header or data) to stdout.
 *
 * Renders the provided cells into the table's configured columns, applying
 * column alignment and width, surrounding cell content with a single space on
 * each side, and separating columns with the table's vertical separator (or a
 * minimal spacing when using MINIMAL style). When rendering a header row and
 * colors are enabled, cell text is emitted in bold.
 *
 * @param row Vector of cell strings for this row; if shorter than the number
 *            of columns the remaining columns are rendered as empty space.
 * @param is_header If true, the row is treated as a header (affects styling
 *                  such as bold when colors are enabled).
 *
 * Side effects: writes formatted output to stdout.
 */
void
Table::print_row(const std::vector<std::string> &row, bool is_header) const
{
  const char *top_left     = nullptr;
  const char *top_right    = nullptr;
  const char *bottom_left  = nullptr;
  const char *bottom_right = nullptr;
  const char *horizontal   = nullptr;
  const char *vertical     = nullptr;
  const char *cross        = nullptr;
  const char *tee_down     = nullptr;
  const char *tee_up       = nullptr;
  const char *top_left_s   = nullptr;
  const char *top_right_s  = nullptr;
  get_box_chars(top_left, top_right, bottom_left, bottom_right, horizontal,
                vertical, cross, tee_down, tee_up, top_left_s, top_right_s);

  const char *sep = (style == TableStyle::MINIMAL) ? "  " : vertical;

  std::cout << sep;
  for (size_t i = 0; i < columns.size(); ++i) {
    std::cout << " ";
    if (i < row.size()) {
      // Apply formatting to the text before alignment
      std::string text = row[i];
      if (is_header && use_colors && !text.empty()) {
        text.insert(0, Colors::BOLD);
        text.append(Colors::RESET);
      }

      std::cout << align_text(text, columns[i].width, columns[i].align);
    } else {
      std::cout << std::string(columns[i].width, ' ');
    }
    std::cout << " " << sep;
  }
  std::cout << "\n";
}

/**
 * @brief Render the table to standard output.
 *
 * Prints the complete table (title, header, rows and separators) using the
 * configured style, column widths, padding and color settings. Behavior:
 * - Returns immediately if the table has no columns.
 * - Selects appropriate box-drawing or ASCII characters based on style and
 *   Unicode/ASCII configuration.
 * - Renders the title (if any), an optional top separator, the header row,
 *   an optional header separator, all data rows (inserting row separators if
 *   enabled), and a final bottom separator when the style is not MINIMAL.
 *
 * Side effects:
 * - Outputs formatted text to standard output (including ANSI color codes when
 *   colors are enabled).
 */
void
Table::print() const
{
  if (columns.empty()) {
    return;
  }

  const char *top_left     = nullptr;
  const char *top_right    = nullptr;
  const char *bottom_left  = nullptr;
  const char *bottom_right = nullptr;
  const char *horizontal   = nullptr;
  const char *vertical     = nullptr;
  const char *cross        = nullptr;
  const char *tee_down     = nullptr;
  const char *tee_up       = nullptr;
  const char *top_left_s   = nullptr;
  const char *top_right_s  = nullptr;
  get_box_chars(top_left, top_right, bottom_left, bottom_right, horizontal,
                vertical, cross, tee_down, tee_up, top_left_s, top_right_s);

  print_title();

  if (style != TableStyle::MINIMAL) {
    if (title.empty()) {
      print_separator(top_left, tee_down, top_right, horizontal);
    }
  }

  std::vector<std::string> headers;
  headers.reserve(columns.size());
  for (const auto &col : columns) {
    headers.emplace_back(col.header);
  }
  print_row(headers, true);

  if (show_header_separator && style != TableStyle::MINIMAL) {
    print_separator(top_left_s, cross, top_right_s, horizontal);
  }

  for (size_t i = 0; i < rows.size(); ++i) {
    print_row(rows[i], false);
    if (show_row_separators && i < rows.size() - 1 &&
        style != TableStyle::MINIMAL) {
      print_separator(top_left_s, cross, top_right_s, horizontal);
    }
  }

  if (style != TableStyle::MINIMAL) {
    print_separator(bottom_left, tee_up, bottom_right, horizontal);
  }
}

void
print_formatted_operations(
    const std::vector<std::tuple<std::string, std::string, std::string,
                                 std::string>> &operations)
{
  std::map<std::string,
           std::vector<std::tuple<std::string, std::string, std::string>>>
      by_category;

  for (const auto &[name, category, desc, params] : operations) {
    by_category[category].emplace_back(name, desc, params);
  }

  Table table;
  table.set_style(TableStyle::ROUNDED);
  table.set_use_colors(true);
  table.set_title("SFCGAL Operations");

  table.add_column("Operation", 26, Column::LEFT);
  table.add_column("Description", 48, Column::LEFT);
  table.add_column("Parameters", 21, Column::LEFT);

  for (const auto &[category, ops] : by_category) {
    std::string cat_display = std::string(Colors::BOLD) + Colors::BRIGHT_CYAN +
                              Box::TRIANGLE + " " + category + Colors::RESET;
    table.add_row({cat_display, "", ""});

    for (const auto &[name, desc, params] : ops) {
      std::string op_display =
          std::string("  ") + Colors::BRIGHT_WHITE + name + Colors::RESET;
      std::string param_display =
          params.empty() ? ""
                         : std::string(Colors::YELLOW) + params + Colors::RESET;
      table.add_row({op_display, desc, param_display});
    }
  }

  table.print();
}

void
print_header(const std::string &title, bool use_colors)
{
  const char *horizontal_char = Box::DOUBLE_HORIZONTAL;
  size_t      width           = 60;

  // Build the formatted title string
  std::string formatted_title;
  if (use_colors) {
    formatted_title = std::string(Colors::BOLD) + Colors::BRIGHT_WHITE + title +
                      Colors::RESET;
  } else {
    formatted_title = title;
  }

  // Calculate visual width (without ANSI codes)
  size_t title_visual_width = visual_width(formatted_title);
  size_t left_padding       = (width - title_visual_width) / 2;
  size_t right_padding      = width - left_padding - title_visual_width;

  if (use_colors) {
    std::cout << Colors::BRIGHT_CYAN;
  }

  std::cout << Box::DOUBLE_TOP_LEFT;
  for (size_t i = 0; i < width + 2; ++i) {
    std::cout << horizontal_char;
  }
  std::cout << Box::DOUBLE_TOP_RIGHT << "\n";

  std::cout << Box::DOUBLE_VERTICAL << " ";
  for (size_t i = 0; i < left_padding; ++i) {
    std::cout << " ";
  }
  std::cout << formatted_title;
  for (size_t i = 0; i < right_padding; ++i) {
    std::cout << " ";
  }
  std::cout << " " << Box::DOUBLE_VERTICAL << "\n";

  std::cout << Box::DOUBLE_BOTTOM_LEFT;
  for (size_t i = 0; i < width + 2; ++i) {
    std::cout << horizontal_char;
  }
  std::cout << Box::DOUBLE_BOTTOM_RIGHT;

  if (use_colors) {
    std::cout << Colors::RESET;
  }
  std::cout << "\n\n";
}

void
print_section(const std::string &title, bool use_colors)
{
  if (use_colors) {
    std::cout << Colors::BOLD << Colors::BRIGHT_BLUE << Box::TRIANGLE << " "
              << title << Colors::RESET << "\n";
  } else {
    std::cout << "=== " << title << " ===" << "\n";
  }
}

void
print_success(const std::string &message)
{
  std::cout << Colors::BRIGHT_GREEN << Box::CHECK_MARK << " " << Colors::GREEN
            << message << Colors::RESET << "\n";
}

void
print_error(const std::string &message)
{
  std::cout << Colors::BRIGHT_RED << Box::CROSS_MARK << " " << Colors::RED
            << message << Colors::RESET << "\n";
}

void
print_warning(const std::string &message)
{
  std::cout << Colors::BRIGHT_YELLOW << "⚠ " << Colors::YELLOW << message
            << Colors::RESET << "\n";
}

void
print_info(const std::string &message)
{
  std::cout << Colors::BRIGHT_BLUE << "ℹ " << Colors::BLUE << message
            << Colors::RESET << "\n";
}

} // namespace TextUI
