/// @file text_ui.hpp
/// @brief Terminal UI components for SFCGALOP
/// @details Provides Unicode tables, colors, and formatting utilities

#ifndef SFCGALOP_TEXT_UI_HPP
#define SFCGALOP_TEXT_UI_HPP

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

/// @namespace TextUI
/// @brief Terminal UI utilities for text-based interfaces
namespace TextUI {

/// @namespace TextUI::Colors
/// @brief ANSI color escape sequences
namespace Colors {
static constexpr const char *RESET     = "\033[0m";
static constexpr const char *BOLD      = "\033[1m";
static constexpr const char *DIM       = "\033[2m";
static constexpr const char *ITALIC    = "\033[3m";
static constexpr const char *UNDERLINE = "\033[4m";

static constexpr const char *BLACK   = "\033[30m";
static constexpr const char *RED     = "\033[31m";
static constexpr const char *GREEN   = "\033[32m";
static constexpr const char *YELLOW  = "\033[33m";
static constexpr const char *BLUE    = "\033[34m";
static constexpr const char *MAGENTA = "\033[35m";
static constexpr const char *CYAN    = "\033[36m";
static constexpr const char *WHITE   = "\033[37m";

static constexpr const char *BRIGHT_BLACK   = "\033[90m";
static constexpr const char *BRIGHT_RED     = "\033[91m";
static constexpr const char *BRIGHT_GREEN   = "\033[92m";
static constexpr const char *BRIGHT_YELLOW  = "\033[93m";
static constexpr const char *BRIGHT_BLUE    = "\033[94m";
static constexpr const char *BRIGHT_MAGENTA = "\033[95m";
static constexpr const char *BRIGHT_CYAN    = "\033[96m";
static constexpr const char *BRIGHT_WHITE   = "\033[97m";

static constexpr const char *BG_BLACK   = "\033[40m";
static constexpr const char *BG_RED     = "\033[41m";
static constexpr const char *BG_GREEN   = "\033[42m";
static constexpr const char *BG_YELLOW  = "\033[43m";
static constexpr const char *BG_BLUE    = "\033[44m";
static constexpr const char *BG_MAGENTA = "\033[45m";
static constexpr const char *BG_CYAN    = "\033[46m";
static constexpr const char *BG_WHITE   = "\033[47m";
} // namespace Colors

/// @namespace TextUI::Box
/// @brief Unicode box-drawing characters
namespace Box {
// Rounded corners
static constexpr const char *ROUNDED_TOP_LEFT     = "╭";
static constexpr const char *ROUNDED_TOP_RIGHT    = "╮";
static constexpr const char *ROUNDED_BOTTOM_LEFT  = "╰";
static constexpr const char *ROUNDED_BOTTOM_RIGHT = "╯";

// Single lines
static constexpr const char *TOP_LEFT     = "┌";
static constexpr const char *TOP_RIGHT    = "┐";
static constexpr const char *BOTTOM_LEFT  = "└";
static constexpr const char *BOTTOM_RIGHT = "┘";
static constexpr const char *HORIZONTAL   = "─";
static constexpr const char *VERTICAL     = "│";
static constexpr const char *CROSS        = "┼";
static constexpr const char *T_DOWN       = "┬";
static constexpr const char *T_UP         = "┴";
static constexpr const char *T_RIGHT      = "├";
static constexpr const char *T_LEFT       = "┤";

// Double lines
static constexpr const char *DOUBLE_TOP_LEFT     = "╔";
static constexpr const char *DOUBLE_TOP_RIGHT    = "╗";
static constexpr const char *DOUBLE_BOTTOM_LEFT  = "╚";
static constexpr const char *DOUBLE_BOTTOM_RIGHT = "╝";
static constexpr const char *DOUBLE_HORIZONTAL   = "═";
static constexpr const char *DOUBLE_VERTICAL     = "║";
static constexpr const char *DOUBLE_CROSS        = "╬";
static constexpr const char *DOUBLE_T_DOWN       = "╦";
static constexpr const char *DOUBLE_T_UP         = "╩";
static constexpr const char *DOUBLE_T_RIGHT      = "╠";
static constexpr const char *DOUBLE_T_LEFT       = "╣";

// Heavy lines
static constexpr const char *HEAVY_HORIZONTAL = "━";
static constexpr const char *HEAVY_VERTICAL   = "┃";

// Symbols
static constexpr const char *BULLET      = "•";
static constexpr const char *ARROW_RIGHT = "→";
static constexpr const char *ARROW_LEFT  = "←";
static constexpr const char *ARROW_UP    = "↑";
static constexpr const char *ARROW_DOWN  = "↓";
static constexpr const char *CHECK_MARK  = "✓";
static constexpr const char *CROSS_MARK  = "✗";
static constexpr const char *STAR        = "★";
static constexpr const char *DOT         = "·";
static constexpr const char *TRIANGLE    = "▶";
} // namespace Box

// Keep old names for compatibility
constexpr const char *BOX_TOP_LEFT =
    Box::TOP_LEFT; ///< Top-left box corner character
constexpr const char *BOX_TOP_RIGHT =
    Box::TOP_RIGHT; ///< Top-right box corner character
constexpr const char *BOX_BOTTOM_LEFT =
    Box::BOTTOM_LEFT; ///< Bottom-left box corner character
constexpr const char *BOX_BOTTOM_RIGHT =
    Box::BOTTOM_RIGHT; ///< Bottom-right box corner character
constexpr const char *BOX_HORIZONTAL =
    Box::HORIZONTAL; ///< Horizontal box line character
constexpr const char *BOX_VERTICAL =
    Box::VERTICAL; ///< Vertical box line character
constexpr const char *BOX_CROSS =
    Box::CROSS; ///< Box cross intersection character
constexpr const char *BOX_T_DOWN =
    Box::T_DOWN; ///< T-junction pointing down character
constexpr const char *BOX_T_UP =
    Box::T_UP; ///< T-junction pointing up character
constexpr const char *BOX_T_RIGHT =
    Box::T_RIGHT; ///< T-junction pointing right character
constexpr const char *BOX_T_LEFT =
    Box::T_LEFT; ///< T-junction pointing left character

/**
 * @struct Column
 * @brief Table column specification
 */
struct Column {
  std::string header; ///< Column header text
  size_t      width;  ///< Column width in characters
  /**
   * @enum Alignment
   * @brief Text alignment options for table columns
   */
  enum Alignment : std::uint8_t {
    LEFT,         ///< Left-align text within column
    CENTER,       ///< Center text within column
    RIGHT         ///< Right-align text within column
  } align = LEFT; ///< Text alignment
};

/// @enum TableStyle
/// @brief Table border drawing styles
enum class TableStyle : std::uint8_t {
  SIMPLE,  ///< Single-line Unicode borders
  DOUBLE,  ///< Double-line Unicode borders
  ROUNDED, ///< Rounded corner Unicode borders
  HEAVY,   ///< Heavy-line Unicode borders
  MINIMAL, ///< Minimal borders (no verticals)
  ASCII    ///< ASCII-only characters for compatibility
};

/// @class Table
/// @brief Unicode table renderer with color support
/// @details Renders formatted tables with various border styles,
///          automatic column width calculation, and ANSI color support
class Table {
private:
  std::vector<Column>                   columns;
  std::vector<std::vector<std::string>> rows;
  bool                                  use_colors  = false;
  bool                                  use_unicode = true;
  TableStyle                            style       = TableStyle::SIMPLE;
  bool                                  show_header_separator = true;
  bool                                  show_row_separators   = false;
  std::string                           title;

public:
  /// @brief Default-constructs an empty Table.
  ///
  /// Constructs a Table with no columns or rows and default rendering options.
  /// By default the table uses Unicode borders, the SIMPLE TableStyle, shows a
  /// header separator, and does not show row separators. The title is empty.
  Table() noexcept;

  /// @brief Destructor
  ~Table() noexcept;

  /**
   * @brief Add column to table
   * @param header Column header text
   * @param width Fixed width (0 for auto)
   * @param align Text alignment
   */
  void
  add_column(const std::string &header, size_t width = 0,
             Column::Alignment align = Column::LEFT);

  /**
   * @brief Add data row to table
   * @param row Vector of cell values
   */
  void
  add_row(const std::vector<std::string> &row);
  /**
   * @brief Enable or disable ANSI color sequences when rendering the table
   * When enabled, the table will include ANSI color/style escape sequences (if
   * supported by the terminal) in headers, cells, and borders where colors are
   * used. When disabled, output will be plain text without ANSI escapes.
   * @param use true to enable colored output; false to force plain
   * (non-colored) output
   */
  void
  set_use_colors(bool use)
  {
    use_colors = use;
  }
  /**
   * @brief Enable or disable Unicode box-drawing characters for table borders
   * When set to true the table renderer will use Unicode box-drawing
   * characters (rounded, double, heavy, etc., depending on style). When set to
   * false the renderer falls back to ASCII-compatible characters suitable for
   * terminals that do not support Unicode.
   * @param use True to use Unicode borders; false to use ASCII fallback
   */
  void
  set_use_unicode(bool use)
  {
    use_unicode = use;
  }
  /**
   * @brief Set the table border drawing style
   * Changes the style used for borders when the table is rendered (affects
   * subsequent calls to print()).
   * @param table_style The TableStyle to apply
   */
  void
  set_style(TableStyle table_style)
  {
    style = table_style;
  }
  /**
   * @brief Set the table title
   * When non-empty, the title will be rendered above the table by print().
   * Passing an empty string clears any previously set title.
   * @param table_title Title text to display
   */
  void
  set_title(const std::string &table_title)
  {
    title = table_title;
  }
  /**
   * @brief Enable or disable drawing a separator line between the header and
   * the table body
   * @param show When true, a horizontal separator will be drawn after the
   * header; when false, no header separator is drawn
   */
  void
  set_show_header_separator(bool show)
  {
    show_header_separator = show;
  }
  /**
   * @brief Enable or disable drawing horizontal separators between table rows
   * When enabled, a separator line is rendered between each data row (in
   * addition to any header separator). Passing `true` turns row separators on;
   * `false` turns them off.
   * @param show Whether to show row separators
   */
  void
  set_show_row_separators(bool show)
  {
    show_row_separators = show;
  }

  void
  print() const;

private:
  void
  print_separator(const char *left, const char *middle, const char *right,
                  const char *fill) const;
  void
  print_row(const std::vector<std::string> &row, bool is_header = false) const;
  void
  print_title() const;
  [[nodiscard]] static auto
  align_text(const std::string &text, size_t width, Column::Alignment align)
      -> std::string;
  void
  get_box_chars(const char *&top_left, const char *&top_right,
                const char *&bottom_left, const char *&bottom_right,
                const char *&horizontal, const char *&vertical,
                const char *&cross, const char *&tee_down, const char *&tee_up,
                const char *&top_left_s, const char *&top_right_s) const;
};

/**
 * @brief Print formatted operations table
 * @param operations Vector of (name, category, description, params) tuples
 */
void
print_formatted_operations(
    const std::vector<std::tuple<std::string, std::string, std::string,
                                 std::string>> &operations);

/**
 * @brief Print decorated header box
 * @param title Header text
 * @param use_colors Enable ANSI colors
 */
void
print_header(const std::string &title, bool use_colors = true);

/**
 * @brief Print section title with formatting
 * @param title Section title
 * @param use_colors Enable ANSI colors
 */
void
print_section(const std::string &title, bool use_colors = true);

/**
 * @brief Print success message with green check mark
 * @param message Success message
 */
void
print_success(const std::string &message);

/**
 * @brief Print error message with red cross mark
 * @param message Error message
 */
void
print_error(const std::string &message);

/**
 * @brief Print warning message with yellow warning sign
 * @param message Warning message
 */
void
print_warning(const std::string &message);

/**
 * @brief Print info message with blue info icon
 * @param message Info message
 */
void
print_info(const std::string &message);

} // namespace TextUI

#endif // SFCGALOP_TEXT_UI_HPP
