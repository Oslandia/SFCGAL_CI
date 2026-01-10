// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "error_handler.hpp"
#include "io.hpp"
#include "operations/operations.hpp"
#include "text_ui.hpp"

#include <SFCGAL/version.h>
#include <boost/version.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <string_view>

#ifndef _WIN32
  #include <getopt.h>
#else
  #include "getopt_win.hpp"
#endif

namespace {

constexpr std::string_view PROGRAM_NAME = "sfcgalop";

/**
 * @brief Get the sfcgalop program version.
 *
 * Returns the SFCGAL library version as the program version, since
 * sfcgalop is tied to the SFCGAL library version.
 *
 * @return std::string The current SFCGAL version.
 */
auto
get_program_version() -> std::string
{
  return SFCGAL_Version();
}

struct Options {
  std::string  source_a;
  std::string  source_b;
  OutputFormat format            = OutputFormat::WKT;
  int          precision         = 6;
  bool         quiet             = false;
  bool         show_time         = false;
  bool         verbose           = false;
  bool         help_requested    = false;
  bool         version_requested = false;
  bool         list_operations   = false;
  bool         validate_mode     = false;
  std::string  operation;
  std::string  op_arg;
  std::string  help_operation;
};

/**
 * @brief Print the CLI help screen for the sfcgalop executable to stdout.
 *
 * Prints usage, available command-line options (grouped), a short overview
 * note, and several example command lines. The output includes formatted
 * tables and colored/styled text via TextUI utilities.
 *
 * Side effects:
 * - Writes the help text to standard output.
 */
void
print_help()
{
  TextUI::print_header("SFCGALOP - SFCGAL Geometry Operations CLI", true);

  std::cout << TextUI::Colors::BOLD << "Usage: " << TextUI::Colors::RESET
            << PROGRAM_NAME << " [OPTIONS...] [operation] [operation_arg]\n\n";

  std::cout << TextUI::Colors::DIM
            << "Note: If no operation is specified, the input geometry will be "
               "displayed in the requested format.\n\n"
            << TextUI::Colors::RESET;

  TextUI::print_section("Options", true);

  TextUI::Table options_table;
  options_table.set_style(TextUI::TableStyle::MINIMAL);
  options_table.set_use_colors(true);
  options_table.add_column("Option", 25, TextUI::Column::LEFT);
  options_table.add_column("Description", 55, TextUI::Column::LEFT);

  options_table.add_row(
      {"-a, --geom-a=ARG", "Source for geometry A (WKT, WKB, file, stdin)"});
  options_table.add_row(
      {"-b, --geom-b=ARG", "Source for geometry B (WKT, WKB, file, stdin)"});
  options_table.add_row(
      {"-f, --format=ARG",
       "Output format: wkt, wkb, txt/ewkt, obj, geojson/json (default: wkt)"});
  options_table.add_row(
      {"-p, --precision=N", "Decimal precision for output (default: 6)"});
  options_table.add_row({""});
  options_table.add_row({std::string(TextUI::Colors::BRIGHT_CYAN) +
                             "Output Control:" + TextUI::Colors::RESET,
                         ""});
  options_table.add_row({"-q, --quiet", "Disable result output"});
  options_table.add_row({"-t, --time", "Print execution time"});
  options_table.add_row({"-v, --verbose", "Verbose output"});
  options_table.add_row({"-l, --list", "List all available operations"});
  options_table.add_row(
      {"--validate", "Validate geometries before operations"});
  options_table.add_row({""});
  options_table.add_row({std::string(TextUI::Colors::BRIGHT_CYAN) +
                             "Information:" + TextUI::Colors::RESET,
                         ""});
  options_table.add_row({"-h, --help", "Show this help message"});
  options_table.add_row({"-V, --version", "Show version information"});
  options_table.add_row(
      {"--help-op=OPERATION", "Show help for a specific operation"});
  options_table.print();

  std::cout << "\n";
  TextUI::print_section("Quick Operations Overview", true);
  std::cout << "  Use '--list' or '-l' to see all available operations with "
               "descriptions\n\n";

  TextUI::print_section("Examples", true);
  std::cout << TextUI::Colors::DIM
            << "  # Calculate distance between two points"
            << TextUI::Colors::RESET << "\n";
  std::cout << "  " << PROGRAM_NAME
            << " -a \"POINT(0 0)\" -b \"POINT(1 1)\" distance\n\n";

  std::cout << TextUI::Colors::DIM << "  # Calculate area of a polygon"
            << TextUI::Colors::RESET << "\n";
  std::cout << "  " << PROGRAM_NAME
            << " -a \"POLYGON((0 0, 0 10, 10 10, 10 0, 0 0))\" area\n\n";

  std::cout << TextUI::Colors::DIM << "  # Compute convex hull from WKT file"
            << TextUI::Colors::RESET << "\n";
  std::cout << "  " << PROGRAM_NAME << " -a geometries.wkt convexhull\n\n";

  std::cout << TextUI::Colors::DIM << "  # Convert geometry to WKB format"
            << TextUI::Colors::RESET << "\n";
  std::cout << "  " << PROGRAM_NAME
            << " -a \"POLYGON((0 0, 0 10, 10 10, 10 0, 0 0))\" -f wkb\n\n";

  std::cout << TextUI::Colors::DIM << "  # Convert geometry to OBJ format"
            << TextUI::Colors::RESET << "\n";
  std::cout << "  " << PROGRAM_NAME
            << " -a \"TRIANGLE((0 0 0, 1 0 0, 0 1 0, 0 0 0))\" -f obj\n\n";

  std::cout << TextUI::Colors::DIM << "  # Display geometry from file"
            << TextUI::Colors::RESET << "\n";
  std::cout << "  " << PROGRAM_NAME << " -a geometry.wkt\n\n";

  std::cout << TextUI::Colors::DIM << "  # Validate geometry and show details"
            << TextUI::Colors::RESET << "\n";
  std::cout << "  " << PROGRAM_NAME << " -a complex.wkt --validate -v\n";
}

/**
 * @brief Print version information for the program and its components.
 *
 * Prints a small table of version metadata to stdout including:
 * - the program version (sfcgalop),
 * - the linked SFCGAL library version,
 * - the C++ standard used by the build (C++17),
 * and, when available, the compiler version.
 *
 * The output is formatted for human-readable terminal display and ends with a
 * trailing newline.
 */
void
print_version()
{
  TextUI::print_header("Version Information", true);

  TextUI::Table version_table;
  version_table.set_style(TextUI::TableStyle::ROUNDED);
  version_table.set_use_colors(true);
  version_table.add_column("Component", 20, TextUI::Column::LEFT);
  version_table.add_column("Version", 30, TextUI::Column::LEFT);

  version_table.add_row({"sfcgalop", get_program_version()});
  version_table.add_row({"SFCGAL Library", SFCGAL_Version()});
  std::string cgal_version = std::to_string(SFCGAL_CGAL_VERSION_MAJOR) + "." +
                             std::to_string(SFCGAL_CGAL_VERSION_MINOR) + "." +
                             std::to_string(SFCGAL_CGAL_VERSION_PATCH);
  version_table.add_row({"CGAL Library", cgal_version});

  // Format Boost version from BOOST_VERSION macro
  std::string boost_version = std::to_string(BOOST_VERSION / 100000) + "." +
                              std::to_string((BOOST_VERSION / 100) % 1000) +
                              "." + std::to_string(BOOST_VERSION % 100);
  version_table.add_row({"Boost Library", boost_version});

  version_table.add_row({"C++ Standard", "C++17"});

#ifdef __VERSION__
  version_table.add_row({"Compiler", __VERSION__});
#endif

  version_table.print();
  std::cout << "\n";
}

/**
 * @brief Print a formatted, categorized list of all available SFCGAL
 * operations.
 *
 * Retrieves operation metadata from Operations::get_all_operations_info() and
 * renders a rounded, colorized table showing each operation name, a short
 * description, and requirements (parameters or default "1 geometry").
 * Categories are printed as section headers and grouped together. Finally,
 * prints a one-line hint about using --help-op for detailed per-operation help.
 *
 * Side effects: writes formatted output to stdout.
 */
void
print_operations_list()
{
  TextUI::print_header("SFCGAL Available Operations", true);

  auto operations = Operations::get_all_operations_info();

  TextUI::Table table;
  table.set_style(TextUI::TableStyle::ROUNDED);
  table.set_use_colors(true);
  table.set_title("Complete Operations Reference");

  table.add_column("Operation", 20, TextUI::Column::LEFT);
  table.add_column("Input", 12, TextUI::Column::LEFT);
  table.add_column("Output", 8, TextUI::Column::LEFT);
  table.add_column("Description", 48, TextUI::Column::LEFT);

  std::string last_category;
  for (const auto &[name, category, desc, input, output] : operations) {
    if (category != last_category) {
      if (!last_category.empty()) {
        table.add_row({""});
      }
      std::string cat_header =
          std::string(TextUI::Colors::BOLD) + TextUI::Colors::BRIGHT_CYAN +
          TextUI::Box::TRIANGLE + " " + category + TextUI::Colors::RESET;
      table.add_row({cat_header, "", "", ""});
      last_category = category;
    }

    std::string op_name = std::string("  ") + TextUI::Colors::BRIGHT_WHITE +
                          name + TextUI::Colors::RESET;

    std::string colored_input =
        std::string(TextUI::Colors::YELLOW) + input + TextUI::Colors::RESET;
    std::string colored_output =
        std::string(TextUI::Colors::GREEN) + output + TextUI::Colors::RESET;

    table.add_row({op_name, colored_input, colored_output, desc});
  }

  table.print();

  std::cout << "\n";
  TextUI::print_info(
      "Use --help-op=<operation> for detailed help on a specific operation");
}

/**
 * @brief Parse command-line arguments into an Options structure.
 *
 * Parses argc/argv for supported long and short options (geom-a, geom-b,
 * format, precision, quiet, time, verbose, list, help, version, help-op,
 * validate) and captures a required operation name with an optional operation
 * argument.
 *
 * On success returns an Options populated from the parsed flags and positional
 * arguments. Returns std::nullopt on parsing or validation errors (invalid
 * precision, unknown output format, unexpected extra arguments, missing
 * operation when one is required). If only the program name is supplied the
 * returned Options will have help_requested set to true.
 *
 * Side effects:
 * - Prints user-facing error/help hints via TextUI::print_error and std::cerr
 *   for invalid inputs or usage issues.
 *
 * @param argc Number of command-line arguments.
 * @param argv Command-line argument vector.
 * @return std::optional<Options> Populated Options on success, or std::nullopt
 * on error.
 */
// NOLINTBEGIN(readability-function-cognitive-complexity)
auto
parse_options(int argc, char **argv) -> std::optional<Options>
{
  if ((argv == nullptr) || argc < 0) {
    return std::nullopt;
  }

  Options options;

  if (argc == 1) {
    options.help_requested = true;
    return options;
  }
  static std::array<struct option, 13> long_options = {
      {{"geom-a", required_argument, nullptr, 'a'},
       {"geom-b", required_argument, nullptr, 'b'},
       {"format", required_argument, nullptr, 'f'},
       {"precision", required_argument, nullptr, 'p'},
       {"quiet", no_argument, nullptr, 'q'},
       {"time", no_argument, nullptr, 't'},
       {"verbose", no_argument, nullptr, 'v'},
       {"list", no_argument, nullptr, 'l'},
       {"help", no_argument, nullptr, 'h'},
       {"version", no_argument, nullptr, 'V'},
       {"help-op", required_argument, nullptr, 'H'},
       {"validate", no_argument, nullptr, 'X'},
       {nullptr, 0, nullptr, 0}}};

  int option_char;
  int option_index = 0;

  while ((option_char = getopt_long(argc, argv, "a:b:f:p:qtvlhV",
                                    long_options.data(), &option_index)) !=
         -1) {
    switch (option_char) {
    case 'a':
      options.source_a = optarg;
      break;
    case 'b':
      options.source_b = optarg;
      break;
    case 'f':
      if (!parse_output_format(optarg, options.format)) {
        TextUI::print_error(std::string("Unknown format: ") + optarg);
        return std::nullopt;
      }
      break;
    case 'p':
      try {
        options.precision = std::stoi(optarg);
        if (options.precision < 0) {
          TextUI::print_error(std::string("Invalid precision: ") + optarg);
          return std::nullopt;
        }
      } catch (...) {
        TextUI::print_error(std::string("Invalid precision: ") + optarg);
        return std::nullopt;
      }
      break;
    case 'q':
      options.quiet = true;
      break;
    case 't':
      options.show_time = true;
      break;
    case 'v':
      options.verbose = true;
      break;
    case 'l':
      options.list_operations = true;
      return options;
    case 'h':
      options.help_requested = true;
      return options;
    case 'V':
      options.version_requested = true;
      return options;
    case 'H':
      options.help_operation = optarg;
      return options;
    case 'X':
      options.validate_mode = true;
      break;
    case '?':
      std::cerr << "Try '" << PROGRAM_NAME
                << " --help' for more information.\n";
      return std::nullopt;
    default:
      return std::nullopt;
    }
  }

  if (optind < argc) {
    options.operation = argv[optind++];

    if (optind < argc) {
      options.op_arg = argv[optind++];

      if (optind < argc) {
        TextUI::print_error("Unexpected extra arguments");
        return std::nullopt;
      }
    }
  } else if (!options.help_requested && !options.version_requested &&
             !options.list_operations && options.help_operation.empty() &&
             options.source_a.empty()) {
    TextUI::print_error("Missing operation or geometry");
    std::cerr << "Try '" << PROGRAM_NAME << " --help' for more information.\n";
    return std::nullopt;
  }

  // Check if operation is a constructor (doesn't need input geometry)
  bool is_constructor =
      !options.operation.empty() && (options.operation.find("make_") == 0);

  if (options.source_a.empty() && !options.operation.empty() &&
      !is_constructor) {
    options.source_a = "stdin";
  }

  return options;
}
// NOLINTEND(readability-function-cognitive-complexity)

/**
 * @brief Display a geometry in the requested format without performing any
 * operation.
 *
 * Loads geometry A, performs optional validation, and outputs it in the
 * specified format. This allows for geometry format conversion and display.
 *
 * @param options Command-line options and flags that control inputs, output
 * format, precision, validation and verbosity.
 * @return true if the geometry was loaded and displayed successfully; false on
 * any error (failed geometry load/validation when in validate mode, or caught
 * exception).
 */
auto
handle_geometry_display(const Options &options) -> bool
{
  if (options.verbose) {
    TextUI::print_section("Geometry Display", true);
    std::cout << "  " << TextUI::Box::BULLET
              << " Source: " << TextUI::Colors::BRIGHT_WHITE << options.source_a
              << TextUI::Colors::RESET << "\n";
    std::cout << "  " << TextUI::Box::BULLET
              << " Format: " << TextUI::Colors::BRIGHT_WHITE;
    switch (options.format) {
    case OutputFormat::WKT:
      std::cout << "WKT";
      break;
    case OutputFormat::WKB:
      std::cout << "WKB";
      break;
    case OutputFormat::TXT:
      std::cout << "TXT";
      break;
    case OutputFormat::OBJ:
      std::cout << "OBJ";
      break;
    case OutputFormat::GEOJSON:
      std::cout << "GEOJSON";
      break;
    }
    std::cout << TextUI::Colors::RESET << "\n\n";
  }

  auto geom_a = load_geometry(options.source_a);
  if (!geom_a) {
    TextUI::print_error("Failed to load geometry");
    return false;
  }

  if (options.validate_mode || options.verbose) {
    auto validation = ErrorHandler::validate_geometry(*geom_a);
    if (!validation.valid) {
      TextUI::print_warning("Geometry validation issues:");
      ErrorHandler::report_validation_error(validation);
      if (options.validate_mode) {
        return false;
      }
    } else if (options.verbose) {
      TextUI::print_success("Geometry is valid");
      ErrorHandler::print_geometry_info(*geom_a);
    }
  }

  if (!options.quiet) {
    // Create an OperationResult from the geometry to use existing print
    // function
    OperationResult                result = geom_a->clone();
    std::optional<OperationResult> optional_result =
        std::make_optional(std::move(result));
    IO::print_result(optional_result, options.format, options.precision);
  }

  return true;
}

/**
 * @brief Execute the requested geometry operation using the provided CLI
 * options.
 *
 * Performs loading of geometry A (and optional geometry B), optional
 * validation, invokes the requested operation, prints the result (unless
 * quiet), and may print execution time. Validation and verbose flags control
 * additional diagnostic output. All runtime exceptions thrown by the operation
 * are caught and reported.
 *
 * @param options Command-line options and flags that control inputs, output
 * format, precision, validation and verbosity.
 * @return true if the operation completed successfully and a result was
 * produced; false on any error (missing operation, failed geometry
 * load/validation when in validate mode, operation error, or caught exception).
 */
// NOLINTBEGIN(readability-function-cognitive-complexity)
auto
handle_operation(const Options &options) -> bool
{
  if (options.operation.empty()) {
    return false;
  }

  if (options.verbose) {
    TextUI::print_section("Operation Details", true);
    std::cout << "  " << TextUI::Box::BULLET
              << " Operation: " << TextUI::Colors::BRIGHT_WHITE
              << options.operation << TextUI::Colors::RESET << "\n";
    if (!options.op_arg.empty()) {
      std::cout << "  " << TextUI::Box::BULLET
                << " Arguments: " << TextUI::Colors::BRIGHT_WHITE
                << options.op_arg << TextUI::Colors::RESET << "\n";
    }
    std::cout << "\n";
  }

  std::unique_ptr<SFCGAL::Geometry> geom_a;

  // Check if operation is a constructor (doesn't need input geometry)
  bool is_constructor = options.operation.find("make_") == 0;

  if (!is_constructor) {
    geom_a = load_geometry(options.source_a);
    if (!geom_a) {
      TextUI::print_error("Failed to load geometry A");
      return false;
    }
  }

  if (geom_a && (options.validate_mode || options.verbose)) {
    auto validation = ErrorHandler::validate_geometry(*geom_a);
    if (!validation.valid) {
      TextUI::print_warning("Geometry A validation issues:");
      ErrorHandler::report_validation_error(validation);
      if (options.validate_mode) {
        return false;
      }
    } else if (options.verbose) {
      TextUI::print_success("Geometry A is valid");
      ErrorHandler::print_geometry_info(*geom_a);
    }
  }

  // Check if the operation requires a second geometry
  bool requires_geom_b =
      Operations::operation_requires_second_geometry(options.operation);

  std::unique_ptr<SFCGAL::Geometry> geom_b;
  if (requires_geom_b) {
    // Fail early if geometry B is required but not provided
    if (options.source_b.empty()) {
      TextUI::print_error("Operation '" + options.operation +
                          "' requires a second geometry (-b option)");
      return false;
    }

    geom_b = load_geometry(options.source_b);
    if (!geom_b) {
      TextUI::print_error("Failed to load geometry B");
      return false;
    }

    // Validate geometry B when required or when verbose/validate mode enabled
    if (options.validate_mode || options.verbose) {
      auto validation = ErrorHandler::validate_geometry(*geom_b);
      if (!validation.valid) {
        TextUI::print_warning("Geometry B validation issues:");
        ErrorHandler::report_validation_error(validation);
        if (options.validate_mode) {
          return false;
        }
      } else if (options.verbose) {
        TextUI::print_success("Geometry B is valid");
        ErrorHandler::print_geometry_info(*geom_b);
      }
    }
  } else if (!options.source_b.empty()) {
    // Load and optionally validate geometry B even if not required
    // (user provided it, so they might want to use it)
    geom_b = load_geometry(options.source_b);
    if (!geom_b) {
      TextUI::print_error("Failed to load geometry B");
      return false;
    }

    if (options.validate_mode || options.verbose) {
      auto validation = ErrorHandler::validate_geometry(*geom_b);
      if (!validation.valid) {
        TextUI::print_warning("Geometry B validation issues:");
        ErrorHandler::report_validation_error(validation);
        if (options.validate_mode) {
          return false;
        }
      } else if (options.verbose) {
        TextUI::print_success("Geometry B is valid");
        ErrorHandler::print_geometry_info(*geom_b);
      }
    }
  }

  auto start = std::chrono::high_resolution_clock::now();

  try {
    auto result = Operations::execute_operation(
        options.operation, options.op_arg, geom_a.get(), geom_b.get());

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;

    if (!result.has_value()) {
      std::string error_msg =
          "Unknown or unsupported operation: " + options.operation;
      if (!options.op_arg.empty()) {
        error_msg += " with argument: " + options.op_arg;
      }
      TextUI::print_error(error_msg);
      return false;
    }

    if (!options.quiet) {
      IO::print_result(result, options.format, options.precision);
    }

    if (options.show_time) {
      std::cout << TextUI::Colors::DIM << "Execution time: " << std::fixed
                << std::setprecision(6) << diff.count() << " seconds"
                << TextUI::Colors::RESET << "\n";
    }

    return true;
  } catch (const std::exception &e) {
    ErrorHandler::report_error(options.operation, e);
    return false;
  }
}
// NOLINTEND(readability-function-cognitive-complexity)

} // namespace

auto
main(int argc, char *argv[]) -> int
{
  auto options = parse_options(argc, argv);
  if (!options) {
    return EXIT_FAILURE;
  }

  if (options->help_requested) {
    print_help();
    return EXIT_SUCCESS;
  }

  if (options->version_requested) {
    print_version();
    return EXIT_SUCCESS;
  }

  if (options->list_operations) {
    print_operations_list();
    return EXIT_SUCCESS;
  }

  if (!options->help_operation.empty()) {
    if (Operations::print_operation_help(options->help_operation.c_str())) {
      return EXIT_SUCCESS;
    }
    TextUI::print_error(std::string("Unknown operation: ") +
                        options->help_operation);
    return EXIT_FAILURE;
  }

  // Handle case with operation
  if (!options->operation.empty()) {
    return handle_operation(*options) ? EXIT_SUCCESS : EXIT_FAILURE;
  }

  // Handle case without operation (geometry display/conversion)
  if (!options->source_a.empty()) {
    return handle_geometry_display(*options) ? EXIT_SUCCESS : EXIT_FAILURE;
  }

  // Should not reach here due to validation in parse_options
  TextUI::print_error("No operation or geometry specified");
  return EXIT_FAILURE;
}
