/**
 * main.c - Command-line tool for SFCGAL operations
 *
 * This utility allows executing SFCGAL geometric operations from the command
 * line.
 */

#include "io.h"
#include "operations/operations.h"
#include "safe_string.h"
#include "sfcgal_c.h"
#include "timeus.h"
#include "util.h"
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PROGRAM_NAME "sfcgalop"

/**
 * Command-line options structure
 */
typedef struct {
  char        *source_a;          // Source for A geometries
  char        *source_b;          // Source for B geometries (optional)
  OutputFormat format;            // Output format
  int          precision;         // Decimal precision for output
  bool         quiet;             // Suppress result output
  bool         show_time;         // Show execution time
  bool         verbose;           // Enable verbose output
  bool         help_requested;    // Help was requested
  bool         version_requested; // Version info was requested
  char        *operation;         // Operation name
  char        *op_arg;            // Operation argument
} Options;

/**
 * Print help information
 */
static void
print_help(void)
{
  printf("Usage: %s [OPTION...] operation [operation_arg]\n\n", PROGRAM_NAME);
  printf("SFCGALOP - Command-line tool for SFCGAL operations\n\n");
  printf("Options:\n");
  printf(
      "  -a, --geom-a=ARG        Source for A geometries (WKT, WKB, TXT, file, "
      "stdin, stdin.wkb)\n");
  printf(
      "  -b, --geom-b=ARG        Source for B geometries (WKT, WKB, TXT, file, "
      "stdin, stdin.wkb)\n");
  printf("  -f, --format=ARG        Output format:\n");
  printf(
      "                            wkt     - WKT with decimal coordinates\n");
  printf("                            wkb     - WKB binary (hex encoded)\n");
  printf("                            txt     - SFCGAL native format (exact "
         "fractions)\n");
  printf("                            geojson - GeoJSON format\n");
  printf("  -p, --precision=N       Set number of decimal places in output "
         "coordinates (wkt/geojson only)\n");
  printf("  -q, --quiet             Disable result output\n");
  printf("  -t, --time              Print execution time\n");
  printf("  -v, --verbose           Verbose output\n");
  printf("  -h, --help              Print help\n");
  printf("  -V, --version           Print version information\n");
  printf("      --help-op=OPERATION Print help for specific operation\n");
  printf("\nAvailable operations:\n");
  print_available_operations();
  printf("\nOutput format examples:\n");
  printf("  wkt:     POINT(0.5 1.25)\n");
  printf("  txt:     POINT(1/2 5/4)\n");
  printf("  wkb:     0101000000000000000000E03F0000000000001440\n");
  printf("\nExamples:\n");
  printf("  %s -a \"POINT(0 0)\" -b \"POINT(1 1)\" distance\n", PROGRAM_NAME);
  printf("  %s -a \"POLYGON((0 0, 0 10, 10 10, 10 0, 0 0))\" area\n",
         PROGRAM_NAME);
  printf("  %s -a \"POLYGON((0 0, 0 10, 10 10, 10 0, 0 0))\" offset_polygon "
         "radius=5\n",
         PROGRAM_NAME);
  printf("  %s -a input.wkt -f txt convexhull\n", PROGRAM_NAME);
  printf("\nOperation arguments can be specified positionally or by name "
         "(name=value).\n");
  printf("For example: offset_polygon radius=5 or offset_polygon 5\n");
}

/**
 * Print version information
 */
static void
print_version(void)
{
  printf("%s version\n", PROGRAM_NAME);
  printf("SFCGAL version: %s\n", sfcgal_version());
  printf("SFCGAL full version: %s\n", sfcgal_full_version());
}

/**
 * Initialize default options
 *
 * @param options Pointer to options structure
 */
static void
init_options(Options *options)
{
  if (!options) {
    return;
  }

  options->source_a          = NULL;
  options->source_b          = NULL;
  options->format            = FORMAT_WKT;
  options->precision         = 6;
  options->quiet             = false;
  options->show_time         = false;
  options->verbose           = false;
  options->help_requested    = false;
  options->version_requested = false;
  options->operation         = NULL;
  options->op_arg            = NULL;
}

/**
 * Parse command-line options
 *
 * @param argc Number of arguments
 * @param argv Array of argument strings
 * @param options Pointer to options structure
 * @return true on success, false on failure
 */
static bool
parse_options(int argc, char *argv[], Options *options)
{
  if (!argv || !options || argc < 0 || argc > 10) {
    return false;
  }

  // Initialize options with defaults
  init_options(options);

  // If no arguments, show help
  if (argc == 1) {
    options->help_requested = true;
    return true;
  }

  // Long options
  static struct option long_options[] = {
      {"geom-a", required_argument, 0, 'a'},
      {"geom-b", required_argument, 0, 'b'},
      {"format", required_argument, 0, 'f'},
      {"precision", required_argument, 0, 'p'},
      {"quiet", no_argument, 0, 'q'},
      {"time", no_argument, 0, 't'},
      {"verbose", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'V'},
      {"help-op", required_argument, 0, 'H'},
      {0, 0, 0, 0}};

  int  c;
  int  option_index = 0;
  bool error        = false;

  if (argc > ARG_MAX || argc < 0) {
    fprintf(stderr, "Invalid argument count\n");
    return EXIT_FAILURE;
  }

  // Validate argv strings length to prevent buffer overflow
  for (int i = 0; i < argc; i++) {
    if (argv[i] && safe_strlen(argv[i], SAFE_MAX_STRING_LENGTH) == SIZE_MAX) {
      fprintf(stderr, "Argument %d too long\n", i);
      return false;
    }
  }

  // Process options
  while ((c = getopt_long(argc, argv, // flawfinder: ignore
                          "a:b:f:p:qtvhV", long_options, &option_index)) !=
         -1) {
    switch (c) {
    case 'a':
      options->source_a = optarg;
      break;
    case 'b':
      options->source_b = optarg;
      break;
    case 'f':
      if (!parse_output_format(optarg, &options->format)) {
        fprintf(stderr, "Unknown format: %s\n", optarg);
        fprintf(stderr, "Valid formats are:\n");
        fprintf(stderr, "  wkt     - WKT with decimal coordinates\n");
        fprintf(stderr, "  wkb     - WKB binary (hex encoded)\n");
        fprintf(stderr, "  txt     - SFCGAL native format (exact fractions)\n");
        fprintf(stderr, "  geojson - GeoJSON format\n");
        error = true;
      }
      break;
    case 'p':
      errno = 0;
      char *endptr;
      long  precision = strtol(optarg, &endptr, 10);

      if (errno != 0 || *endptr != '\0' || precision < 0 ||
          precision > INT_MAX) {
        fprintf(stderr, "Invalid precision: %s\n", optarg);
        fprintf(stderr, "Precision must be a non-negative integer\n");
        error = true;
      } else {
        options->precision = (int)precision;
      }
      break;
    case 'q':
      options->quiet = true;
      break;
    case 't':
      options->show_time = true;
      break;
    case 'v':
      options->verbose = true;
      break;
    case 'h':
      options->help_requested = true;
      break;
    case 'V':
      options->version_requested = true;
      break;
    case 'H':
      if (print_operation_help(optarg)) {
        exit(EXIT_SUCCESS);
      } else {
        exit(EXIT_FAILURE);
      }
      break;
    case '?':
      // getopt_long already printed an error message
      error = true;
      break;
    default:
      fprintf(stderr, "Unexpected getopt return value: %d\n", c);
      error = true;
      break;
    }

    if (error) {
      break;
    }
  }

  if (error) {
    return false;
  }

  // Check remaining arguments (operation and optional operation argument)
  if (optind < argc) {
    options->operation = argv[optind++];

    // Check for operation argument
    if (optind < argc) {
      options->op_arg = argv[optind++];

      // Check for unexpected additional arguments
      if (optind < argc) {
        fprintf(stderr, "Unexpected extra arguments starting with: %s\n",
                argv[optind]);
        return false;
      }
    } else if (operation_requires_arg(options->operation)) {
      // If operation requires an argument but none was provided
      fprintf(stderr, "Error: Operation '%s' requires arguments\n",
              options->operation);
      fprintf(stderr, "Run '%s --help-op=%s' for more information\n",
              PROGRAM_NAME, options->operation);
      return false;
    }
  } else if (!options->help_requested && !options->version_requested) {
    // No operation provided and not asking for help or version
    fprintf(stderr, "Missing required argument: operation\n");
    return false;
  }

  // If no source_a is provided, assume stdin
  if (!options->source_a) {
    options->source_a = safe_strdup("stdin", SAFE_MAX_STRING_LENGTH);
  }

  return true;
}

/**
 * Execute the requested operation with the given options
 *
 * @param options Command-line options
 * @return true on success, false on failure
 */
static bool
handle_operation(const Options *options)
{
  if (!options || !options->operation) {
    return false;
  }

  if (options->verbose) {
    printf("Operation: %s\n", options->operation);
    if (options->op_arg) {
      printf("Arguments: %s\n", options->op_arg);
    }
    printf("Source A: %s\n", options->source_a);
    if (options->source_b) {
      printf("Source B: %s\n", options->source_b);
    }
  }

  double start = 0.0, end = 0.0;

  // Load geometry A
  sfcgal_geometry_t *geom_a = NULL;
  if (!load_geometry(options->source_a, &geom_a)) {
    fprintf(stderr, "Failed to load geometry A from source: %s\n",
            options->source_a);
    return false;
  }

  // Load geometry B if needed
  sfcgal_geometry_t *geom_b = NULL;
  if (options->source_b) {
    if (!load_geometry(options->source_b, &geom_b)) {
      fprintf(stderr, "Failed to load geometry B from source: %s\n",
              options->source_b);
      sfcgal_geometry_delete(geom_a);
      return false;
    }
  }

  // Start timing if requested
  if (options->show_time) {
    start = get_time_us();
  }

  // Execute the operation
  OperationResult result =
      execute_operation(options->operation, options->op_arg, geom_a, geom_b);

  // End timing if requested
  if (options->show_time) {
    end = get_time_us();
    printf("Execution time: %.6f microseconds\n", end - start);
  }

  // Check for operation errors
  if (result.error) {
    fprintf(stderr, "Operation failed: %s\n", result.error_message);
    sfcgal_geometry_delete(geom_a);
    if (geom_b)
      sfcgal_geometry_delete(geom_b);
    return false;
  }

  // Output result if not quiet
  if (!options->quiet) {
    switch (result.type) {
    case RESULT_GEOMETRY:
      if (!output_geometry(result.geometry_result, options->format,
                           options->precision)) {
        fprintf(stderr, "Error outputting geometry\n");
      }
      break;
    case RESULT_BOOLEAN:
      printf("%s\n", result.boolean_result ? "true" : "false");
      break;
    case RESULT_NUMERIC:
      printf("%.10g\n", result.numeric_result);
      break;
    case RESULT_TEXT:
      printf("%s\n", result.text_result);
      break;
    default:
      fprintf(stderr, "Unknown result type\n");
      break;
    }
  }

  // Free resources
  sfcgal_geometry_delete(geom_a);
  if (geom_b)
    sfcgal_geometry_delete(geom_b);

  // Free result resources if necessary
  if (result.type == RESULT_GEOMETRY && result.geometry_result) {
    sfcgal_geometry_delete(result.geometry_result);
  } else if (result.type == RESULT_TEXT && result.text_result) {
    free(result.text_result);
  }

  return true;
}

/**
 * Main entry point
 */
int
main(int argc, char *argv[])
{
  // Initialize SFCGAL
  sfcgal_init();

  // Set up error handlers
  setup_error_handlers();

  // Parse command-line options
  Options options;
  if (!parse_options(argc, argv, &options)) {
    fprintf(stderr, "Try '%s --help' for more information.\n", PROGRAM_NAME);
    return EXIT_FAILURE;
  }

  // Handle special commands
  if (options.help_requested) {
    print_help();
    return EXIT_SUCCESS;
  }

  if (options.version_requested) {
    print_version();
    return EXIT_SUCCESS;
  }

  // Execute the operation
  if (!handle_operation(&options)) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
