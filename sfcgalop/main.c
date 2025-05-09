/**
 * SFCGALOP - Command-line tool for SFCGAL operations
 * 
 * This utility allows executing SFCGAL geometric operations from the command line.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <getopt.h>
#include "sfcgal_c.h"
#include "operations.h"
#include "io.h"
#include "timeus.h"
#include "util.h"

// Command-line options
typedef struct {
    char* sourceA;
    char* sourceB;
    OutputFormat format;
    int precision;
    bool quiet;
    bool showTime;
    bool verbose;
} Options;

// Function prototypes
void printHelp(void);
void printVersion(void);
Options parseOptions(int argc, char* argv[]);
void handleOperation(const char* opName, const char* opArg, Options* options);
bool operationRequiresArg(const char* opName);

/**
 * Print help information
 */
void printHelp(void) {
    printf("SFCGALOP - Command-line tool for SFCGAL operations\n\n");
    printf("Usage: sfcgalop [OPTION...] opName [opArg]\n");
    printf("  -a arg               source for A geometries (WKT, WKB, file, stdin, stdin.wkb)\n");
    printf("  -b arg               source for B geometries (WKT, WKB, file, stdin, stdin.wkb)\n");
    printf("  -f, --format arg     Output format (wkt, wkb, txt or geojson)\n");
    printf("  -p, --precision arg  Set number of decimal places in output coordinates\n");
    printf("  -q, --quiet          Disable result output\n");
    printf("  -t, --time           Print execution time\n");
    printf("  -v, --verbose        Verbose output\n");
    printf("  -h, --help           Print help\n");
    printf("\nSupported operations:\n");
    print_available_operations();
}

/**
 * Print version information
 */
void printVersion(void) {
    printf("SFCGAL version: %s\n", sfcgal_version());
    printf("SFCGAL full version: %s\n", sfcgal_full_version());
}

/**
 * Parse command-line options
 */
Options parseOptions(int argc, char* argv[]) {
    Options options = {
        .sourceA = NULL,
        .sourceB = NULL,
        .format = FORMAT_WKT,
        .precision = 6,
        .quiet = false,
        .showTime = false,
        .verbose = false
    };

    static struct option long_options[] = {
        {"format", required_argument, 0, 'f'},
        {"precision", required_argument, 0, 'p'},
        {"quiet", no_argument, 0, 'q'},
        {"time", no_argument, 0, 't'},
        {"verbose", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 'V'},
        {0, 0, 0, 0}
    };

    int c;
    int option_index = 0;

    while ((c = getopt_long(argc, argv, "a:b:f:p:qtvhV", long_options, &option_index)) != -1) {
        switch (c) {
            case 'a':
                options.sourceA = optarg;
                break;
            case 'b':
                options.sourceB = optarg;
                break;
            case 'f':
                if (strcmp(optarg, "wkt") == 0) {
                    options.format = FORMAT_WKT;
                } else if (strcmp(optarg, "wkb") == 0) {
                    options.format = FORMAT_WKB;
                } else if (strcmp(optarg, "txt") == 0) {
                    options.format = FORMAT_TXT;
                } else if (strcmp(optarg, "geojson") == 0) {
                    options.format = FORMAT_GEOJSON;
                } else {
                    fprintf(stderr, "Unknown format: %s\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'p':
                options.precision = atoi(optarg);
                if (options.precision < 0) {
                    fprintf(stderr, "Precision must be non-negative\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 'q':
                options.quiet = true;
                break;
            case 't':
                options.showTime = true;
                break;
            case 'v':
                options.verbose = true;
                break;
            case 'h':
                printHelp();
                exit(EXIT_SUCCESS);
                break;
            case 'V':
                printVersion();
                exit(EXIT_SUCCESS);
                break;
            case '?':
                // getopt_long already printed an error message
                exit(EXIT_FAILURE);
                break;
            default:
                abort();
        }
    }

    // Ensure we have enough arguments for operation
    if (optind >= argc) {
        fprintf(stderr, "Missing required argument: opName\n");
        exit(EXIT_FAILURE);
    }

    // If no sourceA is provided, assume stdin
    if (options.sourceA == NULL) {
        options.sourceA = "stdin";
    }

    return options;
}

/**
 * Check if an operation requires an argument
 */
bool operationRequiresArg(const char* opName) {
    const Operation* op = find_operation(opName);
    
    if (op == NULL) {
        // Operation not found, assume it doesn't require an argument
        return false;
    }
    
    return op->requires_arg;
}

/**
 * Execute the requested operation with the given options
 */
void handleOperation(const char* opName, const char* opArg, Options* options) {
    double start = 0.0, end = 0.0;
    double elapsed_time = 0.0;
    
    if (options->verbose) {
        printf("Operation: %s\n", opName);
        if (opArg) {
            printf("Argument: %s\n", opArg);
        }
        printf("Source A: %s\n", options->sourceA);
        if (options->sourceB) {
            printf("Source B: %s\n", options->sourceB);
        }
    }

    // Load geometry A
    sfcgal_geometry_t* geomA = NULL;
    if (load_geometry(options->sourceA, &geomA) != 0) {
        fprintf(stderr, "Failed to load geometry A from source: %s\n", options->sourceA);
        exit(EXIT_FAILURE);
    }

    // Load geometry B if needed
    sfcgal_geometry_t* geomB = NULL;
    if (options->sourceB != NULL) {
        if (load_geometry(options->sourceB, &geomB) != 0) {
            fprintf(stderr, "Failed to load geometry B from source: %s\n", options->sourceB);
            sfcgal_geometry_delete(geomA);
            exit(EXIT_FAILURE);
        }
    }

    // Start timing if requested
    if (options->showTime) {
        start = get_time_us();
    }

    // Execute the operation
    OperationResult result = execute_operation(opName, opArg, geomA, geomB);

    // End timing if requested
    if (options->showTime) {
        end = get_time_us();
        elapsed_time = end - start;
        printf("Execution time: %.6f microseconds\n", elapsed_time);
    }

    // Check for operation errors
    if (result.error) {
        fprintf(stderr, "Operation failed: %s\n", result.error_message);
        sfcgal_geometry_delete(geomA);
        if (geomB) sfcgal_geometry_delete(geomB);
        exit(EXIT_FAILURE);
    }

    // Output result if not quiet
    if (!options->quiet) {
        switch (result.type) {
            case RESULT_GEOMETRY:
                output_geometry(result.geometry_result, options->format, options->precision);
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
    sfcgal_geometry_delete(geomA);
    if (geomB) sfcgal_geometry_delete(geomB);
    
    // Free result resources if necessary
    if (result.type == RESULT_GEOMETRY && result.geometry_result) {
        sfcgal_geometry_delete(result.geometry_result);
    } else if (result.type == RESULT_TEXT && result.text_result) {
        free((void*)result.text_result);
    }
}

/**
 * Main entry point
 */
int main(int argc, char* argv[]) {
    // Initialize SFCGAL
    sfcgal_init();
    
    // Parse command-line options
    Options options = parseOptions(argc, argv);
    
    // Get operation name
    const char* opName = argv[optind];
    
    // Get operation argument if provided
    const char* opArg = NULL;
    if (optind + 1 < argc) {
        opArg = argv[optind + 1];
    } else if (operationRequiresArg(opName)) {
        // Error if operation requires an argument but none was provided
        fprintf(stderr, "Error: Operation '%s' requires an argument\n", opName);
        exit(EXIT_FAILURE);
    }
    
    // Set error handlers
    setup_error_handlers();
    
    // Execute the operation
    handleOperation(opName, opArg, &options);
    
    return EXIT_SUCCESS;
}
