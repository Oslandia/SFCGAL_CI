/**
 * @file getopt_win.hpp
 * @brief Simple getopt implementation for Windows compatibility
 *
 * This file provides a minimal getopt/getopt_long implementation for Windows
 * platforms that don't have a native getopt.h. This allows sfcgalop to work
 * on Windows systems.
 */

#ifndef SFCGALOP_GETOPT_WIN_HPP
#define SFCGALOP_GETOPT_WIN_HPP

#ifdef _WIN32

  #include <string>
  #include <vector>

// Global variables used by getopt (standard interface)
extern int   optind; // Index of next argument to process
extern int   opterr; // Print error messages if non-zero
extern int   optopt; // Character that caused error
extern char *optarg; // Pointer to option argument

// Option structure for getopt_long
struct option {
  const char *name;    // Long option name
  int         has_arg; // Argument requirement
  int        *flag;    // Flag to set (if not null)
  int         val;     // Value to return (if flag is null)
};

// Argument requirements for options
constexpr int no_argument       = 0;
constexpr int required_argument = 1;
constexpr int optional_argument = 2;

/**
 * @brief Parse short options from command line arguments.
 *
 * Simple implementation of POSIX getopt() for Windows compatibility.
 * Parses short options of the form "-x" or "-x value".
 *
 * @param argc Number of command line arguments
 * @param argv Array of command line arguments
 * @param optstring String of valid option characters, with ':' after options
 * that require arguments
 * @return Option character on success, '?' on unknown option, ':' on missing
 * argument, -1 when done
 */
int
getopt(int argc, char *const argv[], const char *optstring);

/**
 * @brief Parse long options from command line arguments.
 *
 * Simple implementation of getopt_long() for Windows compatibility.
 * Parses both short options ("-x") and long options ("--option").
 *
 * @param argc Number of command line arguments
 * @param argv Array of command line arguments
 * @param optstring String of valid short option characters
 * @param longopts Array of long option structures
 * @param longindex Pointer to store index of matched long option (can be null)
 * @return Option character on success, '?' on unknown option, ':' on missing
 * argument, -1 when done
 */
int
getopt_long(int argc, char *const argv[], const char *optstring,
            const struct option *longopts, int *longindex);

#endif // _WIN32

#endif // SFCGALOP_GETOPT_WIN_HPP