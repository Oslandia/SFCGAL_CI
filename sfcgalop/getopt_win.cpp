/**
 * @file getopt_win.cpp
 * @brief Simple getopt implementation for Windows compatibility
 */

#ifdef _WIN32

  #include "getopt_win.hpp"
  #include <cstring>
  #include <iostream>

// Global variables for getopt state
int   optind = 1;       // Index of next argument to process
int   opterr = 1;       // Print error messages if non-zero
int   optopt = '?';     // Character that caused error
char *optarg = nullptr; // Pointer to option argument

// Internal state
namespace {
int  nextchar = 0;     // Position within current argument
bool optreset = false; // Reset flag
} // namespace

int
getopt(int argc, char *const argv[], const char *optstring)
{
  // Check for end of arguments
  if (optind >= argc || argv[optind] == nullptr || argv[optind][0] != '-' ||
      argv[optind][1] == '\0') {
    return -1;
  }

  // Handle "--" (end of options)
  if (strcmp(argv[optind], "--") == 0) {
    optind++;
    return -1;
  }

  // Get current character
  char opt = argv[optind][1 + nextchar];
  nextchar++;

  // Check if this is the last character in current argument
  if (argv[optind][1 + nextchar] == '\0') {
    optind++;
    nextchar = 0;
  }

  // Find option in optstring
  const char *opt_ptr = strchr(optstring, opt);
  if (opt_ptr == nullptr || opt == ':') {
    optopt = opt;
    if (opterr) {
      std::cerr << argv[0] << ": invalid option -- '" << static_cast<char>(opt)
                << "'\n";
    }
    return '?';
  }

  // Check if option requires an argument
  if (opt_ptr[1] == ':') {
    if (nextchar == 0) {
      // Argument is in next argv element
      if (optind >= argc) {
        optopt = opt;
        if (opterr) {
          std::cerr << argv[0] << ": option requires an argument -- '"
                    << static_cast<char>(opt) << "'\n";
        }
        return optstring[0] == ':' ? ':' : '?';
      }
      optarg = argv[optind++];
    } else {
      // Argument is remainder of current argv element
      optarg = &argv[optind][1 + nextchar];
      optind++;
      nextchar = 0;
    }
  } else {
    optarg = nullptr;
  }

  return opt;
}

int
getopt_long(int argc, char *const argv[], const char *optstring,
            const struct option *longopts, int *longindex)
{
  // Check for end of arguments
  if (optind >= argc || argv[optind] == nullptr || argv[optind][0] != '-') {
    return -1;
  }

  // Handle "--" (end of options)
  if (strcmp(argv[optind], "--") == 0) {
    optind++;
    return -1;
  }

  // Handle long options (--option)
  if (argv[optind][0] == '-' && argv[optind][1] == '-') {
    const char *name     = &argv[optind][2];
    const char *eq       = strchr(name, '=');
    size_t      name_len = eq ? static_cast<size_t>(eq - name) : strlen(name);

    // Search for matching long option
    for (int i = 0; longopts[i].name != nullptr; i++) {
      if (strncmp(name, longopts[i].name, name_len) == 0 &&
          strlen(longopts[i].name) == name_len) {
        if (longindex) {
          *longindex = i;
        }

        optind++;

        // Handle argument
        if (longopts[i].has_arg == required_argument) {
          if (eq) {
            optarg = const_cast<char *>(eq + 1);
          } else if (optind < argc) {
            optarg = argv[optind++];
          } else {
            if (opterr) {
              std::cerr << argv[0] << ": option '--" << longopts[i].name
                        << "' requires an argument\n";
            }
            return '?';
          }
        } else if (longopts[i].has_arg == optional_argument) {
          optarg = eq ? const_cast<char *>(eq + 1) : nullptr;
        } else {
          if (eq) {
            if (opterr) {
              std::cerr << argv[0] << ": option '--" << longopts[i].name
                        << "' doesn't allow an argument\n";
            }
            return '?';
          }
          optarg = nullptr;
        }

        // Handle flag
        if (longopts[i].flag) {
          *longopts[i].flag = longopts[i].val;
          return 0;
        } else {
          return longopts[i].val;
        }
      }
    }

    if (opterr) {
      std::cerr << argv[0] << ": unrecognized option '--" << name << "'\n";
    }
    return '?';
  }

  // Handle short options using getopt
  return getopt(argc, argv, optstring);
}

#endif // _WIN32