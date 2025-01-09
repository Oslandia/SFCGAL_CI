# Installation

## Source code

SFCGAL source code is available on the corresponding [GitLab repository](https://gitlab.com/sfcgal/SFCGAL).

The `master` branch is the development branch and has a tag for each released version.

<!-- markdownlint-disable MD034 -->
Source code of the {{ get_project_version() }} release can be found with the tag v{{ get_project_version() }}. You can directly download a [zip](https://gitlab.com/sfcgal/SFCGAL/-/archive/v{{ get_project_version() }}/SFCGAL-v{{ get_project_version() }}.zip) or [tarball](https://gitlab.com/sfcgal/SFCGAL/-/archive/v{{ get_project_version() }}/SFCGAL-v{{ get_project_version() }}.tar.gz).
<!-- markdownlint-enable MD034 -->

## Supported platforms

SFCGAL has been successfully compiled and tested on the following platforms:

- Linux 64 bits with gcc and clang
- Windows with mingw
- macOS
- FreeBSD

## Requirements

- A C++ compiler, see above for supported platforms
- [CMake](https://cmake.org) version ≥ 2.8.6
- [CGAL](https://www.cgal.org) version ≥ 5.6
- [Boost](https://www.boost.org) version ≥ 1.74
- [MPFR](https://www.mpfr.org) version ≥ 2.2.1
- [GMP](https://gmplib.org) version ≥ 4.2

**Optional dependencies for viewer and 3D format export:**

- [OpenSceneGraph](https://openscenegraph.github.io/openscenegraph.io/) version ≥ 3.1
- [Qt5](https://contribute.qt-project.org/)

## Compilation

The compilation process is based on CMake. On Linux, run:

```bash
cmake -S . -B build && cmake --build build && sudo cmake --install build
```

You may specify dependencies locations with environment variables in case CMake doesn't find them automatically (see below).

CGAL uses a lot of templated constructions. Therefore, the building process may take a while.

Default building options should work out-of-the-box. You may want to fine-tune the compilation process with build options.

## Dependencies location (environment variables)

| Variable Name | Description         |
| ------------- | ------------------- |
| `GMP_DIR`     | GMP location         |
| `MPFR_DIR`    | MPFR location        |
| `CGAL_DIR`    | CGAL location        |
| `BOOST_ROOT`  | Boost location       |

## Build options (CMake options)

| Option Name              | Default   | Description                                                                                               |
| ------------------------ | --------- | --------------------------------------------------------------------------------------------------------- |
| `CMAKE_INSTALL_PREFIX`    | /usr/local | Specifies where SFCGAL will be installed when `make install` is invoked                                    |
| `CMAKE_BUILD_TYPE`        | Release   | Switches between a `Release` build for speed efficiency and a `Debug` build for development (slower, with assertions) |
| `SFCGAL_BUILD_VIEWER`     | OFF       | Turn to `ON` to build the viewer and conversion tools (viewer-SFCGAL)                                      |
| `SFCGAL_BUILD_TESTS`      | OFF       | Turn to `ON` to build unit and regression tests (run tests with `unit-test-SFCGAL`, etc.)                  |
| `SFCGAL_BUILD_EXAMPLES`   | OFF       | Turn to `ON` to build examples                                                                             |
| `SFCGAL_BUILD_BENCH`      | OFF       | Turn to `ON` to build benchmark tests                                                                      |
| `SFCGAL_WARNING_AS_ERROR` | OFF       | Turn to `ON` to convert build warnings into errors                                                         |
| `SFCGAL_BUILD_WITH_GPROF` | OFF       | Turn to `ON` to build with GNU gprof                                                                       |
| `SFCGAL_USE_STATIC_LIBS`  | OFF       | Turn to `ON` to compile as a static library                                                                |

## Installation with PostGIS

In a PostGIS source tree, an additional option for the `configure` script called `--with-sfcgal` is available. Set it to wherever your SFCGAL installation is. If you used the default directory `/usr/local` and if the `bin` sub-directory is present in your `PATH`, then SFCGAL should be detected without any additional options.
