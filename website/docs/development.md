# Development

The project uses git and [Gitlab](https://about.gitlab.com) to manage the source code, bug reports, and issues.

You are free to fork the project and are welcome to propose evolutions. The preferred way of submitting patches is through Gitlab's merge request mechanism.

Bugs can be reported via the [Issues](https://gitlab.com/sfcgal/SFCGAL/-/issues) section of the Gitlab project.

## Branches

The `master` branch of the repository represents the current development version of SFCGAL. Each commit on this branch is enforced not to break the compilation and unit/regression tests.

For each public release, a new tag is created on this branch.

When different major versions are to be maintained, the most recent major version is hosted on the `master` branch, and older versions have their own branch where fixes can be backported if needed.

In conformance with the git way of development, each new feature is developed on its own temporary branch.

Version numbers follow the [Semantic Versioning 2.0.0](https://semver.org) policy and are tagged as **x.y.z**, where:

- **x** is the major version number. It changes in case of API breaks or major redesign.
- **y** is the minor version number. It changes when new functionalities are added without a major API break.
- **z** is the patch version number. It changes when bugs or packaging issues are resolved.

## Building

The classic way of building with CMake is to first create a build directory and use a CMake client (ccmake or cmake-gui) to configure the build:

```bash
mkdir build_debug
cd build_debug
ccmake -DCMAKE_BUILD_TYPE=Debug .. && make -j 2
```

Environment variables and build options are listed in the installation [section](./installation.md).

## Tests

SFCGAL comes with different layers of tests. The SFCGAL_BUILD_TESTS CMake option allows you to build these tests.

- **Unit tests**: Test each feature independently with hand-crafted datasets that must cover every possible case. Run it with unit-test-SFCGAL.
- **Regression tests**: Consist of function calls with real or near-real datasets. Run it with standalone-regress-test-SFCGAL.
- **Garden test**: Inspired by the eponymous PostGIS test, it ensures that all combinations of parameters are acceptable for each function (i.e., no crashes occur). Run it with garden-test-SFCGAL.
- **Style test**: Ensures the code is properly formatted before committing. Run the script test-style-SFCGAL.sh.
- **Benchmark tests**: Measure the processing speed of SFCGAL algorithms. Compilation is available through the SFCGAL_BUILD_BENCH CMake option.

## Valgrind

[Valgrind](https://valgrind.org/) can automatically detect many memory management and threading bugs, and profile programs in detail. Since SFCGAL uses a CGAL Kernel that depends on floating-point rounding modes that are not supported by the current version of Valgrind, a patched version is necessary:

```sh title "Install Valgrind"
git clone https://github.com/trast/valgrind.git
cd valgrind/
git clone https://github.com/trast/valgrind-VEX.git VEX
./autogen.sh
./configure && make -j 8 && sudo make install
```

## Documentation

The SFCGAL documentation is primarily written in [Doxygen](https://www.doxygen.nl/). If Doxygen is available, you can run `make doc` to generate the code documentation.

## Releasing

To release a new SFCGAL version:

- Update the file NEWS in the master branch.
- Change the version number in the root CMakeLists.txt in the master branch.
- Release from Gitlab (Gitlab will create the tag for you).
- Update links to the new version in the gh-pages branch.

## Bindings

SFCGAL can be accessed from a variety of other languages via bindings to the library.

### Python (official)

[PySFCGAL](https://gitlab.com/sfcgal/pysfcgal) by the SFCGAL Team.

### Rust

[sfcgal-rs](https://github.com/mthh/sfcgal-rs) and [sfcgal-sys](https://github.com/mthh/sfcgal-sys) by Matthieu Viry.

### Nim

[sfcgal.nim](https://gitlab.com/lbartoletti/sfcgal.nim) by Loïc Bartoletti.
