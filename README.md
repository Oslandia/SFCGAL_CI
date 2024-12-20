# SFCGAL

[![GitLab pipeline status](https://gitlab.com/Oslandia/SFCGAL/badges/master/pipeline.svg)](https://gitlab.com/Oslandia/SFCGAL/-/commits/master)
[![GitHub pipeline status](https://github.com/Oslandia/SFCGAL_CI/actions/workflows/msys.yml/badge.svg)](https://github.com/Oslandia/SFCGAL_CI/actions?query=branch%3Amaster)
[![Cirrus pipeline status](https://api.cirrus-ci.com/github/Oslandia/SFCGAL_CI.svg)](http://cirrus-ci.com/github/Oslandia/SFCGAL_CI)

SFCGAL is a C++ wrapper library around [CGAL](http://www.cgal.org) with the aim of supporting ISO 19107:2019 and OGC Simple Features for 3D operations.

Please refer to the <a href="http://sfcgal.gitlab.io/SFCGAL">project page</a> for an updated installation procedure.

[![Packaging status](https://repology.org/badge/vertical-allrepos/sfcgal.svg)](https://repology.org/project/sfcgal/versions)

## Community Resources

* Website: <https://sfcgal.org>
* **git** repository: <https://gitlab.com/sfcgal/SFCGAL>
* **#sfcgal** chat channel (all bridged):
  * IRC: irc://irc.libera.chat/#sfcgal (<https://kiwiirc.com/nextclient/irc.libera.chat/#sfcgal>)

## Build/Install

See the [INSTALL](https://sfcgal.gitlab.io/SFCGAL/installation/) page.

## Reference Docs

* [C API](https://sfcgal.gitlab.io/SFCGAL/API/sfcgal__c_8h/)
* [C++ API](https://sfcgal.gitlab.io/SFCGAL/API/links/)

## Client Applications

### Using the C interface

SFCGAL promises long-term stability of the C API. In general, successive releases
of the C API may add new functions but will not remove or change existing types
or function signatures. The C library uses the C++ interface, but the C library
follows normal ABI-change-sensitive versioning, so programs that link only
against the C library should work without relinking when SFCGAL is upgraded. For
this reason, it is recommended to use the C API for software that is intended
to be dynamically linked to a system install of SFCGAL.

The `sfcgal-config` program can be used to determine appropriate compiler and
linker flags for building against the C library:

    CFLAGS += `sfcgal-config --cflags`
    LDFLAGS += `sfcgal-config --ldflags --libs`

All functionality of the C API is available through the `sfcgal_c.h` header file.

Documentation for the C API is provided via comments in the `sfcgal_c.h` header
file. C API usage examples can be found in the SFCGAL unit tests and in the
source code of software that uses SFCGAL, such as PostGIS and the PySFCGAL package
for Python.

### Using other languages

SFCGAL has bindings in many languages, see the [bindings page](https://sfcgal.gitlab.io/SFCGAL/development/#bindings).

## Documentation

API documentation can be generated using Doxygen. Documentation is not included
in the default build. To build the documentation see the [Development](https://sfcgal.gitlab.io/SFCGAL/development/) page.

## Style

To format your code into the desired style, use the `clang-format` tools.

It can be automatically called with a series of *git hooks* (thanks to the [pre-commit](https://pre-commit.com/)) tool). To install them:  

```bash
pre-commit install
```

To commit without the *git hooks*, add the `--no-verify` option to the `git commit` command.

## Testing

See documentation in [Development](https://sfcgal.gitlab.io/SFCGAL/development/) page.

## Licence

SFCGAL is provided under the following licence LGPL version 2 or later.

:warning: Note that the main dependency for SFCGAL is the CGAL library, and SFCGAL uses CGAL modules which are licenced as GPLv3+. Whenever you compile and distribute SFCGAL with the GPL-licenced CGAL, the full packaged result is automatically considered as GPL version 3 or later, due to GPL "viral" property. **If you link and distribute SFCGAL with another software package, be assured to fully understand the implications and check any legal and technical requirements implied by the licence**.
