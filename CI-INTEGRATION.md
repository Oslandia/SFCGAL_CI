# SFCGAL Continuous Integration Status

## Overview of CI Platforms

| Platform | Badge | Status | Environment | 
|:---------|:------|:-------|:------------|
| **GitLab** | [![GitLab pipeline status](https://gitlab.com/sfcgal/SFCGAL/badges/master/pipeline.svg)](https://gitlab.com/sfcgal/SFCGAL/-/commits/master) | 🚀 Main Pipeline | All environments |

## Detailed CI Status

### 🦊 GitLab Pipelines

#### 📊 Sonar Analysis
| Job | Badge | Environment |
|:----|:------|:------------|
| Sonar | [![Sonar](https://gitlab.com/sfcgal/SFCGAL/badges/master/pipeline.svg?job=sonar-build-test)](https://gitlab.com/sfcgal/SFCGAL/-/commits/master) | Debian |

#### 🧪 Platform Testing
| Job | Badge | Environment |
|:----|:------|:------------|
| Platform Tests | [![Tests](https://gitlab.com/sfcgal/SFCGAL/badges/master/pipeline.svg?job=platform-build-test)](https://gitlab.com/sfcgal/SFCGAL/-/commits/master) | Multiple Environments |

##### 🐧 Debian Environments
| Environment | Compiler | CGAL Version |
|:------------|:---------|:-------------|
| Debian Stable | Clang | CGAL 5.6.2 |
| Debian Testing | Clang | CGAL 5.6.2 |
| Debian Stable | Clang | CGAL 6.0.1 |
| Debian Testing | Clang | CGAL 6.0.1 |
| Debian Stable | GCC | CGAL 5.6.2 |
| Debian Testing | GCC | CGAL 5.6.2 |
| Debian Stable | GCC | CGAL 6.0.1 |
| Debian Testing | GCC | CGAL 6.0.1 |

##### 🎩 Fedora Environments
| Environment | Compiler | CGAL Version |
|:------------|:---------|:-------------|
| Fedora 40 | GCC | CGAL 5.6.2 |
| Fedora 40 | GCC | CGAL 6.0.1 |
| Fedora 40 | Clang | CGAL 5.6.2 |
| Fedora 40 | Clang | CGAL 6.0.1 |
| Fedora 41 | GCC | CGAL 5.6.2 |
| Fedora 41 | GCC | CGAL 6.0.1 |
| Fedora 41 | Clang | CGAL 5.6.2 |
| Fedora 41 | Clang | CGAL 6.0.1 |
| Fedora 42 | GCC | CGAL 5.6.2 |
| Fedora 42 | GCC | CGAL 6.0.1 |
| Fedora 42 | Clang | CGAL 5.6.2 |
| Fedora 42 | Clang | CGAL 6.0.1 |

##### 🦎 OpenSUSE Environment
| Environment | Compiler | CGAL Version |
|:------------|:---------|:-------------|
| OpenSUSE Tumbleweed | GCC | CGAL System* |

#### 🐳 Docker Builds
| Job | Badge | Environment |
|:----|:------|:------------|
| Docker | [![Docker](https://gitlab.com/sfcgal/SFCGAL/badges/master/pipeline.svg?job=docker)](https://gitlab.com/sfcgal/SFCGAL/-/commits/master) | Multiple Environments |

| Environment | Type |
|:------------|:-----|
| Docker Debian | 🐧 Linux |
| Docker Windows | 🪟 Windows |

### 😺 GitHub CI

#### 🔄 Cirrus CI
| Platform | Badge | Environment |
|:---------|:------|:------------|
| Cirrus CI | [![Cirrus pipeline status](https://api.cirrus-ci.com/github/Oslandia/SFCGAL_CI.svg)](http://cirrus-ci.com/github/Oslandia/SFCGAL_CI) | Multiple Environments |

| Environment | Compiler | CGAL Version |
|:------------|:---------|:-------------|
| 🍎 macOS | AppleClang | CGAL System* |
| 😈 FreeBSD 14.2 | Clang | CGAL 6.0.1 |
| 😈 FreeBSD 15.0 | Clang | CGAL 6.0.1 |

#### 🔄 GitHub Actions

##### 🔍 CodeQL Analysis
| Platform | Badge | Environment |
|:---------|:------|:------------|
| CodeQL | [![CodeQL](https://github.com/Oslandia/SFCGAL_CI/actions/workflows/codeql.yml/badge.svg?branch=master)](https://github.com/Oslandia/SFCGAL_CI/actions/workflows/codeql.yml) | Debian GCC CGAL 6.0.1 |

##### 📦 VCPKG Build
| Platform | Badge | Environment |
|:---------|:------|:------------|
| VCPKG | [![Build with vcpkg](https://github.com/Oslandia/SFCGAL_CI/actions/workflows/vcpkg.yml/badge.svg?branch=master)](https://github.com/Oslandia/SFCGAL_CI/actions/workflows/vcpkg.yml) | Multiple Environments |

| Environment | Compiler | CGAL Version |
|:------------|:---------|:-------------|
| 🐧 Ubuntu Latest | System | CGAL System* |
| 🍎 macOS Latest | System | CGAL System* |
| 🪟 Windows Latest | System | CGAL System* |

##### 🪟  MSYS2 Build
| Platform | Badge | Environment |
|:---------|:------|:------------|
| MSYS2 | [![MSYS2](https://github.com/Oslandia/SFCGAL_CI/actions/workflows/msys.yml/badge.svg?branch=master)](https://github.com/Oslandia/SFCGAL_CI/actions/workflows/msys.yml) | Multiple Environments |

| Environment | Compiler | CGAL Version |
|:------------|:---------|:-------------|
| mingw64 | System | CGAL System* |
| ucrt64 | System | CGAL System* |
| clang64 | System | CGAL System* |

##### 🌐 Cross Platform Actions
| Platform | Badge | Environment |
|:---------|:------|:------------|
| Cross Platform | [![CI](https://github.com/Oslandia/SFCGAL_CI/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/Oslandia/SFCGAL_CI/actions/workflows/ci.yml) | Multiple Environments |

| Environment | Compiler | CGAL Version |
|:------------|:---------|:-------------|
| 😈 FreeBSD 14.2 (quarterly) | Clang | CGAL 6.0.1 |
| 🚩 NetBSD 10.1 | GCC | CGAL 6.0.1 |
| 🐡 OpenBSD 7.7 | Clang | CGAL 6.0.1 |
| 🍎 macOS | AppleClang | CGAL System* |

\* *CGAL System refers to the CGAL version available through the system's package manager*
