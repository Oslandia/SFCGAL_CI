# SFCGAL

[![GitLab pipeline status](https://gitlab.com/Oslandia/SFCGAL/badges/master/pipeline.svg)](https://gitlab.com/Oslandia/SFCGAL/-/commits/master)
[![GitHub pipeline status](https://github.com/Oslandia/SFCGAL_CI/actions/workflows/msys.yml/badge.svg)](https://github.com/Oslandia/SFCGAL_CI/actions?query=branch%3Amaster)
[![Cirrus pipeline status](https://api.cirrus-ci.com/github/Oslandia/SFCGAL_CI.svg)](http://cirrus-ci.com/github/Oslandia/SFCGAL_CI)

SFCGAL is a C++ wrapper library around [CGAL](http://www.cgal.org) with the aim of supporting ISO 191007:2013 and OGC Simple Features for 3D operations.

Please refer to the <a href="http://oslandia.gitlab.io/SFCGAL">project page</a> for an updated installation procedure.

## Licence

SFCGAL is provided under the following licence LGPL version 2 or later.

:warning: Note that the main dependency for SFCGAL is the CGAL library, and SFCGAL uses CGAL modules which are licenced as GPLv3+. Whenever you compile and distribute SCFGAL with the GPL-licenced CGAL, the full packaged result is automatically considered as GPL version 3 or later, due to GPL "viral" property. **If you link and distribute SFCGAL with another software package, be assured to fully understand the implications and check any legal and technical requirements implied by the licence**.
