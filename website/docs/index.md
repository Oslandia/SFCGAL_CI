# SFCGAL

![logo SFCGAL](assets/img/logo.svg)

## Downloads

<!-- markdownlint-disable MD034 -->
[Download this project as a .zip file](https://gitlab.com/sfcgal/SFCGAL/-/archive/v{{ get_project_version() }}/SFCGAL-v{{ get_project_version() }}.zip){ .md-button .md-button--primary }
[Download this project as a tar.gz file](https://gitlab.com/sfcgal/SFCGAL/-/archive/v{{ get_project_version() }}/SFCGAL-v{{ get_project_version() }}.tar.gz){ .md-button .md-button--primary }
<!-- markdownlint-enable MD034 -->

## About

SFCGAL is a C++ wrapper library around [CGAL](http://www.cgal.org) with the aim of supporting ISO 19107:2013 and [OGC Simple Features Access 1.2](https://www.ogc.org/publications/standard/sfa/) for 3D operations.

SFCGAL provides standard compliant geometry types and operations. PostGIS uses the C API, to expose some SFCGAL's functions in spatial databases (cf. [PostGIS manual](https://postgis.net/docs/reference_sfcgal.html)).

[PySFCGAL](https://sfcgal.gitlab.io/pysfcgal) is a Python wrapper for SFCGAL, enabling seamless use of the library’s 3D geometric functions and operations directly within Python scripts or applications. It is accessible on the [Python Package Index (PyPI)](https://pypi.org/project/PySFCGAL/) (cf. [PySFCGAL manual](https://sfcgal.gitlab.io/pysfcgal)).

### Supported Geometry Types

Geometry coordinates have an exact rational number representation and can be either 2D or 3D. Among supported geometry types are:

- Points
- LineStrings
- Polygons
- TriangulatedSurfaces
- PolyhedralSurfaces
- GeometryCollections
- Solids

### Supported Operations

Supported operations include:

- WKT reading and writing with exact rational number representation for coordinates
- Intersection operations and predicates
- Convex hull computation
- Tessellation
- Extrusion
- Area and distance computation
- Minkowski sums
- Contour offsets
- Straight skeleton generations

## License

SFCGAL is distributed under the terms of the [GNU Lesser General Public License 2+](http://www.gnu.org/licenses/old-licenses/lgpl-2.0.html).

!!! warning License

**Note** that the main dependency for **SFCGAL** is the **CGAL** library, and SFCGAL uses CGAL modules that are licensed under **GPLv3+**. Whenever you compile and distribute **SFCGAL** with the GPL-licensed version of **CGAL**, the entire packaged result is automatically considered to be under **GPL version 3 or later**, due to the GPL’s “viral” property.

If you link and distribute **SFCGAL** with another software package, make sure you fully understand the implications and verify any legal and technical requirements implied by the license.
