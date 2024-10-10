# SFCGAL

![](assets/img/logo.svg)

## Downloads

[Download this project as a .zip file](https://gitlab.com/SFCGAL/SFCGAL/-/archive/v{{ get_project_version() }}/SFCGAL-v{{ get_project_version() }}.zip){ .md-button .md-button--primary }
[Download this project as a tar.gz file](https://gitlab.com/SFCGAL/SFCGAL/-/archive/v{{ get_project_version() }}/SFCGAL-v{{ get_project_version() }}.tar.gz){ .md-button .md-button--primary }

## About

SFCGAL is a C++ wrapper library around [CGAL](http://www.cgal.org) with the aim of supporting ISO 19107:2013 and [OGC Simple Features Access 1.2](http://www.opengeospatial.org/standards/sfa) for 3D operations.

SFCGAL provides standard compliant geometry types and operations. PostGIS uses the C API, to expose some SFCGAL's functions in spatial databases (cf. [PostGIS manual](https://postgis.net/docs/reference_sfcgal.html)).

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

⚠ **Note** that the main dependency for SFCGAL is the CGAL library, and SFCGAL uses CGAL modules which are licenced as GPLv3+. Whenever you compile and distribute SCFGAL with the GPL-licenced CGAL, the full packaged result is automatically considered as GPL version 3 or later, due to GPL "viral" property.

If you link and distribute SFCGAL with another software package, be assured to fully understand the implications and check any legal and technical requirements implied by the licence.
