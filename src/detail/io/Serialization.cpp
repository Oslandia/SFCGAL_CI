// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/io/Serialization.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"

namespace SFCGAL::io {

BinarySerializer::BinarySerializer(std::ostream &ostr)
    : boost::archive::binary_oarchive(ostr)
{
  using namespace SFCGAL;
  register_type<Point>();
  register_type<LineString>();
  register_type<Triangle>();
  register_type<Polygon>();
  register_type<TriangulatedSurface>();
  register_type<PolyhedralSurface>();
  register_type<Solid>();
  register_type<GeometryCollection>();
  register_type<MultiPoint>();
  register_type<MultiLineString>();
  register_type<MultiPolygon>();
  register_type<MultiSolid>();
}

BinaryUnserializer::BinaryUnserializer(std::istream &ostr)
    : boost::archive::binary_iarchive(ostr)
{
  using namespace SFCGAL;
  register_type<Point>();
  register_type<LineString>();
  register_type<Triangle>();
  register_type<Polygon>();
  register_type<TriangulatedSurface>();
  register_type<PolyhedralSurface>();
  register_type<Solid>();
  register_type<GeometryCollection>();
  register_type<MultiPoint>();
  register_type<MultiLineString>();
  register_type<MultiPolygon>();
  register_type<MultiSolid>();
}

auto
writeBinaryGeometry(const Geometry &g) -> std::string
{
  std::ostringstream ostr;
  BinarySerializer   arc(ostr);
  // use the pointer version to force dynamic type writing
  const Geometry *pg = &g;
  arc << pg;
  return ostr.str();
}

auto
writeBinaryPrepared(const PreparedGeometry &g) -> std::string
{
  std::ostringstream      ostr;
  BinarySerializer        arc(ostr);
  const PreparedGeometry *pg = &g;
  arc << pg;
  return ostr.str();
}

auto
readBinaryGeometry(const std::string &str) -> std::unique_ptr<Geometry>
{
  std::istringstream istr(str);
  BinaryUnserializer iarc(istr);
  Geometry          *g = nullptr;
  iarc >> g;
  return std::unique_ptr<Geometry>(g);
}

auto
readBinaryPrepared(const std::string &str) -> std::unique_ptr<PreparedGeometry>
{
  std::istringstream istr(str);
  BinaryUnserializer iarc(istr);
  PreparedGeometry  *pg = nullptr;
  iarc >> pg;
  return std::unique_ptr<PreparedGeometry>(pg);
}
} // namespace SFCGAL::io
namespace boost::serialization {

void
save(boost::archive::text_oarchive &ar, const CGAL::Gmpz &z,
     const unsigned int /*version*/)
{
  std::ostringstream ostr;
  ostr << z;
  std::string const str = ostr.str();
  ar << str;
}

// specialization for binary archives
void
save(boost::archive::binary_oarchive &ar, const CGAL::Gmpz &z,
     const unsigned int /* version*/)
{
  const mpz_t  &mpz  = z.mpz();
  int32_t const size = mpz->_mp_size;
  ar & size;
  uint32_t const rsize = size >= 0 ? size : -size;

  for (uint32_t i = 0; i < rsize; ++i) {
    ar & mpz->_mp_d[i];
  }
}

void
load(boost::archive::text_iarchive &ar, CGAL::Gmpz &z,
     const unsigned int /*version*/)
{
  std::string line;
  ar >> line;
  std::istringstream istr(line);
  istr >> z;
}

void
load(boost::archive::binary_iarchive &ar, CGAL::Gmpz &z,
     const unsigned int /*version*/)
{
  int32_t  size  = 0;
  uint32_t rsize = 0;
  mpz_t   &mpz   = z.mpz();
  ar & size;
  rsize         = size >= 0 ? size : -size;
  mpz->_mp_size = size;
  _mpz_realloc(mpz, rsize);
  uint32_t i = 0;

  for (i = 0; i < rsize; ++i) {
    ar & mpz->_mp_d[i];
  }
}

#ifdef CGAL_USE_GMPXX
void
save(boost::archive::text_oarchive &ar, const mpz_class &z,
     const unsigned int /*version*/)
{
  std::ostringstream ostr;
  ostr << z;
  std::string const str = ostr.str();
  ar << str;
}

// specialization for binary archives
void
save(boost::archive::binary_oarchive &ar, const mpz_class &z,
     const unsigned int /* version*/)
{
  mpz_srcptr    mpz  = z.get_mpz_t();
  int32_t const size = mpz->_mp_size;
  ar & size;
  uint32_t const rsize = size >= 0 ? size : -size;

  for (uint32_t i = 0; i < rsize; ++i) {
    ar & mpz->_mp_d[i];
  }
}

void
load(boost::archive::text_iarchive &ar, mpz_class &z,
     const unsigned int /*version*/)
{
  std::string line;
  ar >> line;
  std::istringstream istr(line);
  istr >> z;
}

void
load(boost::archive::binary_iarchive &ar, mpz_class &z,
     const unsigned int /*version*/)
{
  int32_t  size  = 0;
  uint32_t rsize = 0;
  mpz_ptr  mpz   = z.get_mpz_t();
  ar & size;
  rsize         = size >= 0 ? size : -size;
  mpz->_mp_size = size;
  _mpz_realloc(mpz, rsize);
  uint32_t i = 0;

  for (i = 0; i < rsize; ++i) {
    ar & mpz->_mp_d[i];
  }
}
#endif

} // namespace boost::serialization
