// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_SERIALIZATION_
#define SFCGAL_SERIALIZATION_

#include "SFCGAL/config.h"

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/PreparedGeometry.h"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/assert.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/split_member.hpp>
#ifdef CGAL_USE_GMPXX
  #include <CGAL/mpq_class.h>
  #include <CGAL/mpz_class.h>
#endif

/*
 * @todo SHOULD BE IN io FOLDER OR ADD MISSING detail NAMESPACE
 */

namespace SFCGAL {

namespace io {

/**
 * @brief Binary serialization wrapper for SFCGAL geometries
 *
 * Extends boost::archive::binary_oarchive to provide binary serialization
 * capabilities for SFCGAL geometry objects.
 */
class SFCGAL_API BinarySerializer : public boost::archive::binary_oarchive {
public:
  /**
   * @brief Construct binary serializer
   * @param ostr Output stream to write binary data to
   */
  BinarySerializer(std::ostream &ostr);
};

/**
 * @brief Binary deserialization wrapper for SFCGAL geometries
 *
 * Extends boost::archive::binary_iarchive to provide binary deserialization
 * capabilities for SFCGAL geometry objects.
 */
class SFCGAL_API BinaryUnserializer : public boost::archive::binary_iarchive {
public:
  /**
   * @brief Construct binary deserializer
   * @param istr Input stream to read binary data from
   */
  BinaryUnserializer(std::istream &istr);
};

/**
 * Convert a Geometry to its binary representation
 * @param geom Geometry to convert
 * @return The binary representation as a string
 * @warning resulting string may contain 0s
 */
SFCGAL_API auto
writeBinaryGeometry(const SFCGAL::Geometry &geom) -> std::string;

/**
 * Convert a PreparedGeometry to its binary representation
 * @param preparedGeom PreparedGeometry to convert
 * @return The binary representation as a string
 * @warning resulting string may contain 0s
 */
SFCGAL_API auto
writeBinaryPrepared(const SFCGAL::PreparedGeometry &preparedGeom)
    -> std::string;

/**
 * Read a Geometry from a binary representation
 * @param str binary representation
 * @return A SFCGAL::Geometry
 */
SFCGAL_API auto
readBinaryGeometry(const std::string &str) -> std::unique_ptr<SFCGAL::Geometry>;

/**
 * Read a PreparedGeometry from a binary representation
 * @param str binary representation
 * @return A SFCGAL::PreparedGeometry
 */
SFCGAL_API auto
readBinaryPrepared(const std::string &str)
    -> std::unique_ptr<SFCGAL::PreparedGeometry>;
} // namespace io
} // namespace SFCGAL

namespace boost {
namespace serialization {

/**
 * Serialization of Gmpz for text archives
 */
SFCGAL_API void
save(boost::archive::text_oarchive &ar, const CGAL::Gmpz &z,
     const unsigned int version);

/**
 * Serialization of Gmpz for binary archives
 */
SFCGAL_API void
save(boost::archive::binary_oarchive &ar, const CGAL::Gmpz &z,
     const unsigned int version);

/**
 * Unserialization of Gmpz for text archives
 */
SFCGAL_API void
load(boost::archive::text_iarchive &ar, CGAL::Gmpz &z,
     const unsigned int version);

/**
 * Unserialization of Gmpz for binary archives
 */
SFCGAL_API void
load(boost::archive::binary_iarchive &ar, CGAL::Gmpz &z,
     const unsigned int version);

template <class Archive>
void
serialize(Archive &ar, CGAL::Gmpz &z, const unsigned int version)
{
  split_free(ar, z, version);
}

/**
 * Serializer of Gmpq
 */
template <class Archive>
void
save(Archive &ar, const CGAL::Gmpq &q, const unsigned int /*version*/)
{
  CGAL::Gmpz n = q.numerator();
  CGAL::Gmpz d = q.denominator();
  ar & n;
  ar & d;
}

/**
 * Unserializer of Gmpq
 */
template <class Archive>
void
load(Archive &ar, CGAL::Gmpq &q, const unsigned int /*version*/)
{
  CGAL::Gmpz n;
  CGAL::Gmpz d;
  ar & n;
  ar & d;
  q = CGAL::Gmpq(n, d);
}
template <class Archive>
void
serialize(Archive &ar, CGAL::Gmpq &q, const unsigned int version)
{
  split_free(ar, q, version);
}

#ifdef CGAL_USE_GMPXX
/**
 * Serialization of mpz_class for text archives
 */
SFCGAL_API void
save(boost::archive::text_oarchive &ar, const mpz_class &z,
     const unsigned int version);

/**
 * Serialization of mpz_class for binary archives
 */
SFCGAL_API void
save(boost::archive::binary_oarchive &ar, const mpz_class &z,
     const unsigned int version);

/**
 * Unserialization of mpz_class for text archives
 */
SFCGAL_API void
load(boost::archive::text_iarchive &ar, mpz_class &z,
     const unsigned int version);

/**
 * Unserialization of mpz_class for binary archives
 */
SFCGAL_API void
load(boost::archive::binary_iarchive &ar, mpz_class &z,
     const unsigned int version);

template <class Archive>
void
serialize(Archive &ar, mpz_class &z, const unsigned int version)
{
  split_free(ar, z, version);
}

/**
 * Serializer of mpq_class
 */
template <class Archive>
void
save(Archive &ar, const mpq_class &q, const unsigned int /*version*/)
{
  mpz_class n = q.get_num();
  mpz_class d = q.get_den();
  ar & n;
  ar & d;
}

/**
 * Unserializer of mpq_class
 */
template <class Archive>
void
load(Archive &ar, mpq_class &q, const unsigned int /*version*/)
{
  mpz_class n;
  mpz_class d;
  ar & n;
  ar & d;
  q = mpq_class(n, d);
}
template <class Archive>
void
serialize(Archive &ar, mpq_class &q, const unsigned int version)
{
  split_free(ar, q, version);
}
#endif

/**
 * Serializer of Kernel::FT
 */
template <class Archive>
void
save(Archive &ar, const SFCGAL::Kernel::FT &q, const unsigned int /*version*/)
{
  SFCGAL::Kernel::Exact_kernel::FT eq = CGAL::exact(q);
  ar << eq;
}

/**
 * Unserializer of Kernel::FT
 */
template <class Archive>
void
load(Archive &ar, SFCGAL::Kernel::FT &q, const unsigned int /*version*/)
{
  SFCGAL::Kernel::Exact_kernel::FT eq;
  ar >> eq;
  q = eq;
}
template <class Archive>
void
serialize(Archive &ar, SFCGAL::Kernel::FT &q, const unsigned int version)
{
  split_free(ar, q, version);
}

} // namespace serialization
} // namespace boost

#endif
