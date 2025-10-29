// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TOOLS_CHARARRAYBUFFER_H_
#define SFCGAL_TOOLS_CHARARRAYBUFFER_H_

#include "SFCGAL/config.h"

#include <streambuf>

namespace SFCGAL::tools {
/// Streambuf on a char*
/// http://www.mr-edd.co.uk/blog/beginners_guide_streambuf
class SFCGAL_API CharArrayBuffer : public std::streambuf {
public:
  /**
   * @brief Constructor with begin and end pointers
   * @param begin Pointer to start of character array
   * @param end Pointer to end of character array
   */
  CharArrayBuffer(const char *begin, const char *end);
  /**
   * @brief Constructor with null-terminated string
   * @param str Null-terminated character string
   */
  explicit CharArrayBuffer(const char *str);

private:
  auto
  underflow() -> int_type override;
  auto
  uflow() -> int_type override;
  auto
  pbackfail(int_type ch) -> int_type override;
  auto
  showmanyc() -> std::streamsize override;

  // copy ctor and assignment not implemented;
  // copying not allowed
  CharArrayBuffer(const CharArrayBuffer &) = delete;
  auto
  operator=(const CharArrayBuffer &) -> CharArrayBuffer & = delete;

  auto
  seekpos(std::streampos pos, std::ios_base::openmode)
      -> std::streampos override;
  auto
  seekoff(std::streamoff off, std::ios_base::seekdir way,
          std::ios_base::openmode) -> std::streampos override;

private:
  const char *const begin_;
  const char *const end_;
  const char       *current_;
};

} // namespace SFCGAL::tools

#endif
