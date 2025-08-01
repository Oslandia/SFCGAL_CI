// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TOOLS_CHARARRAYBUFFER_H_
#define SFCGAL_TOOLS_CHARARRAYBUFFER_H_

#include "SFCGAL/config.h"

#include <streambuf>

namespace SFCGAL {
namespace tools {
/// Streambuf on a char*
/// http://www.mr-edd.co.uk/blog/beginners_guide_streambuf
class SFCGAL_API CharArrayBuffer : public std::streambuf {
public:
  CharArrayBuffer(const char *begin, const char *end);
  explicit CharArrayBuffer(const char *str);

private:
  int_type
  underflow() override;
  int_type
  uflow() override;
  int_type
  pbackfail(int_type ch) override;
  std::streamsize
  showmanyc() override;

  // copy ctor and assignment not implemented;
  // copying not allowed
  CharArrayBuffer(const CharArrayBuffer &);
  CharArrayBuffer &
  operator=(const CharArrayBuffer &);

  std::streampos
  seekpos(std::streampos pos, std::ios_base::openmode) override;
  std::streampos
  seekoff(std::streamoff off, std::ios_base::seekdir way,
          std::ios_base::openmode) override;

private:
  const char *const begin_;
  const char *const end_;
  const char       *current_;
};

} // namespace tools
} // namespace SFCGAL

#endif
