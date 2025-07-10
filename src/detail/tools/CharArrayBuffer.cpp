// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/tools/CharArrayBuffer.h"

namespace SFCGAL {
namespace tools {

CharArrayBuffer::CharArrayBuffer(const char *begin, const char *end)
    : begin_(begin), end_(end), current_(begin_)
{
}

CharArrayBuffer::CharArrayBuffer(const char *str)
    : begin_(str), end_(begin_ + std::string(str).length()), current_(begin_)
{
}

auto
CharArrayBuffer::seekoff(std::streamoff off, std::ios_base::seekdir way,
                         std::ios_base::openmode /*__which*/) -> std::streampos
{
  if (way == std::ios_base::cur) {
    if (current_ + off < end_) {
      current_ += off;
    } else {
      return -1;
    }
  } else if (way == std::ios_base::beg) {
    if (begin_ + off < end_) {
      current_ = begin_ + off;
    } else {
      return -1;
    }
  }

  return current_ - begin_;
}

auto
CharArrayBuffer::seekpos(std::streampos pos,
                         std::ios_base::openmode /*__which*/) -> std::streampos
{
  if (begin_ + pos >= end_) {
    return -1;
  }

  current_ = begin_ + pos;
  return current_ - begin_;
}

auto
CharArrayBuffer::underflow() -> CharArrayBuffer::int_type
{
  if (current_ == end_) {
    return traits_type::eof();
  }

  return traits_type::to_int_type(*current_);
}

auto
CharArrayBuffer::uflow() -> CharArrayBuffer::int_type
{
  if (current_ == end_) {
    return traits_type::eof();
  }

  return traits_type::to_int_type(*current_++);
}

auto
CharArrayBuffer::pbackfail(int_type ch) -> CharArrayBuffer::int_type
{
  if (current_ == begin_ || (ch != traits_type::eof() && ch != current_[-1])) {
    return traits_type::eof();
  }

  return traits_type::to_int_type(*--current_);
}

auto
CharArrayBuffer::showmanyc() -> std::streamsize
{
  return end_ - current_;
}

} // namespace tools
} // namespace SFCGAL
