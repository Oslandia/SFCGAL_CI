// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TOOLS_INPUTSTREAMREADER_H_
#define SFCGAL_TOOLS_INPUTSTREAMREADER_H_

#include "SFCGAL/config.h"

#include <cctype>
#include <sstream>
#include <stack>
#include <string>

namespace SFCGAL::tools {

template <typename CharType>
class BasicInputStreamReader;

/**
 * typedef for std::istream
 */
using InputStreamReader = BasicInputStreamReader<char>;
/**
 * typedef for std::wistream
 */
using WInputStreamReader = BasicInputStreamReader<wchar_t>;

/**
 * Helper class to parse data from stream.
 */
template <typename CharType>
class BasicInputStreamReader {
public:
  using char_type = CharType; ///< Character type for the reader
  using string_type =
      typename std::basic_string<char_type>; ///< String type for the reader
  using istream_type =
      typename std::basic_istream<char_type>; ///< Input stream type
  using pos_type =
      typename std::basic_istream<char_type>::pos_type; ///< Stream position
                                                        ///< type

  /**
   * @brief Constructor with an input stream
   * @param s The input stream to read from
   * @param skipWhiteSpaces Whether to skip white spaces during parsing
   */
  BasicInputStreamReader(istream_type &s, bool skipWhiteSpaces = true)
      : _s(s), _states(), _skipWhiteSpaces(skipWhiteSpaces)
  {
    _s >> std::noskipws;
  }

  /**
   * @brief Try to match a char in the input stream
   * @param c The character to match
   * @return True if character matches
   */
  auto
  match(char_type const &c) -> bool
  {
    begin();

    if (_skipWhiteSpaces) {
      skipWhiteSpaces();
    }

    if (!_s.eof() && _s.get() == c) {
      commit();
      return true;
    } else {
      rollback();
      return false;
    }
  }

  /**
   * @brief Try to match a char in the input stream, case-insensitive variant
   * @param c The character to match
   * @return True if character matches (case-insensitive)
   */
  auto
  imatch(char_type const &c) -> bool
  {
    begin();

    if (_skipWhiteSpaces) {
      skipWhiteSpaces();
    }

    if (!_s.eof() && ::tolower(_s.get()) == ::tolower(c)) {
      commit();
      return true;
    } else {
      rollback();
      return false;
    }
  }

  /**
   * @brief Try to match a string in the input stream
   * @param str The string to match
   * @return True if string matches
   */
  auto
  match(string_type const &str) -> bool
  {
    begin();

    if (_skipWhiteSpaces) {
      skipWhiteSpaces();
    }

    for (typename string_type::const_iterator it = str.begin(); it != str.end();
         ++it) {
      if (!_s.eof() && _s.get() == *it) {
        continue;
      }

      rollback();
      return false;
    }

    commit();
    return true;
  }

  /**
   * @brief Try to match a string in the input stream, case-insensitive variant
   * @param str The string to match
   * @return True if string matches (case-insensitive)
   */
  auto
  imatch(string_type const &str) -> bool
  {
    begin();

    if (_skipWhiteSpaces) {
      skipWhiteSpaces();
    }

    for (typename string_type::const_iterator it = str.begin(); it != str.end();
         ++it) {
      if (!_s.eof() && ::tolower(_s.get()) == ::tolower(*it)) {
        continue;
      }

      rollback();
      return false;
    }

    commit();
    return true;
  }

  /**
   * @brief Try to read a value in the input stream, restore state if read fails
   * @param value Reference to store the read value
   * @return True if read was successful
   */
  template <typename T>
  bool
  read(T &value)
  {
    begin();

    if (_skipWhiteSpaces) {
      skipWhiteSpaces();
    }

    if (_s >> CGAL::iformat(value)) {
      commit();
      return true;
    } else {
      rollback();
      return false;
    }
  }

  /**
   * @brief Read bytes from input stream into buffer
   * @param buffer String buffer to store read bytes
   * @param bytesToRead Number of bytes to read
   */
  auto
  readBytes(std::string &buffer, size_t bytesToRead) -> void
  {
    begin();

    _s.read(&buffer[0], bytesToRead);

    commit();
  }
  /// \brief save input stream state (read position)
  void
  begin()
  {
    _states.push(_s.tellg());
  }
  /// \brief validate read from input stream
  void
  commit()
  {
    _states.pop();
  }
  /// \brief restore stream state (read position)
  void
  rollback()
  {
    _s.seekg(_states.top());
    _s.clear();
    _states.pop();
  }

  /**
   * @brief Test if read is complete (either tried to read after eof, either on
   * eof)
   * @return True if at end of stream
   */
  bool
  eof() const
  {
    return _s.eof() || (_s.peek() == std::char_traits<char_type>::eof());
  }

  /**
   * @brief Returns the wrapped stream
   * @return Reference to the input stream
   */
  inline istream_type &
  s()
  {
    return _s;
  }
  /**
   * @brief Returns the wrapped stream
   * @return Const reference to the input stream
   */
  inline istream_type const &
  s() const
  {
    return _s;
  }

  /**
   * @brief Returns a string corresponding to the current state
   * @param nMax Maximum number of characters to include in context
   * @return String showing current stream context
   */
  string_type
  context(size_t nMax = 20)
  {
    begin();
    std::basic_ostringstream<char_type> oss;

    for (size_t i = 0; i < nMax; i++) {
      if (eof()) {
        break;
      }

      oss << (char_type)s().get();
    }

    if (!eof()) {
      oss << "...";
    }

    rollback();
    return oss.str();
  }

private:
  /// \brief the input stream
  istream_type &_s;
  /// \brief read position saved
  std::stack<pos_type> _states;
  /// \brief indicates if white chars should be skipped
  bool _skipWhiteSpaces;

  /// \brief skip white spaces
  void
  skipWhiteSpaces()
  {
    while (!_s.eof() && std::isspace(_s.peek())) {
      _s.get();
      continue;
    }
  }

  /// \brief no copy
  BasicInputStreamReader(BasicInputStreamReader const &other);
};

} // namespace SFCGAL::tools

#endif
