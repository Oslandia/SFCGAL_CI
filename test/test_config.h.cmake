#ifndef _TEST_CONFIG_H_
#define _TEST_CONFIG_H_

#include <algorithm>
#include <string>
#include <algorithm>

#include <nlohmann/json.hpp>

#cmakedefine SFCGAL_TEST_DIRECTORY "@SFCGAL_TEST_DIRECTORY@"

inline std::string random_string( size_t length = 12 )
{
    auto randchar = []() -> char
    {
        const char charset[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
        const size_t max_index = (sizeof(charset) - 1);
        return charset[ rand() % max_index ];
    };
    std::string str(length,0);
    std::generate_n( str.begin(), length, randchar );
    return str;
}

// Compare JSON objects after sorting to make the test order-independent
inline bool compare_json(const std::string &actualStr, const std::string &expectedStr) {
  nlohmann::json actualJson = nlohmann::json::parse(actualStr);
  nlohmann::json expectedJson = nlohmann::json::parse(expectedStr);

  auto normalize_array = [](nlohmann::json& jsonData)
  {
    if (!jsonData.is_array())
    {
      return;
    }

    std::sort(jsonData.begin(), jsonData.end(), [](const nlohmann::json& a, const nlohmann::json& b)
    {
      if (a.is_object() && b.is_object() && a.contains("name") && b.contains("name"))
      {
        return a.at("name") < b.at("name");
      }
      return false;
    });
  };

  normalize_array(actualJson);
  normalize_array(expectedJson);

  return actualJson == expectedJson;
}

#endif
