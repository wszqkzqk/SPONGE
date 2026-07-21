#pragma once

#include <map>
#include <set>
#include <string>

namespace sponge::toml_wrap
{
bool ParseAndFlatten(const std::string& content, const std::string& source_path,
                     std::map<std::string, std::string>* parsed_commands,
                     std::set<std::string>* scalar_string_keys,
                     std::string* error_message);
}
