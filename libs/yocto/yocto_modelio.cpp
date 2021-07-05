//
// Implementation for Yocto/Ply.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

#include "yocto_modelio.h"

#include <cstdio>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "yocto_color.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::string_view;
using std::unordered_map;
using std::unordered_set;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Formats values to string
inline void format_value(string& str, const string& value) { str += value; }
inline void format_value(string& str, int8_t value) {
  str += std::to_string((int32_t)value);
}
inline void format_value(string& str, int16_t value) {
  str += std::to_string(value);
}
inline void format_value(string& str, int32_t value) {
  str += std::to_string(value);
}
inline void format_value(string& str, int64_t value) {
  str += std::to_string(value);
}
inline void format_value(string& str, uint8_t value) {
  str += std::to_string((uint32_t)value);
}
inline void format_value(string& str, uint16_t value) {
  str += std::to_string(value);
}
inline void format_value(string& str, uint32_t value) {
  str += std::to_string(value);
}
inline void format_value(string& str, uint64_t value) {
  str += std::to_string(value);
}
inline void format_value(string& str, float value) {
  auto buf = array<char, 256>{};
  snprintf(buf.data(), buf.size(), "%g", value);
  str += buf.data();
}
inline void format_value(string& str, double value) {
  auto buf = array<char, 256>{};
  snprintf(buf.data(), buf.size(), "%g", value);
  str += buf.data();
}
inline void format_value(string& str, const vec2f& value) {
  for (auto i = 0; i < 2; i++) {
    if (i != 0) str += " ";
    format_value(str, value[i]);
  }
}
inline void format_value(string& str, const vec3f& value) {
  for (auto i = 0; i < 3; i++) {
    if (i != 0) str += " ";
    format_value(str, value[i]);
  }
}
inline void format_value(string& str, const frame3f& value) {
  for (auto i = 0; i < 4; i++) {
    if (i != 0) str += " ";
    format_value(str, value[i]);
  }
}
inline void format_value(string& str, const vec4f& value) {
  for (auto i = 0; i < 4; i++) {
    if (i != 0) str += " ";
    format_value(str, value[i]);
  }
}
inline void format_value(string& str, const mat4f& value) {
  for (auto i = 0; i < 4; i++) {
    if (i != 0) str += " ";
    format_value(str, value[i]);
  }
}

// Foramt to file
inline void format_values(string& str, const string& fmt) {
  auto pos = fmt.find("{}");
  if (pos != string::npos) throw std::invalid_argument("bad format string");
  str += fmt;
}
template <typename Arg, typename... Args>
inline void format_values(
    string& str, const string& fmt, const Arg& arg, const Args&... args) {
  auto pos = fmt.find("{}");
  if (pos == string::npos) throw std::invalid_argument("bad format string");
  str += fmt.substr(0, pos);
  format_value(str, arg);
  format_values(str, fmt.substr(pos + 2), args...);
}

template <typename... Args>
inline void format_values(
    file_stream& fs, const string& fmt, const Args&... args) {
  auto str = ""s;
  format_values(str, fmt, args...);
  write_text(fs, str);
}
template <typename T>
inline void format_value(file_stream& fs, const T& value) {
  auto str = ""s;
  format_value(str, value);
  write_text(fs, str);
}

inline bool is_newline(char c) { return c == '\r' || c == '\n'; }
inline bool is_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}
inline void skip_whitespace(string_view& str) {
  while (!str.empty() && is_space(str.front())) str.remove_prefix(1);
}

inline void remove_comment(
    string_view& str, char comment_char = '#', bool handle_quotes = false) {
  if (!handle_quotes) {
    while (!str.empty() && is_newline(str.back())) str.remove_suffix(1);
    auto cpy = str;
    while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
    str.remove_suffix(cpy.size());
  } else {
    while (!str.empty() && is_newline(str.back())) str.remove_suffix(1);
    auto cpy       = str;
    auto in_string = false;
    while (!cpy.empty()) {
      if (cpy.front() == '"') in_string = !in_string;
      if (cpy.front() == comment_char && !in_string) break;
      cpy.remove_prefix(1);
    }
    str.remove_suffix(cpy.size());
  }
}

// Parse values from a string
inline void parse_value(string_view& str, string_view& value) {
  skip_whitespace(str);
  if (str.empty()) throw std::invalid_argument{"string expected"};
  if (str.front() != '"') {
    auto cpy = str;
    while (!cpy.empty() && !is_space(cpy.front())) cpy.remove_prefix(1);
    value = str;
    value.remove_suffix(cpy.size());
    str.remove_prefix(str.size() - cpy.size());
  } else {
    if (str.front() != '"') throw std::invalid_argument{"string expected"};
    str.remove_prefix(1);
    if (str.empty()) throw std::invalid_argument{"string expected"};
    auto cpy = str;
    while (!cpy.empty() && cpy.front() != '"') cpy.remove_prefix(1);
    if (cpy.empty()) throw std::invalid_argument{"string expected"};
    value = str;
    value.remove_suffix(cpy.size());
    str.remove_prefix(str.size() - cpy.size());
    str.remove_prefix(1);
  }
}
inline void parse_value(string_view& str, string& value) {
  auto valuev = string_view{};
  parse_value(str, valuev);
  value = string{valuev};
}
inline void parse_value(string_view& str, int8_t& value) {
  char* end = nullptr;
  value     = (int8_t)strtol(str.data(), &end, 10);
  if (str.data() == end) throw std::invalid_argument{"integer expected"};
  str.remove_prefix(end - str.data());
}
inline void parse_value(string_view& str, int16_t& value) {
  char* end = nullptr;
  value     = (int16_t)strtol(str.data(), &end, 10);
  if (str.data() == end) throw std::invalid_argument{"integer expected"};
  str.remove_prefix(end - str.data());
}
inline void parse_value(string_view& str, int32_t& value) {
  char* end = nullptr;
  value     = (int32_t)strtol(str.data(), &end, 10);
  if (str.data() == end) throw std::invalid_argument{"integer expected"};
  str.remove_prefix(end - str.data());
}
inline void parse_value(string_view& str, int64_t& value) {
  char* end = nullptr;
  value     = (int64_t)strtoll(str.data(), &end, 10);
  if (str.data() == end) throw std::invalid_argument{"integer expected"};
  str.remove_prefix(end - str.data());
}
inline void parse_value(string_view& str, uint8_t& value) {
  char* end = nullptr;
  value     = (uint8_t)strtoul(str.data(), &end, 10);
  if (str.data() == end) throw std::invalid_argument{"integer expected"};
  str.remove_prefix(end - str.data());
}
inline void parse_value(string_view& str, uint16_t& value) {
  char* end = nullptr;
  value     = (uint16_t)strtoul(str.data(), &end, 10);
  if (str.data() == end) throw std::invalid_argument{"integer expected"};
  str.remove_prefix(end - str.data());
}
inline void parse_value(string_view& str, uint32_t& value) {
  char* end = nullptr;
  value     = (uint32_t)strtoul(str.data(), &end, 10);
  if (str.data() == end) throw std::invalid_argument{"integer expected"};
  str.remove_prefix(end - str.data());
}
inline void parse_value(string_view& str, uint64_t& value) {
  char* end = nullptr;
  value     = (uint64_t)strtoull(str.data(), &end, 10);
  if (str.data() == end) throw std::invalid_argument{"integer expected"};
  str.remove_prefix(end - str.data());
}
inline void parse_value(string_view& str, float& value) {
  char* end = nullptr;
  value     = strtof(str.data(), &end);
  if (str.data() == end) throw std::invalid_argument{"number expected"};
  str.remove_prefix(end - str.data());
}
inline void parse_value(string_view& str, double& value) {
  char* end = nullptr;
  value     = strtod(str.data(), &end);
  if (str.data() == end) throw std::invalid_argument{"number expected"};
  str.remove_prefix(end - str.data());
}
#ifdef __APPLE__
inline void parse_value(string_view& str, size_t& value) {
  char* end = nullptr;
  value     = (size_t)strtoull(str.data(), &end, 10);
  if (str.data() == end) throw std::invalid_argument{"integer expected"};
  str.remove_prefix(end - str.data());
}
#endif
inline void parse_value(string_view& str, bool& value) {
  auto valuei = 0;
  parse_value(str, valuei);
  value = (bool)valuei;
}

inline void parse_value(string_view& str, vec2f& value) {
  for (auto i = 0; i < 2; i++) parse_value(str, value[i]);
}
inline void parse_value(string_view& str, vec3f& value) {
  for (auto i = 0; i < 3; i++) parse_value(str, value[i]);
}
inline void parse_value(string_view& str, vec4f& value) {
  for (auto i = 0; i < 4; i++) parse_value(str, value[i]);
}
inline void parse_value(string_view& str, mat3f& value) {
  for (auto i = 0; i < 3; i++) parse_value(str, value[i]);
}
inline void parse_value(string_view& str, mat4f& value) {
  for (auto i = 0; i < 4; i++) parse_value(str, value[i]);
}
inline void parse_value(string_view& str, frame3f& value) {
  for (auto i = 0; i < 4; i++) parse_value(str, value[i]);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PLY LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto {

// Load and save ply
ply_model load_ply(const string& filename) {
  auto ply = ply_model{};
  load_ply(filename, ply);
  return ply;
}

// Load ply
void load_ply(const string& filename, ply_model& ply) {
  // ply type names
  static auto type_map = unordered_map<string, ply_type>{{"char", ply_type::i8},
      {"short", ply_type::i16}, {"int", ply_type::i32}, {"long", ply_type::i64},
      {"uchar", ply_type::u8}, {"ushort", ply_type::u16},
      {"uint", ply_type::u32}, {"ulong", ply_type::u64},
      {"float", ply_type::f32}, {"double", ply_type::f64},
      {"int8", ply_type::i8}, {"int16", ply_type::i16},
      {"int32", ply_type::i32}, {"int64", ply_type::i64},
      {"uint8", ply_type::u8}, {"uint16", ply_type::u16},
      {"uint32", ply_type::u32}, {"uint64", ply_type::u64},
      {"float32", ply_type::f32}, {"float64", ply_type::f64}};

  // open file
  auto fs = open_file(filename, "rb");

  // parsing checks
  auto first_line = true;
  auto end_header = false;

  // read header ---------------------------------------------
  auto buffer = array<char, 4096>{};
  while (read_line(fs, buffer)) {
    try {
      // str
      auto str = string_view{buffer.data()};
      remove_comment(str);
      skip_whitespace(str);
      if (str.empty()) continue;

      // get command
      auto cmd = ""s;
      parse_value(str, cmd);
      if (cmd.empty()) continue;

      // check magic number
      if (first_line) {
        if (cmd != "ply") throw io_error::parse_error(filename);
        first_line = false;
        continue;
      }

      // possible token values
      if (cmd == "ply") {
        if (!first_line) throw io_error::parse_error(filename);
      } else if (cmd == "format") {
        auto fmt = ""s;
        parse_value(str, fmt);
        if (fmt == "ascii") {
          ply.format = ply_format::ascii;
        } else if (fmt == "binary_little_endian") {
          ply.format = ply_format::binary_little_endian;
        } else if (fmt == "binary_big_endian") {
          ply.format = ply_format::binary_big_endian;
        } else {
          throw io_error::parse_error(filename);
        }
      } else if (cmd == "comment") {
        skip_whitespace(str);
        ply.comments.emplace_back(str);
      } else if (cmd == "obj_info") {
        skip_whitespace(str);
        // comment is the rest of the str
      } else if (cmd == "element") {
        auto& elem = ply.elements.emplace_back();
        parse_value(str, elem.name);
        parse_value(str, elem.count);
      } else if (cmd == "property") {
        if (ply.elements.empty()) throw io_error::parse_error(filename);
        auto& prop  = ply.elements.back().properties.emplace_back();
        auto  tname = ""s;
        parse_value(str, tname);
        if (tname == "list") {
          prop.is_list = true;
          parse_value(str, tname);
          auto itype = type_map.at(tname);
          if (itype != ply_type::u8) throw io_error::parse_error(filename);
          parse_value(str, tname);
          prop.type = type_map.at(tname);
        } else {
          prop.is_list = false;
          prop.type    = type_map.at(tname);
        }
        parse_value(str, prop.name);
      } else if (cmd == "end_header") {
        end_header = true;
        break;
      } else {
        throw io_error::parse_error(filename);
      }
    } catch (...) {
      throw io_error::parse_error(filename);
    }
  }

  // check exit
  if (!end_header) throw io_error::parse_error(filename);

  // allocate data ---------------------------------
  for (auto& element : ply.elements) {
    for (auto& property : element.properties) {
      auto count = property.is_list ? element.count * 3 : element.count;
      switch (property.type) {
        case ply_type::i8: property.data_i8.reserve(count); break;
        case ply_type::i16: property.data_i16.reserve(count); break;
        case ply_type::i32: property.data_i32.reserve(count); break;
        case ply_type::i64: property.data_i64.reserve(count); break;
        case ply_type::u8: property.data_u8.reserve(count); break;
        case ply_type::u16: property.data_u16.reserve(count); break;
        case ply_type::u32: property.data_u32.reserve(count); break;
        case ply_type::u64: property.data_u64.reserve(count); break;
        case ply_type::f32: property.data_f32.reserve(count); break;
        case ply_type::f64: property.data_f64.reserve(count); break;
      }
      if (property.is_list) property.ldata_u8.reserve(element.count);
    }
  }

  // read data -------------------------------------
  if (ply.format == ply_format::ascii) {
    auto buffer = array<char, 4096>{};
    for (auto& elem : ply.elements) {
      for (auto idx = 0; idx < elem.count; idx++) {
        if (!read_line(fs, buffer)) throw io_error::read_error(filename);
        try {
          auto str = string_view{buffer.data()};
          for (auto& prop : elem.properties) {
            if (prop.is_list) parse_value(str, prop.ldata_u8.emplace_back());
            auto vcount = prop.is_list ? prop.ldata_u8.back() : 1;
            for (auto i = 0; i < vcount; i++) {
              switch (prop.type) {
                case ply_type::i8:
                  parse_value(str, prop.data_i8.emplace_back());
                  break;
                case ply_type::i16:
                  parse_value(str, prop.data_i16.emplace_back());
                  break;
                case ply_type::i32:
                  parse_value(str, prop.data_i32.emplace_back());
                  break;
                case ply_type::i64:
                  parse_value(str, prop.data_i64.emplace_back());
                  break;
                case ply_type::u8:
                  parse_value(str, prop.data_u8.emplace_back());
                  break;
                case ply_type::u16:
                  parse_value(str, prop.data_u16.emplace_back());
                  break;
                case ply_type::u32:
                  parse_value(str, prop.data_u32.emplace_back());
                  break;
                case ply_type::u64:
                  parse_value(str, prop.data_u64.emplace_back());
                  break;
                case ply_type::f32:
                  parse_value(str, prop.data_f32.emplace_back());
                  break;
                case ply_type::f64:
                  parse_value(str, prop.data_f64.emplace_back());
                  break;
              }
            }
          }
        } catch (...) {
          throw io_error::parse_error(filename);
        }
      }
    }
  } else {
    auto big_endian = ply.format == ply_format::binary_big_endian;
    for (auto& elem : ply.elements) {
      for (auto idx = 0; idx < elem.count; idx++) {
        for (auto& prop : elem.properties) {
          if (prop.is_list) {
            read_value(fs, prop.ldata_u8.emplace_back(), big_endian);
          }
          auto vcount = prop.is_list ? prop.ldata_u8.back() : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop.type) {
              case ply_type::i8:
                read_value(fs, prop.data_i8.emplace_back(), big_endian);
                break;
              case ply_type::i16:
                read_value(fs, prop.data_i16.emplace_back(), big_endian);
                break;
              case ply_type::i32:
                read_value(fs, prop.data_i32.emplace_back(), big_endian);
                break;
              case ply_type::i64:
                read_value(fs, prop.data_i64.emplace_back(), big_endian);
                break;
              case ply_type::u8:
                read_value(fs, prop.data_u8.emplace_back(), big_endian);
                break;
              case ply_type::u16:
                read_value(fs, prop.data_u16.emplace_back(), big_endian);
                break;
              case ply_type::u32:
                read_value(fs, prop.data_u32.emplace_back(), big_endian);
                break;
              case ply_type::u64:
                read_value(fs, prop.data_u64.emplace_back(), big_endian);
                break;
              case ply_type::f32:
                read_value(fs, prop.data_f32.emplace_back(), big_endian);
                break;
              case ply_type::f64:
                read_value(fs, prop.data_f64.emplace_back(), big_endian);
                break;
            }
          }
        }
      }
    }
  }
}

// save ply
void save_ply(const string& filename, const ply_model& ply) {
  // ply type names
  static auto type_map = unordered_map<ply_type, string>{{ply_type::i8, "char"},
      {ply_type::i16, "short"}, {ply_type::i32, "int"}, {ply_type::i64, "uint"},
      {ply_type::u8, "uchar"}, {ply_type::u16, "ushort"},
      {ply_type::u32, "uint"}, {ply_type::u64, "ulong"},
      {ply_type::f32, "float"}, {ply_type::f64, "double"}};
  static auto format_map = unordered_map<ply_format, string>{
      {ply_format::ascii, "ascii"},
      {ply_format::binary_little_endian, "binary_little_endian"},
      {ply_format::binary_big_endian, "binary_big_endian"}};

  // open file
  auto fs = open_file(filename, "wb");

  // header
  format_values(fs, "ply\n");
  format_values(fs, "format {} 1.0\n", format_map.at(ply.format));
  format_values(fs, "comment Written by Yocto/GL\n");
  format_values(fs, "comment https://github.com/xelatihy/yocto-gl\n");
  for (auto& comment : ply.comments) format_values(fs, "comment {}\n", comment);
  for (auto& elem : ply.elements) {
    format_values(fs, "element {} {}\n", elem.name, (uint64_t)elem.count);
    for (auto& prop : elem.properties) {
      if (prop.is_list) {
        format_values(
            fs, "property list uchar {} {}\n", type_map[prop.type], prop.name);
      } else {
        format_values(fs, "property {} {}\n", type_map[prop.type], prop.name);
      }
    }
  }

  format_values(fs, "end_header\n");

  // properties
  if (ply.format == ply_format::ascii) {
    for (auto& elem : ply.elements) {
      auto cur = vector<size_t>(elem.properties.size(), 0);
      for (auto idx = 0; idx < elem.count; idx++) {
        for (auto& prop : elem.properties) {
          if (prop.is_list) format_values(fs, "{} ", (int)prop.ldata_u8[idx]);
          auto vcount = prop.is_list ? prop.ldata_u8[idx] : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop.type) {
              case ply_type::i8:
                format_values(fs, "{} ", prop.data_i8[cur[idx]++]);
                break;
              case ply_type::i16:
                format_values(fs, "{} ", prop.data_i16[cur[idx]++]);
                break;
              case ply_type::i32:
                format_values(fs, "{} ", prop.data_i32[cur[idx]++]);
                break;
              case ply_type::i64:
                format_values(fs, "{} ", prop.data_i64[cur[idx]++]);
                break;
              case ply_type::u8:
                format_values(fs, "{} ", prop.data_u8[cur[idx]++]);
                break;
              case ply_type::u16:
                format_values(fs, "{} ", prop.data_u16[cur[idx]++]);
                break;
              case ply_type::u32:
                format_values(fs, "{} ", prop.data_u32[cur[idx]++]);
                break;
              case ply_type::u64:
                format_values(fs, "{} ", prop.data_u64[cur[idx]++]);
                break;
              case ply_type::f32:
                format_values(fs, "{} ", prop.data_f32[cur[idx]++]);
                break;
              case ply_type::f64:
                format_values(fs, "{} ", prop.data_f64[cur[idx]++]);
                break;
            }
          }
          format_values(fs, "\n");
        }
      }
    }
  } else {
    auto big_endian = ply.format == ply_format::binary_big_endian;
    for (auto& elem : ply.elements) {
      auto cur = vector<size_t>(elem.properties.size(), 0);
      for (auto idx = 0; idx < elem.count; idx++) {
        for (auto pidx = 0; pidx < elem.properties.size(); pidx++) {
          auto& prop = elem.properties[pidx];
          if (prop.is_list) write_value(fs, prop.ldata_u8[idx], big_endian);
          auto vcount = prop.is_list ? prop.ldata_u8[idx] : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop.type) {
              case ply_type::i8:
                write_value(fs, prop.data_i8[cur[pidx]++], big_endian);
                break;
              case ply_type::i16:
                write_value(fs, prop.data_i16[cur[pidx]++], big_endian);
                break;
              case ply_type::i32:
                write_value(fs, prop.data_i32[cur[pidx]++], big_endian);
                break;
              case ply_type::i64:
                write_value(fs, prop.data_i64[cur[pidx]++], big_endian);
                break;
              case ply_type::u8:
                write_value(fs, prop.data_u8[cur[pidx]++], big_endian);
                break;
              case ply_type::u16:
                write_value(fs, prop.data_u16[cur[pidx]++], big_endian);
                break;
              case ply_type::u32:
                write_value(fs, prop.data_u32[cur[pidx]++], big_endian);
                break;
              case ply_type::u64:
                write_value(fs, prop.data_u64[cur[pidx]++], big_endian);
                break;
              case ply_type::f32:
                write_value(fs, prop.data_f32[cur[pidx]++], big_endian);
                break;
              case ply_type::f64:
                write_value(fs, prop.data_f64[cur[pidx]++], big_endian);
                break;
            }
          }
        }
      }
    }
  }
}

// Load ply
bool load_ply(const string& filename, ply_model& ply, string& error) {
  try {
    load_ply(filename, ply);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Save ply
bool save_ply(const string& filename, const ply_model& ply, string& error) {
  try {
    save_ply(filename, ply);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Get ply properties
bool has_property(
    const ply_model& ply, const string& element, const string& property) {
  for (auto& elem : ply.elements) {
    if (elem.name != element) continue;
    for (auto& prop : elem.properties) {
      if (prop.name == property) return true;
    }
  }
  return false;
}
ply_property& get_property(
    ply_model& ply, const string& element, const string& property) {
  for (auto& elem : ply.elements) {
    if (elem.name != element) continue;
    for (auto& prop : elem.properties) {
      if (prop.name == property) return prop;
    }
  }
  throw std::runtime_error("property not found");
}
const ply_property& get_property(
    const ply_model& ply, const string& element, const string& property) {
  for (auto& elem : ply.elements) {
    if (elem.name != element) continue;
    for (auto& prop : elem.properties) {
      if (prop.name == property) return prop;
    }
  }
  throw std::runtime_error("property not found");
}
template <typename T, typename T1>
inline bool convert_property(const vector<T1>& prop, vector<T>& values) {
  values = vector<T>(prop.size());
  for (auto i = (size_t)0; i < prop.size(); i++) values[i] = (T)prop[i];
  return true;
}
template <typename T>
inline bool convert_property(const ply_property& prop, vector<T>& values) {
  switch (prop.type) {
    case ply_type::i8: return convert_property(prop.data_i8, values);
    case ply_type::i16: return convert_property(prop.data_i16, values);
    case ply_type::i32: return convert_property(prop.data_i32, values);
    case ply_type::i64: return convert_property(prop.data_i64, values);
    case ply_type::u8: return convert_property(prop.data_u8, values);
    case ply_type::u16: return convert_property(prop.data_u16, values);
    case ply_type::u32: return convert_property(prop.data_u32, values);
    case ply_type::u64: return convert_property(prop.data_u64, values);
    case ply_type::f32: return convert_property(prop.data_f32, values);
    case ply_type::f64: return convert_property(prop.data_f64, values);
  }
  // return here to silence warnings
  throw std::runtime_error{"should not have gotten here"};
  return false;
}
bool get_value(const ply_model& ply, const string& element,
    const string& property, vector<float>& values) {
  values.clear();
  if (!has_property(ply, element, property)) return false;
  auto& prop = get_property(ply, element, property);
  if (prop.is_list) return false;
  if (!convert_property(prop, values)) return false;
  return true;
}
bool get_values(const ply_model& ply, const string& element,
    const array<string, 2>& properties, vector<vec2f>& values) {
  values.clear();
  auto x = vector<float>{}, y = vector<float>{};
  if (!get_value(ply, element, properties[0], x)) return false;
  if (!get_value(ply, element, properties[1], y)) return false;
  values = vector<vec2f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++) values[i] = {x[i], y[i]};
  return true;
}
bool get_values(const ply_model& ply, const string& element,
    const array<string, 3>& properties, vector<vec3f>& values) {
  values.clear();
  auto x = vector<float>{}, y = vector<float>{}, z = vector<float>{};
  if (!get_value(ply, element, properties[0], x)) return false;
  if (!get_value(ply, element, properties[1], y)) return false;
  if (!get_value(ply, element, properties[2], z)) return false;
  values = vector<vec3f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++)
    values[i] = {x[i], y[i], z[i]};
  return true;
}
bool get_values(const ply_model& ply, const string& element,
    const array<string, 4>& properties, vector<vec4f>& values) {
  values.clear();
  auto x = vector<float>{}, y = vector<float>{}, z = vector<float>{},
       w = vector<float>{};
  if (!get_value(ply, element, properties[0], x)) return false;
  if (!get_value(ply, element, properties[1], y)) return false;
  if (!get_value(ply, element, properties[2], z)) return false;
  if (!get_value(ply, element, properties[3], w)) return false;
  values = vector<vec4f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++)
    values[i] = {x[i], y[i], z[i], w[i]};
  return true;
}
bool get_values(const ply_model& ply, const string& element,
    const array<string, 12>& properties, vector<frame3f>& values) {
  values.clear();
  auto coords = array<vector<float>, 12>{};
  for (auto idx = 0; idx < 12; idx++)
    if (!get_value(ply, element, properties[idx], coords[idx])) return false;
  values = vector<frame3f>(coords[0].size());
  for (auto i = (size_t)0; i < values.size(); i++) {
    for (auto c = 0; c < 12; c++) (&values[i].x.x)[c] = coords[c][i];
  }
  return true;
}
bool get_lists(const ply_model& ply, const string& element,
    const string& property, vector<vector<int>>& lists) {
  lists.clear();
  if (!has_property(ply, element, property)) return false;
  auto& prop = get_property(ply, element, property);
  if (!prop.is_list) return false;
  auto& sizes  = prop.ldata_u8;
  auto  values = vector<int>{};
  if (!convert_property(prop, values)) return false;
  lists    = vector<vector<int>>(sizes.size());
  auto cur = (size_t)0;
  for (auto i = (size_t)0; i < lists.size(); i++) {
    lists[i].resize(sizes[i]);
    for (auto c = 0; c < sizes[i]; c++) {
      lists[i][c] = values[cur++];
    }
  }
  return true;
}
bool get_list_sizes(const ply_model& ply, const string& element,
    const string& property, vector<byte>& sizes) {
  if (!has_property(ply, element, property)) return {};
  auto& prop = get_property(ply, element, property);
  if (!prop.is_list) return {};
  sizes = prop.ldata_u8;
  return true;
}
bool get_list_values(const ply_model& ply, const string& element,
    const string& property, vector<int>& values) {
  if (!has_property(ply, element, property)) return {};
  auto& prop = get_property(ply, element, property);
  if (!prop.is_list) return {};
  return convert_property<int>(prop, values);
}

inline vector<vec2f> flip_ply_texcoord(const vector<vec2f>& texcoords) {
  auto flipped = texcoords;
  for (auto& uv : flipped) uv.y = 1 - uv.y;
  return flipped;
}

// Get ply properties for meshes
bool get_positions(const ply_model& ply, vector<vec3f>& positions) {
  return get_values(ply, "vertex", {"x", "y", "z"}, positions);
}
bool get_normals(const ply_model& ply, vector<vec3f>& normals) {
  return get_values(ply, "vertex", {"nx", "ny", "nz"}, normals);
}
bool get_texcoords(const ply_model& ply, vector<vec2f>& texcoords, bool flipv) {
  if (has_property(ply, "vertex", "u")) {
    if (!get_values(ply, "vertex", {"u", "v"}, texcoords)) return false;
  } else {
    if (!get_values(ply, "vertex", {"s", "t"}, texcoords)) return false;
  }
  if (flipv) {
    for (auto& uv : texcoords) uv.y = 1 - uv.y;
  }
  return true;
}
bool get_colors(const ply_model& ply, vector<vec3f>& colors) {
  return get_values(ply, "vertex", {"red", "green", "blue"}, colors);
}
bool get_colors(const ply_model& ply, vector<vec4f>& colors) {
  if (has_property(ply, "vertex", "alpha")) {
    return get_values(ply, "vertex", {"red", "green", "blue", "alpha"}, colors);
  } else {
    auto colors3 = vector<vec3f>{};
    if (!get_values(ply, "vertex", {"red", "green", "blue"}, colors3))
      return false;
    colors.resize(colors3.size());
    for (auto i = 0; i < colors.size(); i++)
      colors[i] = {colors3[i].x, colors3[i].y, colors3[i].z, 1};
    return true;
  }
}
bool get_radius(const ply_model& ply, vector<float>& radius) {
  return get_value(ply, "vertex", "radius", radius);
}
bool get_faces(const ply_model& ply, vector<vector<int>>& faces) {
  return get_lists(ply, "face", "vertex_indices", faces);
}
bool get_triangles(const ply_model& ply, vector<vec3i>& triangles) {
  triangles.clear();
  auto indices = vector<int>{};
  auto sizes   = vector<uint8_t>{};
  if (!get_list_values(ply, "face", "vertex_indices", indices)) return false;
  if (!get_list_sizes(ply, "face", "vertex_indices", sizes)) return false;
  triangles = vector<vec3i>{};
  triangles.reserve(sizes.size());
  auto cur = 0;
  for (auto size : sizes) {
    for (auto c = 2; c < size; c++) {
      triangles.push_back(
          {indices[cur + 0], indices[cur + c - 1], indices[cur + c]});
    }
    cur += size;
  }
  return true;
}
bool get_quads(const ply_model& ply, vector<vec4i>& quads) {
  quads.clear();
  auto indices = vector<int>{};
  auto sizes   = vector<uint8_t>{};
  if (!get_list_values(ply, "face", "vertex_indices", indices)) return false;
  if (!get_list_sizes(ply, "face", "vertex_indices", sizes)) return false;
  quads = vector<vec4i>{};
  quads.reserve(sizes.size());
  auto cur = 0;
  for (auto size : sizes) {
    if (size == 4) {
      quads.push_back({indices[cur + 0], indices[cur + 1], indices[cur + 2],
          indices[cur + 3]});
    } else {
      for (auto c = 2; c < size; c++) {
        quads.push_back({indices[cur + 0], indices[cur + c - 1],
            indices[cur + c], indices[cur + c]});
      }
    }
    cur += size;
  }
  return true;
}
bool get_faces(
    const ply_model& ply, vector<vec3i>& triangles, vector<vec4i>& quads) {
  if (has_quads(ply)) {
    return get_quads(ply, quads);
  } else {
    return get_triangles(ply, triangles);
  }
}
bool get_lines(const ply_model& ply, vector<vec2i>& lines) {
  auto indices = vector<int>{};
  auto sizes   = vector<uint8_t>{};
  if (!get_list_values(ply, "line", "vertex_indices", indices)) return false;
  if (!get_list_sizes(ply, "line", "vertex_indices", sizes)) return false;
  lines = vector<vec2i>{};
  lines.reserve(sizes.size());
  auto cur = 0;
  for (auto size : sizes) {
    for (auto c = 1; c < size; c++) {
      lines.push_back({indices[cur + c - 1], indices[cur + c]});
    }
    cur += size;
  }
  return true;
}
bool get_points(const ply_model& ply, vector<int>& values) {
  return get_list_values(ply, "point", "vertex_indices", values);
}
bool has_quads(const ply_model& ply) {
  auto sizes = vector<uint8_t>{};
  if (!get_list_sizes(ply, "face", "vertex_indices", sizes)) return false;
  for (auto size : sizes)
    if (size == 4) return true;
  return false;
}

// Add ply properties
inline ply_element& add_element(
    ply_model& ply, const string& element_name, size_t count) {
  for (auto& elem : ply.elements) {
    if (elem.name == element_name) return elem;
  }
  auto& elem = ply.elements.emplace_back();
  elem.name  = element_name;
  elem.count = count;
  return elem;
}
inline ply_property& add_property(ply_model& ply, const string& element_name,
    const string& property_name, size_t count, ply_type type, bool is_list) {
  add_element(ply, element_name, count);
  for (auto& elem : ply.elements) {
    if (elem.name != element_name) continue;
    for (auto& prop : elem.properties) {
      if (prop.name == property_name) return prop;
    }
    auto& prop   = elem.properties.emplace_back();
    prop.name    = property_name;
    prop.type    = type;
    prop.is_list = is_list;
    return prop;
  }
  throw std::invalid_argument{"should not have gotten here"};
}
template <typename T>
inline vector<T> make_vector(const T* value, size_t count, int stride) {
  auto ret = vector<T>(count);
  for (auto idx = (size_t)0; idx < count; idx++) ret[idx] = value[idx * stride];
  return ret;
}

inline bool add_values(ply_model& ply, const float* values, size_t count,
    const string& element, const string* properties, int nprops) {
  if (values == nullptr) return false;
  for (auto p = 0; p < nprops; p++) {
    add_property(ply, element, properties[p], count, ply_type::f32, false);
    auto& prop = get_property(ply, element, properties[p]);
    prop.data_f32.resize(count);
    for (auto i = 0; i < count; i++) prop.data_f32[i] = values[p + i * nprops];
  }
  return true;
}

inline bool add_values(ply_model& ply, const int* values, size_t count,
    const string& element, const string* properties, int nprops) {
  if (values == nullptr) return false;
  for (auto p = 0; p < nprops; p++) {
    add_property(ply, element, properties[p], count, ply_type::i32, false);
    auto& prop = get_property(ply, element, properties[p]);
    prop.data_i32.resize(count);
    for (auto i = 0; i < count; i++) prop.data_i32[i] = values[p + i * nprops];
  }
  return true;
}

bool add_value(ply_model& ply, const string& element, const string& property,
    const vector<float>& values) {
  if (values.empty()) return false;
  auto properties = vector{property};
  return add_values(
      ply, values.data(), values.size(), element, properties.data(), 1);
}
bool add_values(ply_model& ply, const string& element,
    const array<string, 2>& properties, const vector<vec2f>& values) {
  if (values.empty()) return false;
  return add_values(
      ply, &values.front().x, values.size(), element, properties.data(), 2);
}
bool add_values(ply_model& ply, const string& element,
    const array<string, 3>& properties, const vector<vec3f>& values) {
  if (values.empty()) return false;
  return add_values(
      ply, &values.front().x, values.size(), element, properties.data(), 3);
}
bool add_values(ply_model& ply, const string& element,
    const array<string, 4>& properties, const vector<vec4f>& values) {
  if (values.empty()) return false;
  return add_values(
      ply, &values.front().x, values.size(), element, properties.data(), 4);
}
bool add_values(ply_model& ply, const string& element,
    const array<string, 12>& properties, const vector<frame3f>& values) {
  if (values.empty()) return false;
  return add_values(ply, &values.front().x.x, values.size(), element,
      properties.data(), (int)properties.size());
}

bool add_value(ply_model& ply, const string& element, const string& property,
    const vector<int>& values) {
  if (values.empty()) return false;
  auto properties = vector{property};
  return add_values(
      ply, values.data(), values.size(), element, properties.data(), 1);
}
bool add_values(ply_model& ply, const string& element,
    const array<string, 2>& properties, const vector<vec2i>& values) {
  if (values.empty()) return false;
  return add_values(
      ply, &values.front().x, values.size(), element, properties.data(), 2);
}
bool add_values(ply_model& ply, const string& element,
    const array<string, 3>& properties, const vector<vec3i>& values) {
  if (values.empty()) return false;
  return add_values(
      ply, &values.front().x, values.size(), element, properties.data(), 3);
}
bool add_values(ply_model& ply, const string& element,
    const array<string, 4>& properties, const vector<vec4i>& values) {
  if (values.empty()) return false;
  return add_values(
      ply, &values.front().x, values.size(), element, properties.data(), 4);
}

bool add_lists(ply_model& ply, const string& element, const string& property,
    const vector<vector<int>>& values) {
  if (values.empty()) return false;
  add_property(ply, element, property, values.size(), ply_type::i32, true);
  auto& prop = get_property(ply, element, property);
  prop.data_i32.reserve(values.size() * 4);
  prop.ldata_u8.reserve(values.size());
  for (auto& value : values) {
    prop.data_i32.insert(prop.data_i32.end(), value.begin(), value.end());
    prop.ldata_u8.push_back((uint8_t)value.size());
  }
  return true;
}
bool add_lists(ply_model& ply, const string& element, const string& property,
    const vector<byte>& sizes, const vector<int>& values) {
  if (values.empty()) return false;
  add_property(ply, element, property, sizes.size(), ply_type::i32, true);
  auto& prop    = get_property(ply, element, property);
  prop.data_i32 = values;
  prop.ldata_u8 = sizes;
  return true;
}
bool add_lists(ply_model& ply, const int* values, size_t count, int size,
    const string& element, const string& property) {
  if (values == nullptr) return false;
  add_property(ply, element, property, count, ply_type::i32, true);
  auto& prop = get_property(ply, element, property);
  prop.data_i32.assign(values, values + count * size);
  prop.ldata_u8.assign(count, (byte)size);
  return true;
}
bool add_lists(ply_model& ply, const string& element, const string& property,
    const vector<int>& values) {
  if (values.empty()) return false;
  return add_lists(ply, values.data(), values.size(), 1, element, property);
}
bool add_lists(ply_model& ply, const string& element, const string& property,
    const vector<vec2i>& values) {
  if (values.empty()) return false;
  return add_lists(ply, &values.front().x, values.size(), 2, element, property);
}
bool add_lists(ply_model& ply, const string& element, const string& property,
    const vector<vec3i>& values) {
  if (values.empty()) return false;
  return add_lists(ply, &values.front().x, values.size(), 3, element, property);
}
bool add_lists(ply_model& ply, const string& element, const string& property,
    const vector<vec4i>& values) {
  if (values.empty()) return false;
  return add_lists(ply, &values.front().x, values.size(), 4, element, property);
}

// Add ply properties for meshes
bool add_positions(ply_model& ply, const vector<vec3f>& values) {
  return add_values(ply, "vertex", {"x", "y", "z"}, values);
}
bool add_normals(ply_model& ply, const vector<vec3f>& values) {
  return add_values(ply, "vertex", {"nx", "ny", "nz"}, values);
}
bool add_texcoords(ply_model& ply, const vector<vec2f>& values, bool flipv) {
  return add_values(
      ply, "vertex", {"u", "v"}, flipv ? flip_ply_texcoord(values) : values);
}
bool add_colors(ply_model& ply, const vector<vec3f>& values) {
  return add_values(ply, "vertex", {"red", "green", "blue"}, values);
}
bool add_colors(ply_model& ply, const vector<vec4f>& values) {
  return add_values(ply, "vertex", {"red", "green", "blue", "alpha"}, values);
}
bool add_radius(ply_model& ply, const vector<float>& values) {
  return add_value(ply, "vertex", "radius", values);
}
bool add_faces(ply_model& ply, const vector<vector<int>>& values) {
  return add_lists(ply, "face", "vertex_indices", values);
}
bool add_faces(ply_model& ply, const vector<vec3i>& triangles,
    const vector<vec4i>& quads) {
  if (triangles.empty() && quads.empty()) return false;
  if (quads.empty()) {
    return add_lists(ply, "face", "vertex_indices", triangles);
  } else if (triangles.empty() &&
             std::all_of(quads.begin(), quads.end(),
                 [](const vec4i& q) { return q.z != q.w; })) {
    return add_lists(ply, "face", "vertex_indices", quads);
  } else {
    auto sizes   = vector<uint8_t>();
    auto indices = vector<int>{};
    sizes.reserve(triangles.size() + quads.size());
    indices.reserve(triangles.size() * 3 + quads.size() * 4);
    for (auto& t : triangles) {
      sizes.push_back(3);
      indices.push_back(t.x);
      indices.push_back(t.y);
      indices.push_back(t.z);
    }
    for (auto& q : quads) {
      sizes.push_back(q.z == q.w ? 3 : 4);
      indices.push_back(q.x);
      indices.push_back(q.y);
      indices.push_back(q.z);
      if (q.z != q.w) indices.push_back(q.w);
    }
    return add_lists(ply, "face", "vertex_indices", sizes, indices);
  }
}
bool add_triangles(ply_model& ply, const vector<vec3i>& values) {
  return add_faces(ply, values, {});
}
bool add_quads(ply_model& ply, const vector<vec4i>& values) {
  return add_faces(ply, {}, values);
}
bool add_lines(ply_model& ply, const vector<vec2i>& values) {
  return add_lists(ply, "line", "vertex_indices", values);
}
bool add_points(ply_model& ply, const vector<int>& values) {
  return add_lists(ply, "point", "vertex_indices", values);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR OBJ LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto {

inline void parse_value(string_view& str, obj_vertex& value) {
  value = obj_vertex{0, 0, 0};
  parse_value(str, value.position);
  if (!str.empty() && str.front() == '/') {
    str.remove_prefix(1);
    if (!str.empty() && str.front() == '/') {
      str.remove_prefix(1);
      parse_value(str, value.normal);
    } else {
      parse_value(str, value.texcoord);
      if (!str.empty() && str.front() == '/') {
        str.remove_prefix(1);
        parse_value(str, value.normal);
      }
    }
  }
}

// Input for OBJ textures
inline void parse_value(string_view& str, obj_texture& info) {
  // initialize
  info = obj_texture();

  // get tokens
  auto tokens = vector<string>();
  skip_whitespace(str);
  while (!str.empty()) {
    auto token = ""s;
    parse_value(str, token);
    tokens.push_back(token);
    skip_whitespace(str);
  }
  if (tokens.empty()) throw std::invalid_argument{"string excepted"};

  // texture name
  info.path = tokens.back();
  for (auto& c : info.path)
    if (c == '\\') c = '/';

  // texture params
  auto last = string();
  for (auto i = 0; i < tokens.size() - 1; i++) {
    if (tokens[i] == "-bm") info.scale = (float)atof(tokens[i + 1].c_str());
    if (tokens[i] == "-clamp") info.clamp = true;
  }
}

// Load and save obj
obj_model load_obj(
    const string& filename, bool face_varying, bool split_materials) {
  auto obj = obj_model{};
  load_obj(filename, obj, face_varying, split_materials);
  return obj;
}

// Load and save obj shape
obj_shape load_sobj(const string& filename, bool face_varying) {
  auto obj = obj_shape{};
  load_obj(filename, obj, face_varying);
  return obj;
}

// Read obj
static void load_mtl(const string& filename, obj_model& obj) {
  // texture map
  auto texture_map = unordered_map<string, int>{};
  auto texture_id  = 0;
  for (auto& texture : obj.textures) texture_map[texture.path] = texture_id++;
  auto parse_texture = [&texture_map, &obj](string_view& str, int& texture_id) {
    auto texture_path = obj_texture{};
    parse_value(str, texture_path);
    auto texture_it = texture_map.find(texture_path.path);
    if (texture_it == texture_map.end()) {
      auto& texture             = obj.textures.emplace_back();
      texture.path              = texture_path.path;
      texture_id                = (int)obj.textures.size() - 1;
      texture_map[texture.path] = texture_id;
    } else {
      texture_id = texture_it->second;
    }
  };

  // open file
  auto fs = open_file(filename, "rt");

  // init parsing
  obj.materials.emplace_back();

  // read the file str by str
  auto buffer = array<char, 4096>{};
  while (read_line(fs, buffer)) {
    try {
      // str
      auto str = string_view{buffer.data()};
      remove_comment(str);
      skip_whitespace(str);
      if (str.empty()) continue;

      // get command
      auto cmd = ""s;
      parse_value(str, cmd);
      if (cmd.empty()) continue;

      // grab material
      auto& material = obj.materials.back();

      // possible token values
      if (cmd == "newmtl") {
        auto& material = obj.materials.emplace_back();
        parse_value(str, material.name);
      } else if (cmd == "illum") {
        parse_value(str, material.illum);
      } else if (cmd == "Ke") {
        parse_value(str, material.emission);
      } else if (cmd == "Ka") {
        parse_value(str, material.ambient);
      } else if (cmd == "Kd") {
        parse_value(str, material.diffuse);
      } else if (cmd == "Ks") {
        parse_value(str, material.specular);
      } else if (cmd == "Kt") {
        parse_value(str, material.transmission);
      } else if (cmd == "Tf") {
        parse_value(str, material.transmission);
        material.transmission = max(1 - material.transmission, 0.0f);
        if (max(material.transmission) < 0.001)
          material.transmission = {0, 0, 0};
      } else if (cmd == "Tr") {
        parse_value(str, material.opacity);
        material.opacity = 1 - material.opacity;
      } else if (cmd == "Ns") {
        parse_value(str, material.exponent);
      } else if (cmd == "d") {
        parse_value(str, material.opacity);
      } else if (cmd == "map_Ke") {
        parse_texture(str, material.emission_tex);
      } else if (cmd == "map_Ka") {
        parse_texture(str, material.ambient_tex);
      } else if (cmd == "map_Kd") {
        parse_texture(str, material.diffuse_tex);
      } else if (cmd == "map_Ks") {
        parse_texture(str, material.specular_tex);
      } else if (cmd == "map_Tr") {
        parse_texture(str, material.transmission_tex);
      } else if (cmd == "map_d" || cmd == "map_Tr") {
        parse_texture(str, material.opacity_tex);
      } else if (cmd == "map_bump" || cmd == "bump") {
        parse_texture(str, material.bump_tex);
      } else if (cmd == "map_disp" || cmd == "disp") {
        parse_texture(str, material.displacement_tex);
      } else if (cmd == "map_norm" || cmd == "norm") {
        parse_texture(str, material.normal_tex);
      } else {
        continue;
      }
    } catch (...) {
      throw io_error::parse_error(filename);
    }
  }

  // remove placeholder material
  obj.materials.erase(obj.materials.begin());
}

// Read obj
static void load_obx(const string& filename, obj_model& obj) {
  // texture map
  auto texture_map = unordered_map<string, int>{};
  auto texture_id  = 0;
  for (auto& texture : obj.textures) texture_map[texture.path] = texture_id++;
  auto parse_texture = [&texture_map, &obj](string_view& str, int& texture_id) {
    auto texture_path = obj_texture{};
    parse_value(str, texture_path);
    auto texture_it = texture_map.find(texture_path.path);
    if (texture_it == texture_map.end()) {
      auto& texture             = obj.textures.emplace_back();
      texture.path              = texture_path.path;
      texture_id                = (int)obj.textures.size() - 1;
      texture_map[texture.path] = texture_id;
    } else {
      texture_id = texture_it->second;
    }
  };

  // open file
  auto fs = open_file(filename, "rt");

  // init parsing
  obj.cameras.emplace_back();
  obj.environments.emplace_back();

  // read the file str by str
  auto buffer = array<char, 4096>{};
  while (read_line(fs, buffer)) {
    try {
      // str
      auto str = string_view{buffer.data()};
      remove_comment(str);
      skip_whitespace(str);
      if (str.empty()) continue;

      // get command
      auto cmd = ""s;
      parse_value(str, cmd);
      if (cmd.empty()) continue;

      // grab elements
      auto& camera      = obj.cameras.back();
      auto& environment = obj.environments.back();

      // read values
      if (cmd == "newCam") {
        auto& camera = obj.cameras.emplace_back();
        parse_value(str, camera.name);
      } else if (cmd == "Co") {
        parse_value(str, camera.ortho);
      } else if (cmd == "Ca") {
        parse_value(str, camera.aspect);
      } else if (cmd == "Cl") {
        parse_value(str, camera.lens);
      } else if (cmd == "Cs") {
        parse_value(str, camera.film);
      } else if (cmd == "Cf") {
        parse_value(str, camera.focus);
      } else if (cmd == "Cp") {
        parse_value(str, camera.aperture);
      } else if (cmd == "Cx") {
        parse_value(str, camera.frame);
      } else if (cmd == "Ct") {
        auto lookat = mat3f{};
        parse_value(str, lookat);
        camera.frame = lookat_frame(lookat.x, lookat.y, lookat.z);
        if (camera.focus == 0) camera.focus = length(lookat.y - lookat.x);
      } else if (cmd == "newEnv") {
        auto& environment = obj.environments.emplace_back();
        parse_value(str, environment.name);
      } else if (cmd == "Ee") {
        parse_value(str, environment.emission);
      } else if (cmd == "map_Ee") {
        parse_texture(str, environment.emission_tex);
      } else if (cmd == "Ex") {
        parse_value(str, environment.frame);
      } else if (cmd == "Et") {
        auto lookat = mat3f{};
        parse_value(str, lookat);
        environment.frame = lookat_frame(lookat.x, lookat.y, lookat.z, true);
      } else {
        // unused
      }
    } catch (...) {
      throw io_error::parse_error(filename);
    }
  }

  // remove placeholders
  obj.cameras.erase(obj.cameras.begin());
  obj.environments.erase(obj.environments.begin());
}

// Read obj
void load_obj(const string& filename, obj_model& obj, bool face_varying,
    bool split_materials) {
  // open file
  auto fs = open_file(filename, "rt");

  // parsing state
  auto opositions   = vector<vec3f>{};
  auto onormals     = vector<vec3f>{};
  auto otexcoords   = vector<vec2f>{};
  auto oname        = ""s;
  auto gname        = ""s;
  auto mtllibs      = vector<string>{};
  auto material_map = unordered_map<string, int>{};
  int  cur_material = -1;
  auto cur_shape    = (obj_shape*)nullptr;
  auto cur_shapes   = unordered_map<int, int>{};

  // initialize obj
  obj       = {};
  cur_shape = &obj.shapes.emplace_back();
  if (split_materials) {
    cur_shapes = {{cur_material, (int)obj.shapes.size() - 1}};
  }

  // read the file str by str
  auto buffer = array<char, 4096>{};
  while (read_line(fs, buffer)) {
    try {
      // str
      auto str = string_view{buffer.data()};
      remove_comment(str);
      skip_whitespace(str);
      if (str.empty()) continue;

      // get command
      auto cmd = ""s;
      parse_value(str, cmd);
      if (cmd.empty()) continue;

      // possible token values
      if (cmd == "v") {
        parse_value(str, opositions.emplace_back());
      } else if (cmd == "vn") {
        parse_value(str, onormals.emplace_back());
      } else if (cmd == "vt") {
        parse_value(str, otexcoords.emplace_back());
      } else if (cmd == "f" || cmd == "l" || cmd == "p") {
        // elemnet type
        auto etype = (cmd == "f")   ? obj_etype::face
                     : (cmd == "l") ? obj_etype::line
                                    : obj_etype::point;
        // grab shape and add element
        auto& shape   = *cur_shape;
        auto& element = shape.elements.emplace_back();
        if (cur_material < 0) {
          auto& material              = obj.materials.emplace_back();
          material.name               = "__default__";
          material.diffuse            = {0.8f, 0.8f, 0.8f};
          cur_material                = 0;
          material_map[material.name] = 0;
        }
        element.material = cur_material;
        element.etype    = etype;
        // parse vertices
        skip_whitespace(str);
        while (!str.empty()) {
          auto vert = obj_vertex{};
          parse_value(str, vert);
          if (vert.position == 0) break;
          if (vert.position < 0)
            vert.position = (int)opositions.size() + vert.position + 1;
          if (vert.texcoord < 0)
            vert.texcoord = (int)otexcoords.size() + vert.texcoord + 1;
          if (vert.normal < 0)
            vert.normal = (int)onormals.size() + vert.normal + 1;
          shape.vertices.push_back(vert);
          element.size += 1;
          skip_whitespace(str);
        }
      } else if (cmd == "o" || cmd == "g") {
        skip_whitespace(str);
        auto& name = cmd == "o" ? oname : gname;
        if (str.empty()) {
          name = "";
        } else {
          parse_value(str, name);
        }
        if (split_materials) {
          cur_shape       = &obj.shapes.emplace_back();
          cur_shapes      = {{cur_material, (int)obj.shapes.size() - 1}};
          cur_shape->name = oname + gname;
        } else {
          if (!cur_shape->vertices.empty()) {
            cur_shape = &obj.shapes.emplace_back();
          }
          cur_shape->name = oname + gname;
        }
      } else if (cmd == "usemtl") {
        auto mname = string{};
        parse_value(str, mname);
        auto material_it = material_map.find(mname);
        if (material_it == material_map.end())
          throw io_error::material_error(filename, mname);
        if (split_materials && cur_material != material_it->second) {
          cur_material  = material_it->second;
          auto shape_it = cur_shapes.find(cur_material);
          if (shape_it == cur_shapes.end()) {
            cur_shape                = &obj.shapes.emplace_back();
            cur_shapes[cur_material] = (int)obj.shapes.size() - 1;
            cur_shape->name          = oname + gname;
          } else {
            cur_shape = &obj.shapes.at(shape_it->second);
          }
        } else {
          cur_material = material_it->second;
        }
      } else if (cmd == "mtllib") {
        auto mtllib = ""s;
        parse_value(str, mtllib);
        if (std::find(mtllibs.begin(), mtllibs.end(), mtllib) ==
            mtllibs.end()) {
          mtllibs.push_back(mtllib);
          try {
            load_mtl(path_join(path_dirname(filename), mtllib), obj);
          } catch (const io_error& error) {
            throw io_error::dependent_error(filename, error);
          }
          auto material_id = 0;
          for (auto& material : obj.materials)
            material_map[material.name] = material_id++;
        }
      } else {
        // unused
      }
    } catch (const io_error&) {
      throw;
    } catch (...) {
      throw io_error::parse_error(filename);
    }
  }

  // remove empty shapes if splitting by materials
  if (split_materials) {
    obj.shapes.erase(
        std::remove_if(obj.shapes.begin(), obj.shapes.end(),
            [](const obj_shape& shape) { return shape.elements.empty(); }),
        obj.shapes.end());
  }

  // convert vertex data
  if (face_varying) {
    auto ipositions = vector<int>{};
    auto inormals   = vector<int>{};
    auto itexcoords = vector<int>{};
    for (auto& shape : obj.shapes) {
      ipositions.assign(opositions.size() + 1, 0);
      inormals.assign(onormals.size() + 1, 0);
      itexcoords.assign(otexcoords.size() + 1, 0);
      for (auto& vertex : shape.vertices) {
        if (vertex.position != 0 && ipositions[vertex.position] == 0) {
          shape.positions.push_back(opositions[vertex.position - 1]);
          ipositions[vertex.position] = (int)shape.positions.size();
        }
        if (vertex.normal != 0 && inormals[vertex.normal] == 0) {
          shape.normals.push_back(onormals[vertex.normal - 1]);
          inormals[vertex.normal] = (int)shape.normals.size();
        }
        if (vertex.texcoord != 0 && itexcoords[vertex.texcoord] == 0) {
          shape.texcoords.push_back(otexcoords[vertex.texcoord - 1]);
          itexcoords[vertex.texcoord] = (int)shape.texcoords.size();
        }
        vertex.position = ipositions[vertex.position];
        vertex.normal   = inormals[vertex.normal];
        vertex.texcoord = itexcoords[vertex.texcoord];
      }
    }
  } else {
    auto vertex_map = unordered_map<obj_vertex, obj_vertex>{};
    for (auto& shape : obj.shapes) {
      vertex_map.clear();
      for (auto& vertex : shape.vertices) {
        auto vertex_it = vertex_map.find(vertex);
        if (vertex_it == vertex_map.end()) {
          auto new_vertex = vertex;
          auto index      = (int)vertex_map.size();
          if (vertex.position > 0) {
            shape.positions.push_back(opositions[vertex.position - 1]);
            new_vertex.position = index + 1;
          }
          if (vertex.normal > 0) {
            shape.normals.push_back(onormals[vertex.normal - 1]);
            new_vertex.normal = index + 1;
          }
          if (vertex.texcoord > 0) {
            shape.texcoords.push_back(otexcoords[vertex.texcoord - 1]);
            new_vertex.texcoord = index + 1;
          }
          vertex_map[vertex] = new_vertex;
          vertex             = new_vertex;
        } else {
          vertex = vertex_it->second;
        }
      }
    }
  }

  // load extensions
  auto extfilename = replace_extension(filename, ".obx");
  if (path_exists(extfilename)) {
    try {
      load_obx(extfilename, obj);
    } catch (const io_error& error) {
      throw io_error::dependent_error(filename, error);
    }
  }
}

// Read obj
void load_obj(const string& filename, obj_shape& shape, bool face_varying) {
  // open file
  auto fs = open_file(filename, "rt");

  // parsing state
  auto material_map = unordered_map<string, int>{};
  int  cur_material = -1;

  // initialize obj
  shape = {};

  // read the file str by str
  auto buffer = array<char, 4096>{};
  while (read_line(fs, buffer)) {
    try {
      // str
      auto str = string_view{buffer.data()};
      remove_comment(str);
      skip_whitespace(str);
      if (str.empty()) continue;

      // get command
      auto cmd = ""s;
      parse_value(str, cmd);
      if (cmd.empty()) continue;

      // possible token values
      if (cmd == "v") {
        parse_value(str, shape.positions.emplace_back());
      } else if (cmd == "vn") {
        parse_value(str, shape.normals.emplace_back());
      } else if (cmd == "vt") {
        parse_value(str, shape.texcoords.emplace_back());
      } else if (cmd == "f" || cmd == "l" || cmd == "p") {
        // elemnet type
        auto etype = (cmd == "f")   ? obj_etype::face
                     : (cmd == "l") ? obj_etype::line
                                    : obj_etype::point;
        // grab shape and add element
        auto& element    = shape.elements.emplace_back();
        element.material = cur_material;
        element.etype    = etype;
        // parse vertices
        skip_whitespace(str);
        while (!str.empty()) {
          auto vert = obj_vertex{};
          parse_value(str, vert);
          if (vert.position == 0) break;
          if (vert.position < 0)
            vert.position = (int)shape.positions.size() + vert.position + 1;
          if (vert.texcoord < 0)
            vert.texcoord = (int)shape.texcoords.size() + vert.texcoord + 1;
          if (vert.normal < 0)
            vert.normal = (int)shape.normals.size() + vert.normal + 1;
          shape.vertices.push_back(vert);
          element.size += 1;
          skip_whitespace(str);
        }
      } else if (cmd == "usemtl") {
        auto mname = string{};
        parse_value(str, mname);
        auto material_it = material_map.find(mname);
        if (material_it == material_map.end()) {
          cur_material        = (int)material_map.size();
          material_map[mname] = cur_material;
        } else {
          cur_material = material_it->second;
        }
      } else {
        // unused
      }
    } catch (...) {
      throw io_error::parse_error(filename);
    }
  }

  // convert vertex data
  if (!face_varying) {
    auto opositions = vector<vec3f>{};
    auto onormals   = vector<vec3f>{};
    auto otexcoords = vector<vec2f>{};
    shape.positions.swap(opositions);
    shape.normals.swap(onormals);
    shape.texcoords.swap(otexcoords);
    auto vertex_map = unordered_map<obj_vertex, obj_vertex>{};
    for (auto& vertex : shape.vertices) {
      auto vertex_it = vertex_map.find(vertex);
      if (vertex_it == vertex_map.end()) {
        auto new_vertex = vertex;
        auto index      = (int)vertex_map.size();
        if (vertex.position > 0) {
          shape.positions.push_back(opositions[vertex.position - 1]);
          new_vertex.position = index + 1;
        }
        if (vertex.normal > 0) {
          shape.normals.push_back(onormals[vertex.normal - 1]);
          new_vertex.normal = index + 1;
        }
        if (vertex.texcoord > 0) {
          shape.texcoords.push_back(otexcoords[vertex.texcoord - 1]);
          new_vertex.texcoord = index + 1;
        }
        vertex_map[vertex] = new_vertex;
        vertex             = new_vertex;
      } else {
        vertex = vertex_it->second;
      }
    }
  }
}

// Format values
inline void format_value(string& str, const obj_texture& value) {
  str += value.path.empty() ? "" : value.path;
}
inline void format_value(string& str, const obj_vertex& value) {
  format_value(str, value.position);
  if (value.texcoord != 0) {
    str += "/";
    format_value(str, value.texcoord);
    if (value.normal != 0) {
      str += "/";
      format_value(str, value.normal);
    }
  } else if (value.normal != 0) {
    str += "//";
    format_value(str, value.normal);
  }
}

// Save obj
inline void save_mtl(const string& filename, const obj_model& obj) {
  // open file
  auto fs = open_file(filename, "wt");

  // save comments
  format_values(fs, "#\n");
  format_values(fs, "# Written by Yocto/GL\n");
  format_values(fs, "# https://github.com/xelatihy/yocto-gl\n");
  format_values(fs, "#\n\n");
  for (auto& comment : obj.comments) {
    format_values(fs, "# {}\n", comment);
  }
  format_values(fs, "\n");

  // write material
  for (auto& material : obj.materials) {
    format_values(fs, "newmtl {}\n", material.name);
    format_values(fs, "illum {}\n", material.illum);
    if (material.emission != zero3f)
      format_values(fs, "Ke {}\n", material.emission);
    if (material.ambient != zero3f)
      format_values(fs, "Ka {}\n", material.ambient);
    format_values(fs, "Kd {}\n", material.diffuse);
    format_values(fs, "Ks {}\n", material.specular);
    if (material.reflection != zero3f)
      format_values(fs, "Kr {}\n", material.reflection);
    if (material.transmission != zero3f)
      format_values(fs, "Kt {}\n", material.transmission);
    format_values(fs, "Ns {}\n", (int)material.exponent);
    if (material.opacity != 1) format_values(fs, "d {}\n", material.opacity);
    if (material.emission_tex >= 0)
      format_values(
          fs, "map_Ke {}\n", obj.textures[material.emission_tex].path);
    if (material.diffuse_tex >= 0)
      format_values(fs, "map_Kd {}\n", obj.textures[material.diffuse_tex].path);
    if (material.specular_tex >= 0)
      format_values(
          fs, "map_Ks {}\n", obj.textures[material.specular_tex].path);
    if (material.transmission_tex >= 0)
      format_values(
          fs, "map_Kt {}\n", obj.textures[material.transmission_tex].path);
    if (material.reflection_tex >= 0)
      format_values(
          fs, "map_Kr {}\n", obj.textures[material.reflection_tex].path);
    if (material.exponent_tex >= 0)
      format_values(
          fs, "map_Ns {}\n", obj.textures[material.exponent_tex].path);
    if (material.opacity_tex >= 0)
      format_values(fs, "map_d {}\n", obj.textures[material.opacity_tex].path);
    if (material.bump_tex >= 0)
      format_values(fs, "map_bump {}\n", obj.textures[material.bump_tex].path);
    if (material.displacement_tex >= 0)
      format_values(
          fs, "map_disp {}\n", obj.textures[material.displacement_tex].path);
    if (material.normal_tex >= 0)
      format_values(
          fs, "map_norm {}\n", obj.textures[material.normal_tex].path);
    format_values(fs, "\n");
  }
}

// Save obj
inline void save_obx(const string& filename, const obj_model& obj) {
  // open file
  auto fs = open_file(filename, "wt");

  // save comments
  format_values(fs, "#\n");
  format_values(fs, "# Written by Yocto/GL\n");
  format_values(fs, "# https://github.com/xelatihy/yocto-gl\n");
  format_values(fs, "#\n\n");
  for (auto& comment : obj.comments) {
    format_values(fs, "# {}\n", comment);
  }
  format_values(fs, "\n");

  // cameras
  for (auto& camera : obj.cameras) {
    format_values(fs, "newCam {}\n", camera.name);
    format_values(fs, "  Co {}\n", camera.ortho);
    format_values(fs, "  Ca {}\n", camera.aspect);
    format_values(fs, "  Cl {}\n", camera.lens);
    format_values(fs, "  Cs {}\n", camera.film);
    format_values(fs, "  Cf {}\n", camera.focus);
    format_values(fs, "  Cp {}\n", camera.aperture);
    format_values(fs, "  Cx {}\n", camera.frame);
  }

  // environments
  for (auto& environment : obj.environments) {
    format_values(fs, "newEnv {}\n", environment.name);
    format_values(fs, "  Ee {}\n", environment.emission);
    if (environment.emission_tex >= 0) {
      format_values(
          fs, "  map_Ee {}\n", obj.textures[environment.emission_tex].path);
    }
    format_values(fs, "  Ex {}\n", environment.frame);
  }
}

// Save obj
void save_obj(const string& filename, const obj_model& obj) {
  // open file
  auto fs = open_file(filename, "wt");

  // save comments
  format_values(fs, "#\n");
  format_values(fs, "# Written by Yocto/GL\n");
  format_values(fs, "# https://github.com/xelatihy/yocto-gl\n");
  format_values(fs, "#\n\n");
  for (auto& comment : obj.comments) {
    format_values(fs, "# {}\n", comment);
  }
  format_values(fs, "\n");

  // save material library
  if (!obj.materials.empty()) {
    format_values(fs, "mtllib {}\n\n",
        replace_extension(path_filename(filename), ".mtl"));
  }

  // save objects
  auto vert_size = obj_vertex{0, 0, 0};
  for (auto& shape : obj.shapes) {
    format_values(fs, "o {}\n", shape.name);
    for (auto& p : shape.positions) format_values(fs, "v {}\n", p);
    for (auto& n : shape.normals) format_values(fs, "vn {}\n", n);
    for (auto& t : shape.texcoords) format_values(fs, "vt {}\n", t);
    auto cur_material = -1, cur_vertex = 0;
    for (auto& element : shape.elements) {
      if (!obj.materials.empty() && cur_material != element.material) {
        format_values(fs, "usemtl {}\n", obj.materials[element.material].name);
        cur_material = element.material;
      }
      if (element.etype == obj_etype::face) {
        format_values(fs, "{}", "f");
      } else if (element.etype == obj_etype::line) {
        format_values(fs, "{}", "l");
      } else if (element.etype == obj_etype::point) {
        format_values(fs, "{}", "p");
      }
      for (auto c = 0; c < element.size; c++) {
        auto vert = shape.vertices[cur_vertex++];
        if (vert.position != 0) vert.position += vert_size.position;
        if (vert.normal != 0) vert.normal += vert_size.normal;
        if (vert.texcoord != 0) vert.texcoord += vert_size.texcoord;
        format_values(fs, " {}", vert);
      }
      format_values(fs, "\n");
    }
    format_values(fs, "\n");
    vert_size.position += (int)shape.positions.size();
    vert_size.normal += (int)shape.normals.size();
    vert_size.texcoord += (int)shape.texcoords.size();
  }

  // save mtl
  if (!obj.materials.empty()) {
    try {
      save_mtl(replace_extension(filename, ".mtl"), obj);
    } catch (const io_error& error) {
      throw io_error::dependent_error(filename, error);
    }
  }

  // save obx
  if (!obj.cameras.empty() || !obj.environments.empty()) {
    try {
      save_obx(replace_extension(filename, ".obx"), obj);
    } catch (const io_error& error) {
      throw io_error::dependent_error(filename, error);
    }
  }
}

// Save obj
void save_obj(const string& filename, const obj_shape& shape) {
  // open file
  auto fs = open_file(filename, "wt");

  // save comments
  format_values(fs, "#\n");
  format_values(fs, "# Written by Yocto/GL\n");
  format_values(fs, "# https://github.com/xelatihy/yocto-gl\n");
  format_values(fs, "#\n\n");
  format_values(fs, "\n");

  // save objects
  format_values(fs, "o {}\n", shape.name);
  for (auto& p : shape.positions) format_values(fs, "v {}\n", p);
  for (auto& n : shape.normals) format_values(fs, "vn {}\n", n);
  for (auto& t : shape.texcoords) format_values(fs, "vt {}\n", t);
  auto cur_material = -1, cur_vertex = 0;
  for (auto& element : shape.elements) {
    if (cur_material != element.material) {
      format_values(
          fs, "usemtl {}\n", "material" + std::to_string(element.material));
      cur_material = element.material;
    }
    if (element.etype == obj_etype::face) {
      format_values(fs, "{}", "f");
    } else if (element.etype == obj_etype::line) {
      format_values(fs, "{}", "l");
    } else if (element.etype == obj_etype::point) {
      format_values(fs, "{}", "p");
    }
    for (auto c = 0; c < element.size; c++) {
      auto& vert = shape.vertices[cur_vertex++];
      format_values(fs, " {}", vert);
    }
    format_values(fs, "\n");
  }
}

// Load and save obj
bool load_obj(const string& filename, obj_model& obj, string& error,
    bool face_varying, bool split_materials) {
  try {
    load_obj(filename, obj, face_varying, split_materials);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}
bool save_obj(const string& filename, const obj_model& obj, string& error) {
  try {
    save_obj(filename, obj);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Load and save obj shape
bool load_obj(
    const string& filename, obj_shape& obj, string& error, bool face_varying) {
  try {
    load_obj(filename, obj, face_varying);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}
bool save_obj(const string& filename, const obj_shape& obj, string& error) {
  try {
    save_obj(filename, obj);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Get obj shape.
void get_positions(const obj_shape& shape, vector<vec3f>& positions) {
  positions = shape.positions;
}
void get_normals(const obj_shape& shape, vector<vec3f>& normals) {
  normals = shape.normals;
}
void get_texcoords(
    const obj_shape& shape, vector<vec2f>& texcoords, bool flipv) {
  texcoords = shape.texcoords;
  if (flipv) {
    for (auto& texcoord : texcoords) texcoord.y = 1 - texcoord.y;
  }
}
void get_faces(const obj_shape& shape, vector<vec3i>& triangles,
    vector<vec4i>& quads, vector<int>& materials) {
  if (has_quads(shape)) {
    get_quads(shape, quads, materials);
  } else {
    get_triangles(shape, triangles, materials);
  }
}
void get_triangles(
    const obj_shape& shape, vector<vec3i>& triangles, vector<int>& materials) {
  triangles.clear();
  materials.clear();
  triangles.reserve(shape.elements.size());
  materials.reserve(shape.elements.size());
  auto cur = 0;
  for (auto& element : shape.elements) {
    if (element.etype != obj_etype::face) continue;
    for (auto c = 2; c < element.size; c++) {
      triangles.push_back({shape.vertices[cur + 0].position - 1,
          shape.vertices[cur + c - 1].position - 1,
          shape.vertices[cur + c].position - 1});
      materials.push_back(element.material);
    }
    cur += element.size;
  }
}
void get_quads(
    const obj_shape& shape, vector<vec4i>& quads, vector<int>& materials) {
  quads.clear();
  materials.clear();
  quads.reserve(shape.elements.size());
  materials.reserve(shape.elements.size());
  auto cur = 0;
  for (auto& element : shape.elements) {
    if (element.etype != obj_etype::face) continue;
    if (element.size == 4) {
      quads.push_back({shape.vertices[cur + 0].position - 1,
          shape.vertices[cur + 1].position - 1,
          shape.vertices[cur + 2].position - 1,
          shape.vertices[cur + 3].position - 1});
      materials.push_back(element.material);
    } else {
      for (auto c = 2; c < element.size; c++) {
        quads.push_back({shape.vertices[cur + 0].position - 1,
            shape.vertices[cur + c - 1].position - 1,
            shape.vertices[cur + c].position - 1,
            shape.vertices[cur + c].position - 1});
        materials.push_back(element.material);
      }
    }
    cur += element.size;
  }
}
void get_lines(
    const obj_shape& shape, vector<vec2i>& lines, vector<int>& materials) {
  lines.clear();
  materials.clear();
  lines.reserve(shape.elements.size());
  materials.reserve(shape.elements.size());
  auto cur = 0;
  for (auto& element : shape.elements) {
    if (element.etype != obj_etype::line) continue;
    for (auto c = 1; c < element.size; c++) {
      lines.push_back({shape.vertices[cur + c - 1].position - 1,
          shape.vertices[cur + c].position - 1});
      materials.push_back(element.material);
    }
    cur += element.size;
  }
}
void get_points(
    const obj_shape& shape, vector<int>& points, vector<int>& materials) {
  points.clear();
  materials.clear();
  points.reserve(shape.elements.size());
  materials.reserve(shape.elements.size());
  auto cur = 0;
  for (auto& element : shape.elements) {
    if (element.etype != obj_etype::point) continue;
    for (auto c = 0; c < element.size; c++) {
      points.push_back({shape.vertices[cur + 0].position - 1});
      materials.push_back(element.material);
    }
    cur += element.size;
  }
}
void get_fvquads(const obj_shape& shape, vector<vec4i>& quadspos,
    vector<vec4i>& quadsnorm, vector<vec4i>& quadstexcoord,
    vector<int>& materials) {
  quadspos.clear();
  quadsnorm.clear();
  quadstexcoord.clear();
  materials.clear();
  quadspos.reserve(shape.elements.size());
  quadsnorm.reserve(shape.elements.size());
  quadstexcoord.reserve(shape.elements.size());
  materials.reserve(shape.elements.size());
  auto cur = 0;
  for (auto& element : shape.elements) {
    if (element.etype != obj_etype::face) continue;
    if (element.size == 4) {
      if (shape.vertices[0].position != 0)
        quadspos.push_back({shape.vertices[cur + 0].position - 1,
            shape.vertices[cur + 1].position - 1,
            shape.vertices[cur + 2].position - 1,
            shape.vertices[cur + 3].position - 1});
      if (shape.vertices[0].normal != 0)
        quadsnorm.push_back({shape.vertices[cur + 0].normal - 1,
            shape.vertices[cur + 1].normal - 1,
            shape.vertices[cur + 2].normal - 1,
            shape.vertices[cur + 3].normal - 1});
      if (shape.vertices[0].texcoord != 0)
        quadstexcoord.push_back({shape.vertices[cur + 0].texcoord - 1,
            shape.vertices[cur + 1].texcoord - 1,
            shape.vertices[cur + 2].texcoord - 1,
            shape.vertices[cur + 3].texcoord - 1});
      materials.push_back(element.material);
    } else {
      for (auto c = 2; c < element.size; c++) {
        if (shape.vertices[0].position != 0)
          quadspos.push_back({shape.vertices[cur + 0].position - 1,
              shape.vertices[cur + c - 1].position - 1,
              shape.vertices[cur + c].position - 1,
              shape.vertices[cur + c].position - 1});
        if (shape.vertices[0].normal != 0)
          quadsnorm.push_back({shape.vertices[cur + 0].normal - 1,
              shape.vertices[cur + c - 1].normal - 1,
              shape.vertices[cur + c].normal - 1,
              shape.vertices[cur + c].normal - 1});
        if (shape.vertices[0].texcoord != 0)
          quadstexcoord.push_back({shape.vertices[cur + 0].texcoord - 1,
              shape.vertices[cur + c - 1].texcoord - 1,
              shape.vertices[cur + c].texcoord - 1,
              shape.vertices[cur + c].texcoord - 1});
        materials.push_back(element.material);
      }
    }
    cur += element.size;
  }
}
void get_faces(const obj_shape& shape, int material, vector<vec3i>& triangles,
    vector<vec4i>& quads) {
  if (has_quads(shape)) {
    get_quads(shape, material, quads);
  } else {
    get_triangles(shape, material, triangles);
  }
}
void get_triangles(
    const obj_shape& shape, int material, vector<vec3i>& triangles) {
  triangles.clear();
  if (shape.elements.empty()) return;
  triangles.reserve(shape.elements.size());
  auto cur = 0;
  for (auto& element : shape.elements) {
    if (element.etype != obj_etype::face) continue;
    if (element.material != material) continue;
    for (auto c = 2; c < element.size; c++) {
      triangles.push_back({shape.vertices[cur + 0].position - 1,
          shape.vertices[cur + c - 1].position - 1,
          shape.vertices[cur + c].position - 1});
    }
    cur += element.size;
  }
}
void get_quads(const obj_shape& shape, int material, vector<vec4i>& quads) {
  quads.clear();
  if (shape.elements.empty()) return;
  quads.reserve(shape.elements.size());
  auto cur = 0;
  for (auto& element : shape.elements) {
    if (element.etype != obj_etype::face) continue;
    if (element.material != material) continue;
    if (element.size == 4) {
      quads.push_back({shape.vertices[cur + 0].position - 1,
          shape.vertices[cur + 1].position - 1,
          shape.vertices[cur + 2].position - 1,
          shape.vertices[cur + 3].position - 1});
    } else {
      for (auto c = 2; c < element.size; c++) {
        quads.push_back({shape.vertices[cur + 0].position - 1,
            shape.vertices[cur + c - 1].position - 1,
            shape.vertices[cur + c].position - 1,
            shape.vertices[cur + c].position - 1});
      }
    }
    cur += element.size;
  }
}
void get_lines(const obj_shape& shape, int material, vector<vec2i>& lines) {
  lines.clear();
  if (shape.elements.empty()) return;
  lines.reserve(shape.elements.size());
  auto cur = 0;
  for (auto& element : shape.elements) {
    if (element.etype != obj_etype::line) continue;
    if (element.material != material) continue;
    for (auto c = 1; c < element.size; c++) {
      lines.push_back({shape.vertices[cur + c - 1].position - 1,
          shape.vertices[cur + c].position - 1});
    }
    cur += element.size;
  }
}
void get_points(const obj_shape& shape, int material, vector<int>& points) {
  points.clear();
  if (shape.elements.empty()) return;
  points.reserve(shape.elements.size());
  auto cur = 0;
  for (auto& element : shape.elements) {
    if (element.etype != obj_etype::point) continue;
    if (element.material != material) continue;
    for (auto c = 0; c < element.size; c++) {
      points.push_back({shape.vertices[cur + 0].position - 1});
    }
    cur += element.size;
  }
}

bool has_quads(const obj_shape& shape) {
  for (auto& element : shape.elements)
    if (element.etype == obj_etype::face && element.size == 4) return true;
  return false;
}

vector<int> get_materials(const obj_shape& shape) {
  auto materials    = vector<int>{};
  auto material_set = unordered_set<int>{};
  for (auto& element : shape.elements) {
    if (material_set.find(element.material) == material_set.end()) {
      material_set.insert(element.material);
      materials.push_back(element.material);
    }
  }
  return materials;
}

// Add obj shape
void add_positions(obj_shape& shape, const vector<vec3f>& positions) {
  shape.positions.insert(
      shape.positions.end(), positions.begin(), positions.end());
}
void add_normals(obj_shape& shape, const vector<vec3f>& normals) {
  shape.normals.insert(shape.normals.end(), normals.begin(), normals.end());
}
void add_texcoords(
    obj_shape& shape, const vector<vec2f>& texcoords, bool flipv) {
  shape.texcoords.insert(
      shape.texcoords.end(), texcoords.begin(), texcoords.end());
  if (flipv) {
    for (auto idx = shape.texcoords.size() - texcoords.size();
         idx < shape.texcoords.size(); idx++)
      shape.texcoords[idx].y = 1 - shape.texcoords[idx].y;
  }
}
void add_triangles(obj_shape& shape, const vector<vec3i>& triangles,
    int material, bool has_normals, bool has_texcoord) {
  for (auto idx = 0; idx < triangles.size(); idx++) {
    auto& triangle = triangles[idx];
    for (auto c = 0; c < 3; c++) {
      shape.vertices.push_back({
          triangle[c] + 1,
          !has_texcoord ? 0 : triangle[c] + 1,
          !has_normals ? 0 : triangle[c] + 1,
      });
    }
    shape.elements.push_back({3, obj_etype::face, material});
  }
}
void add_quads(obj_shape& shape, const vector<vec4i>& quads, int material,
    bool has_normals, bool has_texcoord) {
  for (auto idx = 0; idx < quads.size(); idx++) {
    auto& quad = quads[idx];
    auto  nv   = quad.z == quad.w ? 3 : 4;
    for (auto c = 0; c < nv; c++) {
      shape.vertices.push_back({
          quad[c] + 1,
          !has_texcoord ? 0 : quad[c] + 1,
          !has_normals ? 0 : quad[c] + 1,
      });
    }
    shape.elements.push_back({(uint16_t)nv, obj_etype::face, material});
  }
}
void add_lines(obj_shape& shape, const vector<vec2i>& lines, int material,
    bool has_normals, bool has_texcoord) {
  for (auto idx = 0; idx < lines.size(); idx++) {
    auto& line = lines[idx];
    for (auto c = 0; c < 2; c++) {
      shape.vertices.push_back({
          line[c] + 1,
          !has_texcoord ? 0 : line[c] + 1,
          !has_normals ? 0 : line[c] + 1,
      });
    }
    shape.elements.push_back({2, obj_etype::line, material});
  }
}
void add_points(obj_shape& shape, const vector<int>& points, int material,
    bool has_normals, bool has_texcoord) {
  for (auto idx = 0; idx < points.size(); idx++) {
    auto& point = points[idx];
    shape.vertices.push_back({
        point + 1,
        !has_texcoord ? 0 : point + 1,
        !has_normals ? 0 : point + 1,
    });
    shape.elements.push_back({1, obj_etype::point, material});
  }
}
void add_fvquads(obj_shape& shape, const vector<vec4i>& quadspos,
    const vector<vec4i>& quadsnorm, const vector<vec4i>& quadstexcoord,
    int material) {
  for (auto idx = 0; idx < quadspos.size(); idx++) {
    auto nv = quadspos[idx].z == quadspos[idx].w ? 3 : 4;
    for (auto c = 0; c < nv; c++) {
      shape.vertices.push_back({
          quadspos.empty() ? 0 : quadspos[idx][c] + 1,
          quadstexcoord.empty() ? 0 : quadstexcoord[idx][c] + 1,
          quadsnorm.empty() ? 0 : quadsnorm[idx][c] + 1,
      });
    }
    shape.elements.push_back({(uint16_t)nv, obj_etype::face, material});
  }
}
void add_quads(obj_shape& shape, const vector<vec4i>& quads,
    const vector<int>& materials, bool has_normals, bool has_texcoord) {
  for (auto idx = 0; idx < quads.size(); idx++) {
    auto& quad = quads[idx];
    auto  nv   = quad.z == quad.w ? 3 : 4;
    for (auto c = 0; c < nv; c++) {
      shape.vertices.push_back({
          quad[c] + 1,
          !has_texcoord ? 0 : quad[c] + 1,
          !has_normals ? 0 : quad[c] + 1,
      });
    }
    shape.elements.push_back({(uint16_t)nv, obj_etype::face, materials[idx]});
  }
}
void add_lines(obj_shape& shape, const vector<vec2i>& lines,
    const vector<int>& materials, bool has_normals, bool has_texcoord) {
  for (auto idx = 0; idx < lines.size(); idx++) {
    auto& line = lines[idx];
    for (auto c = 0; c < 2; c++) {
      shape.vertices.push_back({
          line[c] + 1,
          !has_texcoord ? 0 : line[c] + 1,
          !has_normals ? 0 : line[c] + 1,
      });
    }
    shape.elements.push_back({2, obj_etype::line, materials[idx]});
  }
}
void add_points(obj_shape& shape, const vector<int>& points,
    const vector<int>& materials, bool has_normals, bool has_texcoord) {
  for (auto idx = 0; idx < points.size(); idx++) {
    auto& point = points[idx];
    shape.vertices.push_back({
        point + 1,
        !has_texcoord ? 0 : point + 1,
        !has_normals ? 0 : point + 1,
    });
    shape.elements.push_back({1, obj_etype::point, materials[idx]});
  }
}
void add_fvquads(obj_shape& shape, const vector<vec4i>& quadspos,
    const vector<vec4i>& quadsnorm, const vector<vec4i>& quadstexcoord,
    const vector<int>& materials) {
  for (auto idx = 0; idx < quadspos.size(); idx++) {
    auto nv = quadspos[idx].z == quadspos[idx].w ? 3 : 4;
    for (auto c = 0; c < nv; c++) {
      shape.vertices.push_back({
          quadspos.empty() ? 0 : quadspos[idx][c] + 1,
          quadstexcoord.empty() ? 0 : quadstexcoord[idx][c] + 1,
          quadsnorm.empty() ? 0 : quadsnorm[idx][c] + 1,
      });
    }
    shape.elements.push_back({(uint16_t)nv, obj_etype::face, materials[idx]});
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// HELPER FOR StL
// -----------------------------------------------------------------------------
namespace std {

// Hash functor for vector for use with hash_map
template <>
struct hash<yocto::vec3f> {
  size_t operator()(const yocto::vec3f& v) const {
    const std::hash<float> hasher = std::hash<float>();
    auto                   h      = (size_t)0;
    h ^= hasher(v.x) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hasher(v.y) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hasher(v.z) + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
  }
};

}  // namespace std

// -----------------------------------------------------------------------------
// STL PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Load and save stl
stl_model load_stl(const string& filename, bool unique_vertices) {
  auto stl = stl_model{};
  load_stl(filename, stl);
  return stl;
}

// Load/save stl
void load_stl(const string& filename, stl_model& stl, bool unique_vertices) {
  stl.shapes.clear();

  // open file
  auto fs = open_file(filename, "rb");

  // assume it is binary and read hader
  auto header = array<char, 80>{};
  read_value(fs, header);

  // check if binary
  auto binary = header[0] != 's' || header[1] != 'o' || header[2] != 'l' ||
                header[3] != 'i' || header[4] != 'd';

  // check size in case the binary had a bad header
  if (!binary) {
    auto ntriangles = (uint32_t)0;
    read_value(fs, ntriangles);
    fseek(fs.fs, 0, SEEK_SET);
    fseek(fs.fs, 0, SEEK_END);
    auto length = ftell(fs.fs);
    fseek(fs.fs, 0, SEEK_SET);
    auto size = 80 + 4 + (4 * 12 + 2) * (size_t)ntriangles;
    binary    = length == size;
  }

  // close file
  close_file(fs);

  // switch on type
  if (binary) {
    // open file
    auto fs = open_file(filename, "rb");

    // skip header
    auto header = array<char, 80>{};
    read_value(fs, header);

    // read shapes until the end
    auto ntriangles = (uint32_t)0;
    while (is_eof(fs)) {
      // append shape
      auto& shape = stl.shapes.emplace_back();

      // read num triangles
      read_value(fs, ntriangles);

      // resize buffers
      shape.fnormals.resize(ntriangles);
      shape.triangles.resize(ntriangles);
      shape.positions.resize(ntriangles * 3);

      // read all data
      for (auto triangle_id = 0; triangle_id < (int)ntriangles; triangle_id++) {
        // read triangle data
        read_value(fs, shape.fnormals[triangle_id]);
        read_value(fs, shape.positions[triangle_id * 3 + 0]);
        read_value(fs, shape.positions[triangle_id * 3 + 1]);
        read_value(fs, shape.positions[triangle_id * 3 + 2]);
        shape.triangles[triangle_id] = {
            triangle_id * 3 + 0, triangle_id * 3 + 1, triangle_id * 3 + 2};
        // read unused attrobute count
        auto attribute_count = (uint16_t)0;
        read_value(fs, attribute_count);
      }
    }

    // check if read at least one
    if (stl.shapes.empty()) throw io_error::read_error(filename);
  } else {
    // if ascii, re-open the file as text
    auto fs = open_file(filename, "rt");

    // parse state
    auto in_solid = false, in_facet = false, in_loop = false;
    // raed all lines
    auto buffer = array<char, 4096>{};
    while (read_line(fs, buffer)) {
      try {
        // str
        auto str = string_view{buffer.data()};
        remove_comment(str);
        skip_whitespace(str);
        if (str.empty()) continue;

        // get command
        auto cmd = ""s;
        parse_value(str, cmd);
        if (cmd.empty()) continue;

        // switch over command
        if (cmd == "solid") {
          if (in_solid) throw io_error::parse_error(filename);
          in_solid = true;
          stl.shapes.emplace_back();
        } else if (cmd == "endsolid") {
          if (!in_solid) throw io_error::parse_error(filename);
          in_solid = false;
        } else if (cmd == "facet") {
          if (!in_solid || in_facet) throw io_error::parse_error(filename);
          in_facet = true;
          // next command
          parse_value(str, cmd);
          if (cmd != "normal") throw io_error::parse_error(filename);
          // vertex normal
          parse_value(str, stl.shapes.back().fnormals.emplace_back());
        } else if (cmd == "endfacet") {
          if (!in_solid || !in_facet || in_loop)
            throw io_error::parse_error(filename);
          in_facet = false;
          // check that it was a triangle
          auto last_pos = (int)stl.shapes.back().positions.size() - 3;
          if (stl.shapes.back().triangles.empty() && last_pos != 0)
            throw io_error::parse_error(filename);
          if (!stl.shapes.back().triangles.empty() &&
              last_pos != stl.shapes.back().triangles.back().z + 1)
            throw io_error::parse_error(filename);
          // add triangle
          stl.shapes.back().triangles.push_back(
              {last_pos + 0, last_pos + 1, last_pos + 2});
        } else if (cmd == "outer") {
          if (!in_solid || !in_facet || in_loop)
            throw io_error::parse_error(filename);
          in_loop = true;
          // next command
          parse_value(str, cmd);
          if (cmd != "loop") throw io_error::parse_error(filename);
        } else if (cmd == "endloop") {
          if (!in_solid || !in_facet || !in_loop)
            throw io_error::parse_error(filename);
          in_loop = false;
        } else if (cmd == "vertex") {
          // vertex position
          parse_value(str, stl.shapes.back().positions.emplace_back());
        } else {
          throw io_error::parse_error(filename);
        }
      } catch (const io_error&) {
        throw;
      } catch (...) {
        throw io_error::parse_error(filename);
      }
    }
  }

  // make unique vertices
  if (unique_vertices) {
    for (auto& shape : stl.shapes) {
      auto vertex_map       = unordered_map<vec3f, int>{};
      auto unique_positions = vector<vec3f>{};
      for (auto& triangle : shape.triangles) {
        for (auto& vertex_id : triangle) {
          auto vertex_it = vertex_map.find(shape.positions[vertex_id]);
          if (vertex_it == vertex_map.end()) {
            auto new_vertex_id = (int)unique_positions.size();
            unique_positions.push_back(shape.positions[vertex_id]);
            vertex_map.insert(
                vertex_it, {unique_positions.back(), new_vertex_id});
            vertex_id = new_vertex_id;
          } else {
            vertex_id = vertex_it->second;
          }
        }
      }
      std::swap(unique_positions, shape.positions);
    }
  }
}

void save_stl(const string& filename, const stl_model& stl, bool ascii) {
  // helper
  auto triangle_normal = [](const vec3f& p0, const vec3f& p1, const vec3f& p2) {
    return normalize(cross(p1 - p0, p2 - p0));
  };

  // open file
  auto fs = open_file(filename, ascii ? "wt" : "wb");

  // switch on format
  if (!ascii) {
    // header
    auto header = array<char, 80>{0};
    snprintf(header.data(), header.size(), "Binary STL - Written by Yocto/GL");
    write_value(fs, header);

    // write shapes
    for (auto& shape : stl.shapes) {
      auto ntriangles = (uint32_t)shape.triangles.size();
      write_value(fs, ntriangles);
      for (auto triangle_idx = 0; triangle_idx < shape.triangles.size();
           triangle_idx++) {
        auto& triangle = shape.triangles[triangle_idx];
        auto  fnormal  = !shape.fnormals.empty()
                             ? shape.fnormals[triangle_idx]
                             : triangle_normal(shape.positions[triangle.x],
                                 shape.positions[triangle.y],
                                 shape.positions[triangle.z]);
        write_value(fs, fnormal);
        write_value(fs, shape.positions[triangle.x]);
        write_value(fs, shape.positions[triangle.y]);
        write_value(fs, shape.positions[triangle.z]);
        auto attribute_count = (uint16_t)0;
        write_value(fs, attribute_count);
      }
    }
  } else {
    for (auto& shape : stl.shapes) {
      format_values(fs, "solid \n");
      for (auto triangle_idx = 0; triangle_idx < shape.triangles.size();
           triangle_idx++) {
        auto& triangle = shape.triangles[triangle_idx];
        auto  fnormal  = !shape.fnormals.empty()
                             ? shape.fnormals[triangle_idx]
                             : triangle_normal(shape.positions[triangle.x],
                                 shape.positions[triangle.y],
                                 shape.positions[triangle.z]);
        format_values(fs, "facet normal {}\n", fnormal);
        format_values(fs, "outer loop\n");
        format_values(fs, "vertex {}\n", shape.positions[triangle.x]);
        format_values(fs, "vertex {}\n", shape.positions[triangle.y]);
        format_values(fs, "vertex {}\n", shape.positions[triangle.z]);
        format_values(fs, "endloop\n");
        format_values(fs, "endfacet\n");
      }
      format_values(fs, "endsolid \n");
    }
  }
}

// Load and save stl
bool load_stl(const string& filename, stl_model& stl, string& error,
    bool unique_vertices) {
  try {
    load_stl(filename, stl, unique_vertices);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}
bool save_stl(
    const string& filename, const stl_model& stl, string& error, bool ascii) {
  try {
    save_stl(filename, stl, ascii);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Get/set data
bool get_triangles(const stl_model& stl, int shape_id, vector<vec3i>& triangles,
    vector<vec3f>& positions, vector<vec3f>& fnormals) {
  if (shape_id < 0 || shape_id >= stl.shapes.size()) return false;
  auto& shape = stl.shapes.at(shape_id);
  triangles   = shape.triangles;
  positions   = shape.positions;
  fnormals    = shape.fnormals;
  return true;
}
void add_triangles(stl_model& stl, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& fnormals) {
  auto& shape     = stl.shapes.emplace_back();
  shape.triangles = triangles;
  shape.positions = positions;
  shape.fnormals  = fnormals;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PBRT PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Pbrt type
enum struct pbrt_type {
  // clang-format off
  real, integer, boolean, string, point, normal, vector, texture, color,
  point2, vector2, spectrum
  // clang-format on
};

// Pbrt value
struct pbrt_value {
  string        name     = "";
  pbrt_type     type     = pbrt_type::real;
  int           value1i  = 0;
  float         value1f  = 0;
  vec2f         value2f  = {0, 0};
  vec3f         value3f  = {0, 0, 0};
  bool          value1b  = false;
  string        value1s  = "";
  vector<float> vector1f = {};
  vector<vec2f> vector2f = {};
  vector<vec3f> vector3f = {};
  vector<int>   vector1i = {};
};

// Pbrt command
struct pbrt_command {
  string             name   = "";
  string             type   = "";
  vector<pbrt_value> values = {};
  frame3f            frame  = identity3x4f;
  frame3f            frend  = identity3x4f;
};

// get pbrt value
inline bool get_pbrt_value(const pbrt_value& pbrt, string& val) {
  if (pbrt.type == pbrt_type::string || pbrt.type == pbrt_type::texture) {
    val = pbrt.value1s;
    return true;
  } else {
    return false;
  }
}
inline bool get_pbrt_value(const pbrt_value& pbrt, bool& val) {
  if (pbrt.type == pbrt_type::boolean) {
    val = pbrt.value1b;
    return true;
  } else {
    return false;
  }
}
inline bool get_pbrt_value(const pbrt_value& pbrt, int& val) {
  if (pbrt.type == pbrt_type::integer) {
    val = pbrt.value1i;
    return true;
  } else {
    return false;
  }
}
inline bool get_pbrt_value(const pbrt_value& pbrt, float& val) {
  if (pbrt.type == pbrt_type::real) {
    val = pbrt.value1f;
    return true;
  } else {
    return false;
  }
}
inline bool get_pbrt_value(const pbrt_value& pbrt, vec2f& val) {
  if (pbrt.type == pbrt_type::point2 || pbrt.type == pbrt_type::vector2) {
    val = pbrt.value2f;
    return true;
  } else {
    return false;
  }
}
inline bool get_pbrt_value(const pbrt_value& pbrt, vec3f& val) {
  if (pbrt.type == pbrt_type::point || pbrt.type == pbrt_type::vector ||
      pbrt.type == pbrt_type::normal || pbrt.type == pbrt_type::color) {
    val = pbrt.value3f;
    return true;
  } else if (pbrt.type == pbrt_type::real) {
    val = vec3f{pbrt.value1f, pbrt.value1f, pbrt.value1f};
    return true;
  } else {
    return false;
  }
}
inline bool get_pbrt_value(const pbrt_value& pbrt, vector<float>& val) {
  if (pbrt.type == pbrt_type::real) {
    if (!pbrt.vector1f.empty()) {
      val = pbrt.vector1f;
    } else {
      val = {pbrt.value1f};
    }
    return true;
  } else {
    return false;
  }
}
inline bool get_pbrt_value(const pbrt_value& pbrt, vector<vec2f>& val) {
  if (pbrt.type == pbrt_type::point2 || pbrt.type == pbrt_type::vector2) {
    if (!pbrt.vector2f.empty()) {
      val = pbrt.vector2f;
    } else {
      val = {pbrt.value2f};
    }
    return true;
  } else if (pbrt.type == pbrt_type::real) {
    if (pbrt.vector1f.empty() || (pbrt.vector1f.size() % 2) != 0)
      throw std::runtime_error("bad pbrt type");
    val.resize(pbrt.vector1f.size() / 2);
    for (auto i = 0; i < val.size(); i++)
      val[i] = {pbrt.vector1f[i * 2 + 0], pbrt.vector1f[i * 2 + 1]};
    return true;
  } else {
    return false;
  }
}
inline bool get_pbrt_value(const pbrt_value& pbrt, vector<vec3f>& val) {
  if (pbrt.type == pbrt_type::point || pbrt.type == pbrt_type::vector ||
      pbrt.type == pbrt_type::normal || pbrt.type == pbrt_type::color) {
    if (!pbrt.vector3f.empty()) {
      val = pbrt.vector3f;
    } else {
      val = {pbrt.value3f};
    }
    return true;
  } else if (pbrt.type == pbrt_type::real) {
    if (pbrt.vector1f.empty() || (pbrt.vector1f.size() % 3) != 0)
      throw std::invalid_argument{"expected float3 array"};
    val.resize(pbrt.vector1f.size() / 3);
    for (auto i = 0; i < val.size(); i++)
      val[i] = {pbrt.vector1f[i * 3 + 0], pbrt.vector1f[i * 3 + 1],
          pbrt.vector1f[i * 3 + 2]};
    return true;
  } else {
    return false;
  }
}

inline bool get_pbrt_value(const pbrt_value& pbrt, vector<int>& val) {
  if (pbrt.type == pbrt_type::integer) {
    if (!pbrt.vector1i.empty()) {
      val = pbrt.vector1i;
    } else {
      val = {pbrt.vector1i};
    }
    return true;
  } else {
    return false;
  }
}
inline bool get_pbrt_value(const pbrt_value& pbrt, vector<vec3i>& val) {
  if (pbrt.type == pbrt_type::integer) {
    if (pbrt.vector1i.empty() || (pbrt.vector1i.size() % 3) != 0)
      throw std::invalid_argument{"expected int3 array"};
    val.resize(pbrt.vector1i.size() / 3);
    for (auto i = 0; i < val.size(); i++)
      val[i] = {pbrt.vector1i[i * 3 + 0], pbrt.vector1i[i * 3 + 1],
          pbrt.vector1i[i * 3 + 2]};
    return true;
  } else {
    return false;
  }
}
inline bool get_pbrt_value(const pbrt_value& pbrt, pair<float, string>& val) {
  if (pbrt.type == pbrt_type::string || pbrt.type == pbrt_type::texture) {
    val.first = 0;
    return get_pbrt_value(pbrt, val.second);
  } else {
    val.second = "";
    return get_pbrt_value(pbrt, val.first);
  }
}
inline bool get_pbrt_value(const pbrt_value& pbrt, pair<vec3f, string>& val) {
  if (pbrt.type == pbrt_type::string || pbrt.type == pbrt_type::texture) {
    val.first = zero3f;
    return get_pbrt_value(pbrt, val.second);
  } else {
    val.second = "";
    return get_pbrt_value(pbrt, val.first);
  }
}
template <typename T>
inline bool get_pbrt_value(
    const vector<pbrt_value>& pbrt, const string& name, T& val) {
  for (auto& p : pbrt) {
    if (p.name == name) {
      return get_pbrt_value(p, val);
    }
  }
  return true;
}

// pbrt value construction
inline pbrt_value make_pbrt_value(
    const string& name, const string& val, pbrt_type type = pbrt_type::string) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1s = val;
  return pbrt;
}
inline pbrt_value make_pbrt_value(
    const string& name, bool val, pbrt_type type = pbrt_type::boolean) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1b = val;
  return pbrt;
}
inline pbrt_value make_pbrt_value(
    const string& name, int val, pbrt_type type = pbrt_type::integer) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1i = val;
  return pbrt;
}
inline pbrt_value make_pbrt_value(
    const string& name, float val, pbrt_type type = pbrt_type::real) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1f = val;
  return pbrt;
}
inline pbrt_value make_pbrt_value(
    const string& name, const vec2f& val, pbrt_type type = pbrt_type::point2) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value2f = val;
  return pbrt;
}
inline pbrt_value make_pbrt_value(
    const string& name, const vec3f& val, pbrt_type type = pbrt_type::color) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value3f = val;
  return pbrt;
}
inline pbrt_value make_pbrt_value(const string& name, const vector<vec2f>& val,
    pbrt_type type = pbrt_type::point2) {
  auto pbrt     = pbrt_value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector2f = val;
  return pbrt;
}
inline pbrt_value make_pbrt_value(const string& name, const vector<vec3f>& val,
    pbrt_type type = pbrt_type::point) {
  auto pbrt     = pbrt_value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector3f = val;
  return pbrt;
}
inline pbrt_value make_pbrt_value(const string& name, const vector<vec3i>& val,
    pbrt_type type = pbrt_type::integer) {
  auto pbrt     = pbrt_value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector1i = {&val.front().x, &val.front().x + val.size() * 3};
  return pbrt;
}

inline void remove_pbrt_comment(string_view& str, char comment_char = '#') {
  while (!str.empty() && is_newline(str.back())) str.remove_suffix(1);
  auto cpy       = str;
  auto in_string = false;
  while (!cpy.empty()) {
    if (cpy.front() == '"') in_string = !in_string;
    if (cpy.front() == comment_char && !in_string) break;
    cpy.remove_prefix(1);
  }
  str.remove_suffix(cpy.size());
}

// Read a pbrt command from file
inline bool read_pbrt_cmdline(file_stream& fs, string& cmd) {
  auto buffer = array<char, 4096>{};
  cmd.clear();
  auto found = false;
  auto pos   = ftell(fs.fs);
  while (read_line(fs, buffer)) {
    // line
    auto line = string_view{buffer.data()};
    remove_comment(line, '#', true);
    skip_whitespace(line);
    if (line.empty()) continue;

    // check if command
    auto is_cmd = line[0] >= 'A' && line[0] <= 'Z';
    if (is_cmd) {
      if (found) {
        fseek(fs.fs, pos, SEEK_SET);
        // line_num -= 1;
        return true;
      } else {
        found = true;
      }
    } else if (!found) {
      return false;
    }
    cmd += line;
    cmd += " ";
    pos = ftell(fs.fs);
  }
  return found;
}

// parse a quoted string
inline void parse_command(string_view& str, string& value) {
  skip_whitespace(str);
  if (!isalpha((int)str.front()))
    throw std::invalid_argument{"string expected"};
  auto pos = str.find_first_not_of(
      "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
  if (pos == string_view::npos) {
    value.assign(str);
    str.remove_prefix(str.size());
  } else {
    value.assign(str.substr(0, pos));
    str.remove_prefix(pos + 1);
  }
}

// parse pbrt value with optional parens
template <typename T>
inline void parse_param(string_view& str, T& value) {
  skip_whitespace(str);
  auto parens = !str.empty() && str.front() == '[';
  if (parens) str.remove_prefix(1);
  parse_value(str, value);
  if (!str.data()) throw std::invalid_argument{"value expected"};
  if (parens) {
    skip_whitespace(str);
    if (!str.empty() && str.front() == '[')
      throw std::invalid_argument{"array expected"};
    str.remove_prefix(1);
  }
}

// parse a quoted string
inline void parse_nametype(string_view& str_, string& name, string& type) {
  auto value = ""s;
  parse_value(str_, value);
  if (str_.empty()) throw std::invalid_argument{"string expected"};
  auto str  = string_view{value};
  auto pos1 = str.find(' ');
  if (pos1 == string_view::npos) throw std::invalid_argument{"string expected"};
  type = string(str.substr(0, pos1));
  str.remove_prefix(pos1);
  auto pos2 = str.find_first_not_of(' ');
  if (pos2 == string_view::npos) throw std::invalid_argument{"string expected"};
  str.remove_prefix(pos2);
  name = string(str);
}

inline pair<vec3f, vec3f> get_etak(const string& name) {
  static const unordered_map<string, pair<vec3f, vec3f>> metal_ior_table = {
      {"a-C", {{2.9440999183f, 2.2271502925f, 1.9681668794f},
                  {0.8874329109f, 0.7993216383f, 0.8152862927f}}},
      {"Ag", {{0.1552646489f, 0.1167232965f, 0.1383806959f},
                 {4.8283433224f, 3.1222459278f, 2.1469504455f}}},
      {"Al", {{1.6574599595f, 0.8803689579f, 0.5212287346f},
                 {9.2238691996f, 6.2695232477f, 4.8370012281f}}},
      {"AlAs", {{3.6051023902f, 3.2329365777f, 2.2175611545f},
                   {0.0006670247f, -0.0004999400f, 0.0074261204f}}},
      {"AlSb", {{-0.0485225705f, 4.1427547893f, 4.6697691348f},
                   {-0.0363741915f, 0.0937665154f, 1.3007390124f}}},
      {"Au", {{0.1431189557f, 0.3749570432f, 1.4424785571f},
                 {3.9831604247f, 2.3857207478f, 1.6032152899f}}},
      {"Be", {{4.1850592788f, 3.1850604423f, 2.7840913457f},
                 {3.8354398268f, 3.0101260162f, 2.8690088743f}}},
      {"Cr", {{4.3696828663f, 2.9167024892f, 1.6547005413f},
                 {5.2064337956f, 4.2313645277f, 3.7549467933f}}},
      {"CsI", {{2.1449030413f, 1.7023164587f, 1.6624194173f},
                  {0.0000000000f, 0.0000000000f, 0.0000000000f}}},
      {"Cu", {{0.2004376970f, 0.9240334304f, 1.1022119527f},
                 {3.9129485033f, 2.4528477015f, 2.1421879552f}}},
      {"Cu2O", {{3.5492833755f, 2.9520622449f, 2.7369202137f},
                   {0.1132179294f, 0.1946659670f, 0.6001681264f}}},
      {"CuO", {{3.2453822204f, 2.4496293965f, 2.1974114493f},
                  {0.5202739621f, 0.5707372756f, 0.7172250613f}}},
      {"d-C", {{2.7112524747f, 2.3185812849f, 2.2288565009f},
                  {0.0000000000f, 0.0000000000f, 0.0000000000f}}},
      {"Hg", {{2.3989314904f, 1.4400254917f, 0.9095512090f},
                 {6.3276269444f, 4.3719414152f, 3.4217899270f}}},
      {"HgTe", {{4.7795267752f, 3.2309984581f, 2.6600252401f},
                   {1.6319827058f, 1.5808189339f, 1.7295753852f}}},
      {"Ir", {{3.0864098394f, 2.0821938440f, 1.6178866805f},
                 {5.5921510077f, 4.0671757150f, 3.2672611269f}}},
      {"K", {{0.0640493070f, 0.0464100621f, 0.0381842017f},
                {2.1042155920f, 1.3489364357f, 0.9132113889f}}},
      {"Li", {{0.2657871942f, 0.1956102432f, 0.2209198538f},
                 {3.5401743407f, 2.3111306542f, 1.6685930000f}}},
      {"MgO", {{2.0895885542f, 1.6507224525f, 1.5948759692f},
                  {0.0000000000f, -0.0000000000f, 0.0000000000f}}},
      {"Mo", {{4.4837010280f, 3.5254578255f, 2.7760769438f},
                 {4.1111307988f, 3.4208716252f, 3.1506031404f}}},
      {"Na", {{0.0602665320f, 0.0561412435f, 0.0619909494f},
                 {3.1792906496f, 2.1124800781f, 1.5790940266f}}},
      {"Nb", {{3.4201353595f, 2.7901921379f, 2.3955856658f},
                 {3.4413817900f, 2.7376437930f, 2.5799132708f}}},
      {"Ni", {{2.3672753521f, 1.6633583302f, 1.4670554172f},
                 {4.4988329911f, 3.0501643957f, 2.3454274399f}}},
      {"Rh", {{2.5857954933f, 1.8601866068f, 1.5544279524f},
                 {6.7822927110f, 4.7029501026f, 3.9760892461f}}},
      {"Se-e", {{5.7242724833f, 4.1653992967f, 4.0816099264f},
                   {0.8713747439f, 1.1052845009f, 1.5647788766f}}},
      {"Se", {{4.0592611085f, 2.8426947380f, 2.8207582835f},
                 {0.7543791750f, 0.6385150558f, 0.5215872029f}}},
      {"SiC", {{3.1723450205f, 2.5259677964f, 2.4793623897f},
                  {0.0000007284f, -0.0000006859f, 0.0000100150f}}},
      {"SnTe", {{4.5251865890f, 1.9811525984f, 1.2816819226f},
                   {0.0000000000f, 0.0000000000f, 0.0000000000f}}},
      {"Ta", {{2.0625846607f, 2.3930915569f, 2.6280684948f},
                 {2.4080467973f, 1.7413705864f, 1.9470377016f}}},
      {"Te-e", {{7.5090397678f, 4.2964603080f, 2.3698732430f},
                   {5.5842076830f, 4.9476231084f, 3.9975145063f}}},
      {"Te", {{7.3908396088f, 4.4821028985f, 2.6370708478f},
                 {3.2561412892f, 3.5273908133f, 3.2921683116f}}},
      {"ThF4", {{1.8307187117f, 1.4422274283f, 1.3876488528f},
                   {0.0000000000f, 0.0000000000f, 0.0000000000f}}},
      {"TiC", {{3.7004673762f, 2.8374356509f, 2.5823030278f},
                  {3.2656905818f, 2.3515586388f, 2.1727857800f}}},
      {"TiN", {{1.6484691607f, 1.1504482522f, 1.3797795097f},
                  {3.3684596226f, 1.9434888540f, 1.1020123347f}}},
      {"TiO2-e", {{3.1065574823f, 2.5131551146f, 2.5823844157f},
                     {0.0000289537f, -0.0000251484f, 0.0001775555f}}},
      {"TiO2", {{3.4566203131f, 2.8017076558f, 2.9051485020f},
                   {0.0001026662f, -0.0000897534f, 0.0006356902f}}},
      {"VC", {{3.6575665991f, 2.7527298065f, 2.5326814570f},
                 {3.0683516659f, 2.1986687713f, 1.9631816252f}}},
      {"VN", {{2.8656011588f, 2.1191817791f, 1.9400767149f},
                 {3.0323264950f, 2.0561075580f, 1.6162930914f}}},
      {"V", {{4.2775126218f, 3.5131538236f, 2.7611257461f},
                {3.4911844504f, 2.8893580874f, 3.1116965117f}}},
      {"W", {{4.3707029924f, 3.3002972445f, 2.9982666528f},
                {3.5006778591f, 2.6048652781f, 2.2731930614f}}},
  };
  return metal_ior_table.at(name);
}

// Pbrt measure subsurface parameters (sigma_prime_s, sigma_a in mm^-1)
// from pbrt code at pbrt/code/medium.cpp
inline pair<vec3f, vec3f> get_subsurface(const string& name) {
  static const unordered_map<string, pair<vec3f, vec3f>> params = {
      // From "A Practical Model for Subsurface Light Transport"
      // Jensen, Marschner, Levoy, Hanrahan
      // Proc SIGGRAPH 2001
      {"Apple", {{2.29f, 2.39f, 1.97f}, {0.0030f, 0.0034f, 0.046f}}},
      {"Chicken1", {{0.15f, 0.21f, 0.38f}, {0.015f, 0.077f, 0.19f}}},
      {"Chicken2", {{0.19f, 0.25f, 0.32f}, {0.018f, 0.088f, 0.20f}}},
      {"Cream", {{7.38f, 5.47f, 3.15f}, {0.0002f, 0.0028f, 0.0163f}}},
      {"Ketchup", {{0.18f, 0.07f, 0.03f}, {0.061f, 0.97f, 1.45f}}},
      {"Marble", {{2.19f, 2.62f, 3.00f}, {0.0021f, 0.0041f, 0.0071f}}},
      {"Potato", {{0.68f, 0.70f, 0.55f}, {0.0024f, 0.0090f, 0.12f}}},
      {"Skimmilk", {{0.70f, 1.22f, 1.90f}, {0.0014f, 0.0025f, 0.0142f}}},
      {"Skin1", {{0.74f, 0.88f, 1.01f}, {0.032f, 0.17f, 0.48f}}},
      {"Skin2", {{1.09f, 1.59f, 1.79f}, {0.013f, 0.070f, 0.145f}}},
      {"Spectralon", {{11.6f, 20.4f, 14.9f}, {0.00f, 0.00f, 0.00f}}},
      {"Wholemilk", {{2.55f, 3.21f, 3.77f}, {0.0011f, 0.0024f, 0.014f}}},
      // From "Acquiring Scattering Properties of Participating Media by
      // Dilution",
      // Narasimhan, Gupta, Donner, Ramamoorthi, Nayar, Jensen
      // Proc SIGGRAPH 2006
      {"Lowfat Milk",
          {{0.89187f, 1.5136f, 2.532f}, {0.002875f, 0.00575f, 0.0115f}}},
      {"Reduced Milk",
          {{2.4858f, 3.1669f, 4.5214f}, {0.0025556f, 0.0051111f, 0.012778f}}},
      {"Regular Milk",
          {{4.5513f, 5.8294f, 7.136f}, {0.0015333f, 0.0046f, 0.019933f}}},
      {"Espresso",
          {{0.72378f, 0.84557f, 1.0247f}, {4.7984f, 6.5751f, 8.8493f}}},
      {"Mint Mocha Coffee",
          {{0.31602f, 0.38538f, 0.48131f}, {3.772f, 5.8228f, 7.82f}}},
      {"Lowfat Soy Milk", {{0.30576f, 0.34233f, 0.61664f},
                              {0.0014375f, 0.0071875f, 0.035937f}}},
      {"Regular Soy Milk",
          {{0.59223f, 0.73866f, 1.4693f}, {0.0019167f, 0.0095833f, 0.065167f}}},
      {"Lowfat Chocolate Milk",
          {{0.64925f, 0.83916f, 1.1057f}, {0.0115f, 0.0368f, 0.1564f}}},
      {"Regular Chocolate Milk",
          {{1.4585f, 2.1289f, 2.9527f}, {0.010063f, 0.043125f, 0.14375f}}},
      {"Coke",
          {{8.9053e-05f, 8.372e-05f, 0.0f}, {0.10014f, 0.16503f, 0.2468f}}},
      {"Pepsi",
          {{6.1697e-05f, 4.2564e-05f, 0.0f}, {0.091641f, 0.14158f, 0.20729f}}},
      {"Sprite", {{6.0306e-06f, 6.4139e-06f, 6.5504e-06f},
                     {0.001886f, 0.0018308f, 0.0020025f}}},
      {"Gatorade", {{0.0024574f, 0.003007f, 0.0037325f},
                       {0.024794f, 0.019289f, 0.008878f}}},
      {"Chardonnay", {{1.7982e-05f, 1.3758e-05f, 1.2023e-05f},
                         {0.010782f, 0.011855f, 0.023997f}}},
      {"White Zinfandel", {{1.7501e-05f, 1.9069e-05f, 1.288e-05f},
                              {0.012072f, 0.016184f, 0.019843f}}},
      {"Merlot", {{2.1129e-05f, 0.0f, 0.0f}, {0.11632f, 0.25191f, 0.29434f}}},
      {"Budweiser Beer", {{2.4356e-05f, 2.4079e-05f, 1.0564e-05f},
                             {0.011492f, 0.024911f, 0.057786f}}},
      {"Coors Light Beer",
          {{5.0922e-05f, 4.301e-05f, 0.0f}, {0.006164f, 0.013984f, 0.034983f}}},
      {"Clorox", {{0.0024035f, 0.0031373f, 0.003991f},
                     {0.0033542f, 0.014892f, 0.026297f}}},
      {"Apple Juice", {{0.00013612f, 0.00015836f, 0.000227f},
                          {0.012957f, 0.023741f, 0.052184f}}},
      {"Cranberry Juice", {{0.00010402f, 0.00011646f, 7.8139e-05f},
                              {0.039437f, 0.094223f, 0.12426f}}},
      {"Grape Juice",
          {{5.382e-05f, 0.0f, 0.0f}, {0.10404f, 0.23958f, 0.29325f}}},
      {"Ruby Grapefruit Juice",
          {{0.011002f, 0.010927f, 0.011036f}, {0.085867f, 0.18314f, 0.25262f}}},
      {"White Grapefruit Juice",
          {{0.22826f, 0.23998f, 0.32748f}, {0.0138f, 0.018831f, 0.056781f}}},
      {"Shampoo", {{0.0007176f, 0.0008303f, 0.0009016f},
                      {0.014107f, 0.045693f, 0.061717f}}},
      {"Strawberry Shampoo", {{0.00015671f, 0.00015947f, 1.518e-05f},
                                 {0.01449f, 0.05796f, 0.075823f}}},
      {"Head & Shoulders Shampoo",
          {{0.023805f, 0.028804f, 0.034306f}, {0.084621f, 0.15688f, 0.20365f}}},
      {"Lemon Tea Powder",
          {{0.040224f, 0.045264f, 0.051081f}, {2.4288f, 4.5757f, 7.2127f}}},
      {"Orange Powder", {{0.00015617f, 0.00017482f, 0.0001762f},
                            {0.001449f, 0.003441f, 0.007863f}}},
      {"Pink Lemonade Powder", {{0.00012103f, 0.00013073f, 0.00012528f},
                                   {0.001165f, 0.002366f, 0.003195f}}},
      {"Cappuccino Powder",
          {{1.8436f, 2.5851f, 2.1662f}, {35.844f, 49.547f, 61.084f}}},
      {"Salt Powder",
          {{0.027333f, 0.032451f, 0.031979f}, {0.28415f, 0.3257f, 0.34148f}}},
      {"Sugar Powder", {{0.00022272f, 0.00025513f, 0.000271f},
                           {0.012638f, 0.031051f, 0.050124f}}},
      {"Suisse Mocha Powder",
          {{2.7979f, 3.5452f, 4.3365f}, {17.502f, 27.004f, 35.433f}}},
      {"Pacific Ocean Surface Water", {{0.0001764f, 0.00032095f, 0.00019617f},
                                          {0.031845f, 0.031324f, 0.030147f}}},
  };
  return params.at(name);
}

inline void parse_params(string_view& str, vector<pbrt_value>& values) {
  auto parse_pvalues = [](string_view& str, auto& value, auto& values) {
    values.clear();
    skip_whitespace(str);
    if (str.empty()) throw std::invalid_argument{"param expected"};
    if (str.front() == '[') {
      str.remove_prefix(1);
      skip_whitespace(str);
      if (str.empty()) throw std::invalid_argument{"param expected"};
      while (!str.empty()) {
        auto& val = values.empty() ? value : values.emplace_back();
        parse_value(str, val);
        if (!str.data()) throw std::invalid_argument{"param expected"};
        skip_whitespace(str);
        if (str.empty()) break;
        if (str.front() == ']') break;
        if (values.empty()) values.push_back(value);
      }
      if (str.empty()) throw std::invalid_argument{"param expected"};
      if (str.front() != ']') throw std::invalid_argument{"param expected"};
      str.remove_prefix(1);
    } else {
      return parse_value(str, value);
    }
  };

  auto starts_with = [](string_view value, string_view prefix) {
    if (prefix.size() > value.size()) return false;
    return value.rfind(prefix, 0) == 0;
  };
  auto ends_with = [](string_view value, string_view postfix) {
    if (postfix.size() > value.size()) return false;
    return std::equal(postfix.rbegin(), postfix.rend(), value.rbegin());
  };

  values.clear();
  skip_whitespace(str);
  while (!str.empty()) {
    auto& value = values.emplace_back();
    auto  type  = ""s;
    parse_nametype(str, value.name, type);
    skip_whitespace(str);
    if (str.empty()) throw std::invalid_argument{"param expected"};
    if (type == "float") {
      value.type = pbrt_type::real;
      parse_pvalues(str, value.value1f, value.vector1f);
    } else if (type == "integer") {
      value.type = pbrt_type::integer;
      parse_pvalues(str, value.value1i, value.vector1i);
    } else if (type == "string") {
      auto vector1s = vector<string>{};
      value.type    = pbrt_type::string;
      parse_pvalues(str, value.value1s, vector1s);
      if (!vector1s.empty()) throw std::invalid_argument{"param expected"};
    } else if (type == "bool") {
      auto value1s  = ""s;
      auto vector1s = vector<string>{};
      value.type    = pbrt_type::boolean;
      parse_pvalues(str, value1s, vector1s);
      if (!vector1s.empty()) throw std::invalid_argument{"param expected"};
      value.value1b = value1s == "true";
    } else if (type == "texture") {
      auto vector1s = vector<string>{};
      value.type    = pbrt_type::texture;
      parse_pvalues(str, value.value1s, vector1s);
      if (!vector1s.empty()) throw std::invalid_argument{"param expected"};
    } else if (type == "point" || type == "point3") {
      value.type = pbrt_type::point;
      parse_pvalues(str, value.value3f, value.vector3f);
    } else if (type == "normal" || type == "normal3") {
      value.type = pbrt_type::normal;
      parse_pvalues(str, value.value3f, value.vector3f);
    } else if (type == "vector" || type == "vector3") {
      value.type = pbrt_type::vector;
      parse_pvalues(str, value.value3f, value.vector3f);
    } else if (type == "point2") {
      value.type = pbrt_type::point2;
      parse_pvalues(str, value.value2f, value.vector2f);
    } else if (type == "vector2") {
      value.type = pbrt_type::vector2;
      parse_pvalues(str, value.value2f, value.vector2f);
    } else if (type == "blackbody") {
      value.type = pbrt_type::color;
      // auto blackbody = zero2f;
      // auto vec tor2f  = vector<vec2f>{};
      // parse_pvalues(str, blackbody, vector2f);
      // if (!vector2f.empty()) return false;
      // value.value3f = blackbody_to_rgb(blackbody.x) * blackbody.y;
      auto blackbody = 0.0f;
      auto vector1f  = vector<float>{};
      parse_pvalues(str, blackbody, vector1f);
      if (vector1f.size() < 2) {
        value.value3f = blackbody_to_rgb(blackbody);
      } else {
        value.value3f = blackbody_to_rgb(vector1f[0]) * vector1f[1];
      }
    } else if (type == "color" || type == "rgb") {
      value.type = pbrt_type::color;
      parse_pvalues(str, value.value3f, value.vector3f);
    } else if (type == "xyz") {
      value.type = pbrt_type::color;
      parse_pvalues(str, value.value3f, value.vector3f);
      // xyz conversion
      throw std::invalid_argument{"xyz not supported"};
    } else if (type == "spectrum") {
      auto is_string = false;
      auto str1      = str;
      skip_whitespace(str1);
      if (!str1.empty() && str1.front() == '"') {
        is_string = true;
      } else if (!str1.empty() && str1.front() == '[') {
        str1.remove_prefix(1);
        skip_whitespace(str1);
        if (!str1.empty() && str1.front() == '"') is_string = true;
      }
      if (is_string) {
        value.type     = pbrt_type::color;
        auto filename  = ""s;
        auto filenames = vector<string>{};
        skip_whitespace(str);
        auto has_parens = str.front() == '[';
        if (has_parens) str.remove_prefix(1);
        parse_value(str, filename);
        if (has_parens) {
          skip_whitespace(str);
          if (str.front() != ']') throw std::invalid_argument{"param expected"};
          str.remove_prefix(1);
        }
        if (str.empty()) throw std::invalid_argument{"param expected"};
        auto filenamep = path_filename(filename);
        auto name      = string_view{filenamep};
        if (ends_with(name, ".spd")) {
          name.remove_suffix(4);
          if (name == "SHPS") {
            value.value3f = {1, 1, 1};
          } else if (ends_with(name, ".eta")) {
            name.remove_suffix(4);
            auto eta      = get_etak(string{name}).first;
            value.value3f = {eta.x, eta.y, eta.z};
          } else if (ends_with(name, ".k")) {
            name.remove_suffix(2);
            auto k        = get_etak(string{name}).second;
            value.value3f = {k.x, k.y, k.z};
          } else {
            throw std::invalid_argument{"param expected"};
          }
        } else if (starts_with(name, "metal-")) {
          name.remove_prefix(6);
          if (ends_with(name, "-eta")) {
            name.remove_suffix(4);
            auto eta      = get_etak(string{name}).first;
            value.value3f = {eta.x, eta.y, eta.z};
          } else if (ends_with(name, "-k")) {
            name.remove_suffix(2);
            auto k        = get_etak(string{name}).second;
            value.value3f = {k.x, k.y, k.z};
          } else {
            throw std::invalid_argument{"param expected"};
          }
        } else if (starts_with(name, "glass-")) {
          value.value3f = {1.5, 1.5, 1.5};
        } else {
          throw std::invalid_argument{"param expected"};
        }
      } else {
        value.type = pbrt_type::spectrum;
        parse_pvalues(str, value.value1f, value.vector1f);
      }
    } else {
      throw std::invalid_argument{"param expected"};
    }
    skip_whitespace(str);
  }
}

// Other pbrt elements
struct pbrt_film {
  // film approximation
  string filename   = "";
  vec2i  resolution = {0, 0};
};

// Pbrt area light
struct pbrt_arealight {
  // arealight parameters
  string name     = "";
  vec3f  emission = {0, 0, 0};
};

// Pbrt medium. Not parsed at the moment.
struct pbrt_medium {
  // medium parameters
  string name = "";
};

// convert pbrt films
inline void convert_film(pbrt_film& film, const pbrt_command& command,
    const string& filename, bool verbose = false) {
  try {
    if (command.type == "image") {
      film.resolution = {512, 512};
      get_pbrt_value(command.values, "xresolution", film.resolution.x);
      get_pbrt_value(command.values, "yresolution", film.resolution.y);
      film.filename = "out.png"s;
      get_pbrt_value(command.values, "filename", film.filename);
    } else if (command.type == "rgb") {
      film.resolution = {512, 512};
      get_pbrt_value(command.values, "xresolution", film.resolution.x);
      get_pbrt_value(command.values, "yresolution", film.resolution.y);
      film.filename = "out.png"s;
      get_pbrt_value(command.values, "filename", film.filename);
    } else {
      throw io_error::type_error(filename, command.type);
    }
  } catch (const io_error& error) {
    throw;
  } catch (...) {
    throw io_error::parse_error(filename);
  }
}

// convert pbrt elements
inline void convert_camera(pbrt_camera& pcamera, const pbrt_command& command,
    const vec2i& resolution, const string& filename, bool verbose = false) {
  pcamera.frame      = command.frame;
  pcamera.frend      = command.frend;
  pcamera.frame      = inverse((frame3f)pcamera.frame);
  pcamera.frame.z    = -pcamera.frame.z;
  pcamera.resolution = resolution;
  auto film_aspect =
      (resolution == zero2i) ? 1 : (float)resolution.x / (float)resolution.y;
  try {
    if (command.type == "perspective") {
      auto fov = 90.0f;
      get_pbrt_value(command.values, "fov", fov);
      // auto lensradius = if(!get_pbrt_value(values, "lensradius", 0.0f);
      pcamera.aspect = film_aspect;
      if (pcamera.aspect >= 1) {
        pcamera.lens = (0.036f / pcamera.aspect) / (2 * tan(radians(fov) / 2));
      } else {
        pcamera.lens = (0.036f * pcamera.aspect) / (2 * tan(radians(fov) / 2));
      }
      get_pbrt_value(command.values, "frameaspectratio", pcamera.aspect);
      pcamera.focus = 10.0f;
      get_pbrt_value(command.values, "focaldistance", pcamera.focus);
    } else if (command.type == "realistic") {
      auto lensfile = ""s;
      get_pbrt_value(command.values, "lensfile", lensfile);
      lensfile     = lensfile.substr(0, lensfile.size() - 4);
      lensfile     = lensfile.substr(lensfile.find('.') + 1);
      lensfile     = lensfile.substr(0, lensfile.size() - 2);
      auto lens    = max((float)std::atof(lensfile.c_str()), 35.0f) * 0.001f;
      pcamera.lens = 2 * atan(0.036f / (2 * lens));
      pcamera.aperture = 0.0f;
      get_pbrt_value(command.values, "aperturediameter", pcamera.aperture);
      pcamera.focus = 10.0f;
      get_pbrt_value(command.values, "focusdistance", pcamera.focus);
      pcamera.aspect = film_aspect;
    } else {
      throw io_error::type_error(filename, command.type);
    }
  } catch (const io_error& error) {
    throw;
  } catch (...) {
    throw io_error::parse_error(filename);
  }
}

// convert pbrt textures
inline void convert_texture(pbrt_texture& ptexture, const pbrt_command& command,
    unordered_map<string, pbrt_texture>& texture_map, const string& filename,
    bool verbose = false) {
  auto make_filename = [&texture_map](const string& name) {
    if (name.empty()) return ""s;
    auto pos = texture_map.find(name);
    if (pos == texture_map.end()) return ""s;
    return pos->second.filename;
  };

  try {
    ptexture.name = command.name;
    if (command.type == "imagemap") {
      ptexture.filename = "";
      get_pbrt_value(command.values, "filename", ptexture.filename);
    } else if (command.type == "constant") {
      ptexture.constant = vec3f{1, 1, 1};
      get_pbrt_value(command.values, "value", ptexture.constant);
    } else if (command.type == "bilerp") {
      ptexture.constant = {1, 0, 0};
    } else if (command.type == "checkerboard") {
      // auto tex1     = if(!get_pbrt_value(command.values, "tex1",
      // pair{vec3f{1,1,1},
      // ""s}); auto tex2     = if(!get_pbrt_value(command.values, "tex2",
      //  pair{vec3f{0}, ""s}); auto rgb1     = tex1.second == "" ?
      //  tex1.first :
      // vec3f{0.4f, 0.4f, 0.4f}; auto rgb2     = tex1.second == "" ? tex2.first
      // : vec3f{0.6f, 0.6f, 0.6f}; auto params   = proc_image_params{};
      // params.type = proc_image_params::type_t::checker; params.color0 =
      // {rgb1.x, rgb1.y, rgb1.z, 1}; params.color1 = {rgb2.x, rgb2.y, rgb2.z,
      // 1}; params.scale = 2; make_proc_image(texture.hdr, params);
      // float_to_byte(texture.ldr, texture.hdr); texture.hdr = {};
      ptexture.constant = {0.5, 0.5, 0.5};
    } else if (command.type == "dots") {
      ptexture.constant = {0.5, 0.5, 0.5};
    } else if (command.type == "fbm") {
      ptexture.constant = {0.5, 0.5, 0.5};
    } else if (command.type == "marble") {
      ptexture.constant = {0.5, 0.5, 0.5};
    } else if (command.type == "mix") {
      auto tex1 = pair{vec3f{0, 0, 0}, ""s}, tex2 = pair{vec3f{1, 1, 1}, ""s};
      get_pbrt_value(command.values, "tex1", tex1);
      get_pbrt_value(command.values, "tex2", tex2);
      if (!make_filename(tex1.second).empty()) {
        ptexture.filename = make_filename(tex1.second);
      } else if (!make_filename(tex2.second).empty()) {
        ptexture.filename = make_filename(tex2.second);
      } else {
        ptexture.constant = {1, 0, 0};
      }
    } else if (command.type == "scale") {
      auto tex1 = pair{vec3f{1, 1, 1}, ""s}, tex2 = pair{vec3f{1, 1, 1}, ""s};
      get_pbrt_value(command.values, "tex1", tex2);
      get_pbrt_value(command.values, "tex2", tex1);
      if (!make_filename(tex1.second).empty()) {
        ptexture.filename = make_filename(tex1.second);
      } else if (!make_filename(tex2.second).empty()) {
        ptexture.filename = make_filename(tex2.second);
      } else {
        ptexture.constant = {1, 0, 0};
      }
    } else if (command.type == "uv") {
      ptexture.constant = {1, 0, 0};
    } else if (command.type == "windy") {
      ptexture.constant = {1, 0, 0};
    } else if (command.type == "wrinkled") {
      ptexture.constant = {1, 0, 0};
    } else {
      throw io_error::type_error(filename, command.type);
    }
  } catch (const io_error& error) {
    throw;
  } catch (...) {
    throw io_error::parse_error(filename);
  }
}

// convert pbrt materials
inline void convert_material(pbrt_material& pmaterial,
    const pbrt_command& command, unordered_map<string, int>& texture_map,
    const unordered_map<string, pbrt_material>& named_materials,
    const unordered_map<string, pbrt_texture>&  named_textures,
    const string& filename, bool verbose = false) {
  // helpers
  auto get_texture_id = [&texture_map](const string& path) {
    if (path.empty()) return -1;
    auto texture_it = texture_map.find(path);
    if (texture_it == texture_map.end()) {
      auto texture_id   = (int)texture_map.size();
      texture_map[path] = texture_id;
      return texture_id;
    } else {
      return texture_it->second;
    }
  };
  auto get_texture = [&](const vector<pbrt_value>& values, const string& name,
                         vec3f& color, int& texture_id, const vec3f& def) {
    auto textured = pair{def, ""s};
    get_pbrt_value(values, name, textured);
    if (textured.second.empty()) {
      color      = textured.first;
      texture_id = -1;
    } else {
      auto& texture = named_textures.at(textured.second);
      if (texture.filename.empty()) {
        color      = texture.constant;
        texture_id = -1;
      } else {
        color      = {1, 1, 1};
        texture_id = get_texture_id(texture.filename);
      }
    }
  };
  auto get_scalar = [&](const vector<pbrt_value>& values, const string& name,
                        float& scalar, float def) {
    auto textured = pair{vec3f{def, def, def}, ""s};
    get_pbrt_value(values, name, textured);
    if (textured.second.empty()) {
      scalar = mean(textured.first);
    } else {
      auto& texture = named_textures.at(textured.second);
      if (texture.filename.empty()) {
        scalar = mean(texture.constant);
      } else {
        scalar = def;
      }
    }
  };
  auto get_color = [&](const vector<pbrt_value>& values, const string& name,
                       vec3f& color, const vec3f& def) {
    auto textured = pair{def, ""s};
    get_pbrt_value(values, name, textured);
    if (textured.second.empty()) {
      color = textured.first;
    } else {
      auto& texture = named_textures.at(textured.second);
      if (texture.filename.empty()) {
        color = texture.constant;
      } else {
        color = def;
      }
    }
  };

  auto get_roughness = [&](const vector<pbrt_value>& values, float& roughness,
                           float def = 0.1) {
    auto roughness_ = pair{vec3f{def, def, def}, ""s};
    get_pbrt_value(values, "roughness", roughness_);
    auto uroughness = roughness_, vroughness = roughness_;
    auto remaproughness = true;
    get_pbrt_value(values, "uroughness", uroughness);
    get_pbrt_value(values, "vroughness", vroughness);
    get_pbrt_value(values, "remaproughness", remaproughness);

    roughness = 0;
    if (uroughness.first == zero3f || vroughness.first == zero3f) return;
    roughness = mean(vec2f{mean(uroughness.first), mean(vroughness.first)});
    // from pbrt code
    if (remaproughness) {
      roughness = max(roughness, 1e-3f);
      auto x    = log(roughness);
      roughness = 1.62142f + 0.819955f * x + 0.1734f * x * x +
                  0.0171201f * x * x * x + 0.000640711f * x * x * x * x;
    }
    roughness = sqrt(roughness);
  };

  auto eta_to_reflectivity = [](const vec3f&  eta,
                                 const vec3f& etak = zero3f) -> vec3f {
    return ((eta - 1) * (eta - 1) + etak * etak) /
           ((eta + 1) * (eta + 1) + etak * etak);
  };

  try {
    pmaterial.name = command.name;
    if (command.type == "uber") {
      auto diffuse = zero3f, specular = zero3f, transmission = zero3f;
      auto diffuse_map = -1, specular_map = -1, transmission_map = -1;
      get_texture(
          command.values, "Kd", diffuse, diffuse_map, vec3f{0.25, 0.25, 0.25});
      get_texture(command.values, "Ks", specular, specular_map,
          vec3f{0.25, 0.25, 0.25});
      get_texture(
          command.values, "Kt", transmission, transmission_map, vec3f{0, 0, 0});
      if (max(transmission) > 0.1) {
        pmaterial.type      = pbrt_mtype::thinglass;
        pmaterial.color     = transmission;
        pmaterial.color_tex = transmission_map;
      } else if (max(specular) > 0.1) {
        pmaterial.type      = pbrt_mtype::plastic;
        pmaterial.color     = diffuse;
        pmaterial.color_tex = diffuse_map;
      } else {
        pmaterial.type      = pbrt_mtype::plastic;
        pmaterial.color     = diffuse;
        pmaterial.color_tex = diffuse_map;
      }
      get_scalar(command.values, "opacity", pmaterial.opacity, 1);
      get_scalar(command.values, "eta", pmaterial.ior, 1.5f);
      get_roughness(command.values, pmaterial.roughness, 0.1f);
    } else if (command.type == "plastic") {
      pmaterial.type = pbrt_mtype::plastic;
      get_texture(command.values, "Kd", pmaterial.color, pmaterial.color_tex,
          vec3f{0.25, 0.25, 0.25});
      // get_scalar(command.values, "Ks", pmaterial.specular, 0.25))
      //   return parse_error();
      get_scalar(command.values, "eta", pmaterial.ior, 1.5f);
      pmaterial.roughness = 0.1f;
      get_roughness(command.values, pmaterial.roughness, 0.1f);
    } else if (command.type == "coateddiffuse") {
      pmaterial.type = pbrt_mtype::plastic;
      get_texture(command.values, "reflectance", pmaterial.color,
          pmaterial.color_tex, vec3f{0.25, 0.25, 0.25});
      get_scalar(command.values, "eta", pmaterial.ior, 1.5f);
      pmaterial.roughness = 0.1f;
      get_roughness(command.values, pmaterial.roughness, 0.1f);
    } else if (command.type == "translucent") {
      // not well supported yet
      pmaterial.type = pbrt_mtype::matte;
      get_texture(command.values, "Kd", pmaterial.color, pmaterial.color_tex,
          vec3f{0.25, 0.25, 0.25});
      // get_scalar(command.values, "Ks", pmaterial.specular, 0.25))
      //   return parse_error();
      // get_scalar(command.values, "eta", pmaterial.ior, 1.5))
      //   return parse_error();
      // get_roughness(command.values, pmaterial.roughness, 0.1))
      //   return parse_error();
    } else if (command.type == "diffusetransmission") {
      // not well supported yet
      pmaterial.type = pbrt_mtype::matte;
      get_texture(command.values, "reflectance", pmaterial.color,
          pmaterial.color_tex, vec3f{0.25f, 0.25f, 0.25f});
      // get_texture(command.values, "transmittance", pmaterial.color,
      //         pmaterial.color_tex, vec3f{0.25, 0.25, 0.25}))
      //   return parse_error();
    } else if (command.type == "matte") {
      pmaterial.type = pbrt_mtype::matte;
      get_texture(command.values, "Kd", pmaterial.color, pmaterial.color_tex,
          vec3f{0.5, 0.5, 0.5});
    } else if (command.type == "diffuse") {
      pmaterial.type = pbrt_mtype::matte;
      get_texture(command.values, "reflectance", pmaterial.color,
          pmaterial.color_tex, vec3f{0.5f, 0.5f, 0.5f});
    } else if (command.type == "mirror") {
      pmaterial.type = pbrt_mtype::metal;
      get_texture(command.values, "Kr", pmaterial.color, pmaterial.color_tex,
          vec3f{0.9f, 0.9f, 0.9f});
      pmaterial.roughness = 0;
    } else if (command.type == "metal") {
      pmaterial.type = pbrt_mtype::metal;
      // get_texture(
      //     values, "Kr", material->specular, material->specular_tex,
      //     vec3f{1,1,1});
      auto eta = zero3f, etak = zero3f;
      get_color(command.values, "eta", eta,
          vec3f{0.2004376970f, 0.9240334304f, 1.1022119527f});
      get_color(command.values, "k", etak,
          vec3f{3.9129485033f, 2.4528477015f, 2.1421879552f});
      pmaterial.color     = eta_to_reflectivity(eta, etak);
      pmaterial.roughness = 0.01f;
      get_roughness(command.values, pmaterial.roughness, 0.01f);
    } else if (command.type == "conductor") {
      pmaterial.type = pbrt_mtype::metal;
      auto eta = zero3f, etak = zero3f;
      get_color(command.values, "eta", eta,
          vec3f{0.2004376970f, 0.9240334304f, 1.1022119527f});
      get_color(command.values, "k", etak,
          vec3f{3.9129485033f, 2.4528477015f, 2.1421879552f});
      pmaterial.color     = eta_to_reflectivity(eta, etak);
      pmaterial.roughness = 0.01f;
      get_roughness(command.values, pmaterial.roughness, 0.01f);
    } else if (command.type == "coatedconductor") {
      pmaterial.type = pbrt_mtype::metal;
      auto eta = zero3f, etak = zero3f;
      get_color(command.values, "conductor.eta", eta,
          vec3f{0.2004376970f, 0.9240334304f, 1.1022119527f});
      get_color(command.values, "conductor.k", etak,
          vec3f{3.9129485033f, 2.4528477015f, 2.1421879552f});
      pmaterial.color     = eta_to_reflectivity(eta, etak);
      pmaterial.roughness = 0.01f;
      get_roughness(command.values, pmaterial.roughness, 0.01f);
    } else if (command.type == "substrate") {
      // not well supported
      pmaterial.type = pbrt_mtype::plastic;
      get_texture(command.values, "Kd", pmaterial.color, pmaterial.color_tex,
          vec3f{0.5f, 0.5f, 0.5f});
      auto specular = 0.0f;
      get_scalar(command.values, "Ks", specular, 0.5f);
      get_scalar(command.values, "eta", pmaterial.ior, 1.5f);
      pmaterial.roughness = 0.1f;
      get_roughness(command.values, pmaterial.roughness, 0.1f);
    } else if (command.type == "glass") {
      pmaterial.type = pbrt_mtype::glass;
      get_texture(command.values, "Kt", pmaterial.color, pmaterial.color_tex,
          vec3f{1, 1, 1});
      get_scalar(command.values, "eta", pmaterial.ior, 1.5f);
      pmaterial.roughness = 0;
      get_roughness(command.values, pmaterial.roughness, 0.0f);
    } else if (command.type == "dielectric") {
      pmaterial.type  = pbrt_mtype::glass;
      pmaterial.color = {1, 1, 1};
      get_scalar(command.values, "eta", pmaterial.ior, 1.5f);
      pmaterial.roughness = 0;
      get_roughness(command.values, pmaterial.roughness, 0.0f);
    } else if (command.type == "thindielectric") {
      pmaterial.type  = pbrt_mtype::thinglass;
      pmaterial.color = {1, 1, 1};
      get_scalar(command.values, "eta", pmaterial.ior, 1.5f);
      pmaterial.roughness = 0;
      get_roughness(command.values, pmaterial.roughness, 0.0f);
    } else if (command.type == "hair") {
      pmaterial.type = pbrt_mtype::matte;
      get_texture(command.values, "color", pmaterial.color, pmaterial.color_tex,
          vec3f{0, 0, 0});
      pmaterial.roughness = 1;
      if (verbose) printf("hair material not properly supported\n");
    } else if (command.type == "disney") {
      pmaterial.type = pbrt_mtype::matte;
      get_texture(command.values, "color", pmaterial.color, pmaterial.color_tex,
          vec3f{0.5f, 0.5f, 0.5f});
      pmaterial.roughness = 1;
      if (verbose) printf("disney material not properly supported\n");
    } else if (command.type == "kdsubsurface") {
      pmaterial.type = pbrt_mtype::plastic;
      get_texture(command.values, "Kd", pmaterial.color, pmaterial.color_tex,
          vec3f{0.5f, 0.5f, 0.5f});
      // get_scalar(command.values, "Kr", pmaterial.specular, 1))
      //   return parse_error();
      get_scalar(command.values, "eta", pmaterial.ior, 1.5f);
      pmaterial.roughness = 0;
      get_roughness(command.values, pmaterial.roughness, 0);
      if (verbose) printf("kdsubsurface material not properly supported\n");
    } else if (command.type == "subsurface") {
      pmaterial.type = pbrt_mtype::subsurface;
      // get_scalar(command.values, "Kr", pmaterial.specular, 1))
      //   return parse_error();
      // get_scalar(command.values, "Kt", pmaterial.transmission, 1))
      //   return parse_error();
      pmaterial.color = {1, 1, 1};
      get_scalar(command.values, "eta", pmaterial.ior, 1.5);
      pmaterial.roughness = 0;
      get_roughness(command.values, pmaterial.roughness, 0);
      auto scale = 1.0f;
      get_pbrt_value(command.values, "scale", scale);
      pmaterial.volscale = 1 / scale;
      auto sigma_a = zero3f, sigma_s = zero3f;
      auto sigma_a_tex = -1, sigma_s_tex = -1;
      get_texture(command.values, "sigma_a", sigma_a, sigma_a_tex,
          vec3f{0.011f, .0024f, .014f});
      get_texture(command.values, "sigma_prime_s", sigma_s, sigma_s_tex,
          vec3f{2.55f, 3.12f, 3.77f});
      pmaterial.volmeanfreepath = 1 / (sigma_a + sigma_s);
      pmaterial.volscatter      = sigma_s / (sigma_a + sigma_s);
      if (verbose) printf("subsurface material not properly supported\n");
    } else if (command.type == "mix") {
      auto namedmaterial1 = ""s, namedmaterial2 = ""s;
      get_pbrt_value(command.values, "namedmaterial1", namedmaterial1);
      get_pbrt_value(command.values, "namedmaterial2", namedmaterial2);
      auto matname = (!namedmaterial1.empty()) ? namedmaterial1
                                               : namedmaterial2;
      auto matit   = named_materials.find(matname);
      if (matit == named_materials.end())
        throw io_error::material_error(filename, matname);
      auto saved_name = pmaterial.name;
      pmaterial       = matit->second;
      pmaterial.name  = saved_name;
      if (verbose) printf("mix material not properly supported\n");
    } else if (command.type == "fourier") {
      auto bsdffile = ""s;
      get_pbrt_value(command.values, "bsdffile", bsdffile);
      if (bsdffile.rfind('/') != string::npos)
        bsdffile = bsdffile.substr(bsdffile.rfind('/') + 1);
      if (bsdffile == "paint.bsdf") {
        pmaterial.type      = pbrt_mtype::plastic;
        pmaterial.color     = {0.6f, 0.6f, 0.6f};
        pmaterial.ior       = 1.5f;
        pmaterial.roughness = 0.2f;
      } else if (bsdffile == "ceramic.bsdf") {
        pmaterial.type      = pbrt_mtype::plastic;
        pmaterial.color     = {0.6f, 0.6f, 0.6f};
        pmaterial.ior       = 1.5f;
        pmaterial.roughness = 0.25f;
      } else if (bsdffile == "leather.bsdf") {
        pmaterial.type      = pbrt_mtype::plastic;
        pmaterial.color     = {0.6f, 0.57f, 0.48f};
        pmaterial.ior       = 1.5f;
        pmaterial.roughness = 0.3f;
      } else if (bsdffile == "coated_copper.bsdf") {
        pmaterial.type  = pbrt_mtype::metal;
        auto eta        = vec3f{0.2004376970f, 0.9240334304f, 1.1022119527f};
        auto etak       = vec3f{3.9129485033f, 2.4528477015f, 2.1421879552f};
        pmaterial.color = eta_to_reflectivity(eta, etak);
        pmaterial.roughness = 0.01f;
      } else if (bsdffile == "roughglass_alpha_0.2.bsdf") {
        pmaterial.type      = pbrt_mtype::glass;
        pmaterial.color     = {1, 1, 1};
        pmaterial.ior       = 1.5f;
        pmaterial.roughness = 0.2f;
      } else if (bsdffile == "roughgold_alpha_0.2.bsdf") {
        pmaterial.type  = pbrt_mtype::metal;
        auto eta        = vec3f{0.1431189557f, 0.3749570432f, 1.4424785571f};
        auto etak       = vec3f{3.9831604247f, 2.3857207478f, 1.6032152899f};
        pmaterial.color = eta_to_reflectivity(eta, etak);
        pmaterial.roughness = 0.2f;
      } else {
        throw io_error::type_error(filename, bsdffile);
      }
    } else {
      throw io_error::type_error(filename, command.type);
    }
  } catch (const io_error& error) {
    throw;
  } catch (...) {
    throw io_error::parse_error(filename);
  }
}

// Make a triangle shape from a quad grid
template <typename PositionFunc, typename NormalFunc>
inline void make_shape(vector<vec3i>& triangles, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    const PositionFunc& position_func, const NormalFunc& normal_func) {
  auto vid = [steps](int i, int j) { return j * (steps.x + 1) + i; };
  auto tid = [steps](int i, int j, int c) { return (j * steps.x + i) * 2 + c; };
  positions.resize((steps.x + 1) * (steps.y + 1));
  normals.resize((steps.x + 1) * (steps.y + 1));
  texcoords.resize((steps.x + 1) * (steps.y + 1));
  for (auto j = 0; j < steps.y + 1; j++) {
    for (auto i = 0; i < steps.x + 1; i++) {
      auto uv              = vec2f{i / (float)steps.x, j / (float)steps.y};
      positions[vid(i, j)] = position_func(uv);
      normals[vid(i, j)]   = normal_func(uv);
      texcoords[vid(i, j)] = uv;
    }
  }
  triangles.resize(steps.x * steps.y * 2);
  for (auto j = 0; j < steps.y; j++) {
    for (auto i = 0; i < steps.x; i++) {
      triangles[tid(i, j, 0)] = {vid(i, j), vid(i + 1, j), vid(i + 1, j + 1)};
      triangles[tid(i, j, 1)] = {vid(i, j), vid(i + 1, j + 1), vid(i, j + 1)};
    }
  }
}

// pbrt sphere
inline void make_sphere(vector<vec3i>& triangles, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    float radius) {
  make_shape(
      triangles, positions, normals, texcoords, steps,
      [radius](const vec2f& uv) {
        auto pt = vec2f{2 * pif * uv.x, pif * (1 - uv.y)};
        return radius *
               vec3f{cos(pt.x) * sin(pt.y), sin(pt.x) * sin(pt.y), cos(pt.y)};
      },
      [](const vec2f& uv) {
        auto pt = vec2f{2 * pif * uv.x, pif * (1 - uv.y)};
        return vec3f{cos(pt.x) * sin(pt.y), sin(pt.x) * sin(pt.y), cos(pt.y)};
      });
}
inline void make_disk(vector<vec3i>& triangles, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    float radius) {
  make_shape(
      triangles, positions, normals, texcoords, steps,
      [radius](const vec2f& uv) {
        auto a = 2 * pif * uv.x;
        return radius * (1 - uv.y) * vec3f{cos(a), sin(a), 0};
      },
      [](const vec2f& uv) {
        return vec3f{0, 0, 1};
      });
}
inline void make_quad(vector<vec3i>& triangles, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    float radius) {
  make_shape(
      triangles, positions, normals, texcoords, steps,
      [radius](const vec2f& uv) {
        return vec3f{(uv.x - 0.5f) * radius, (uv.y - 0.5f) * radius, 0};
      },
      [](const vec2f& uv) {
        return vec3f{0, 0, 1};
      });
}

// Convert pbrt shapes
inline void convert_shape(pbrt_shape& pshape, const pbrt_command& command,
    string& alphamap, const unordered_map<string, pbrt_texture>& named_textures,
    const string& ply_dirname, bool ply_meshes, const string& filename,
    bool verbose = false) {
  // helpers
  auto get_alpha = [&](const vector<pbrt_value>& values, const string& name,
                       string& filename) -> bool {
    auto def      = 1.0f;
    auto textured = pair{def, ""s};
    get_pbrt_value(values, name, textured);
    if (textured.second.empty()) {
      filename = "";
    } else {
      filename = named_textures.at(textured.second).filename;
    }
    return true;
  };

  try {
    pshape.frame = command.frame;
    pshape.frend = command.frend;
    if (command.type == "trianglemesh") {
      pshape.positions = {};
      pshape.normals   = {};
      pshape.texcoords = {};
      pshape.triangles = {};
      get_pbrt_value(command.values, "P", pshape.positions);
      get_pbrt_value(command.values, "N", pshape.normals);
      get_pbrt_value(command.values, "uv", pshape.texcoords);
      for (auto& uv : pshape.texcoords) uv.y = (1 - uv.y);
      get_pbrt_value(command.values, "indices", pshape.triangles);
    } else if (command.type == "loopsubdiv") {
      pshape.positions = {};
      pshape.triangles = {};
      get_pbrt_value(command.values, "P", pshape.positions);
      get_pbrt_value(command.values, "indices", pshape.triangles);
      pshape.normals.resize(pshape.positions.size());
      // compute_normals(pshape.normals, pshape.triangles, pshape.positions);
    } else if (command.type == "plymesh") {
      pshape.filename_ = ""s;
      get_pbrt_value(command.values, "filename", pshape.filename_);
      get_alpha(command.values, "alpha", alphamap);
      if (ply_meshes) {
        auto ply = ply_model{};
        try {
          load_ply(path_join(ply_dirname, pshape.filename_), ply);
        } catch (const io_error& error) {
          throw io_error::dependent_error(filename, error);
        }
        get_positions(ply, pshape.positions);
        get_normals(ply, pshape.normals);
        get_texcoords(ply, pshape.texcoords);
        get_triangles(ply, pshape.triangles);
      }
    } else if (command.type == "sphere") {
      auto radius = 1.0f;
      get_pbrt_value(command.values, "radius", radius);
      make_sphere(pshape.triangles, pshape.positions, pshape.normals,
          pshape.texcoords, {32, 16}, radius);
    } else if (command.type == "disk") {
      auto radius = 1.0f;
      get_pbrt_value(command.values, "radius", radius);
      make_disk(pshape.triangles, pshape.positions, pshape.normals,
          pshape.texcoords, {32, 1}, radius);
    } else {
      throw io_error::type_error(filename, command.type);
    }
  } catch (const io_error& error) {
    throw;
  } catch (...) {
    throw io_error::parse_error(filename);
  }
}

// Convert pbrt arealights
inline void convert_arealight(pbrt_arealight& parealight,
    const pbrt_command& command, const string& filename, bool verbose = false) {
  try {
    parealight.name = command.name;
    if (command.type == "diffuse") {
      auto l = vec3f{1, 1, 1}, scale = vec3f{1, 1, 1};
      get_pbrt_value(command.values, "L", l);
      get_pbrt_value(command.values, "scale", scale);
      parealight.emission = l * scale;
    } else {
      throw io_error::type_error(filename, command.type);
    }
  } catch (const io_error& error) {
    throw;
  } catch (...) {
    throw io_error::parse_error(filename);
  }
}

// Convert pbrt lights
inline void convert_light(pbrt_light& plight, const pbrt_command& command,
    const string& filename, bool verbose = false) {
  try {
    plight.frame = command.frame;
    plight.frend = command.frend;
    if (command.type == "distant") {
      auto l = vec3f{1, 1, 1}, scale = vec3f{1, 1, 1};
      get_pbrt_value(command.values, "L", l);
      get_pbrt_value(command.values, "scale", scale);
      plight.emission = l * scale;
      plight.from     = vec3f{0, 0, 0};
      plight.to       = vec3f{0, 0, 1};
      get_pbrt_value(command.values, "from", plight.from);
      get_pbrt_value(command.values, "to", plight.to);
      plight.distant       = true;
      auto distant_dist    = 100.0f;
      auto size            = distant_dist * sin(5 * pif / 180);
      plight.area_emission = plight.emission * (distant_dist * distant_dist) /
                             (size * size);
      plight.area_frame = plight.frame *
                          lookat_frame(
                              normalize(plight.from - plight.to) * distant_dist,
                              {0, 0, 0}, {0, 1, 0}, true);
      plight.area_frend = plight.frend *
                          lookat_frame(
                              normalize(plight.from - plight.to) * distant_dist,
                              {0, 0, 0}, {0, 1, 0}, true);
      auto texcoords = vector<vec2f>{};
      make_quad(plight.area_triangles, plight.area_positions,
          plight.area_normals, texcoords, {4, 2}, size);
    } else if (command.type == "point" || command.type == "goniometric" ||
               command.type == "spot") {
      auto i = vec3f{1, 1, 1}, scale = vec3f{1, 1, 1};
      get_pbrt_value(command.values, "I", i);
      get_pbrt_value(command.values, "scale", scale);
      plight.emission = i * scale;
      plight.from     = zero3f;
      get_pbrt_value(command.values, "from", plight.from);
      plight.area_emission = plight.emission;
      plight.area_frame    = plight.frame * translation_frame(plight.from);
      plight.area_frend    = plight.frend * translation_frame(plight.from);
      auto texcoords       = vector<vec2f>{};
      make_sphere(plight.area_triangles, plight.area_positions,
          plight.area_normals, texcoords, {4, 2}, 0.0025f);
    } else {
      throw io_error::type_error(filename, command.type);
    }
  } catch (const io_error& error) {
    throw;
  } catch (...) {
    throw io_error::parse_error(filename);
  }
}

inline void convert_environment(pbrt_environment& penvironment,
    const pbrt_command& command, unordered_map<string, int>& texture_map,
    const string& filename, bool verbose = false) {
  penvironment.frame = command.frame;
  penvironment.frend = command.frend;
  penvironment.frame = penvironment.frame *
                       frame3f{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 0, 0}};
  penvironment.frend = penvironment.frend *
                       frame3f{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 0, 0}};
  try {
    if (command.type == "infinite") {
      auto l = vec3f{1, 1, 1}, scale = vec3f{1, 1, 1};
      get_pbrt_value(command.values, "L", l);
      get_pbrt_value(command.values, "scale", scale);
      penvironment.emission     = scale * l;
      penvironment.emission_tex = -1;
      auto mapname              = ""s;
      get_pbrt_value(command.values, "mapname", mapname);
      if (!mapname.empty()) {
        if (texture_map.find(mapname) == texture_map.end()) {
          auto texture_id      = (int)texture_map.size();
          texture_map[mapname] = texture_id;
        }
        penvironment.emission_tex = texture_map.at(mapname);
      }
    } else {
      throw io_error::type_error(filename, command.type);
    }
  } catch (const io_error& error) {
    throw;
  } catch (...) {
    throw io_error::parse_error(filename);
  }
}

// pbrt stack ctm
struct pbrt_stack_element {
  frame3f        transform_start        = identity3x4f;
  frame3f        transform_end          = identity3x4f;
  pbrt_material  material               = {};
  pbrt_arealight arealight              = {};
  pbrt_medium    interior               = {};
  pbrt_medium    exterior               = {};
  bool           reverse                = false;
  bool           active_transform_start = true;
  bool           active_transform_end   = true;
};

// pbrt parsing context
struct pbrt_context {
  vector<pbrt_stack_element>                stack           = {};
  unordered_map<string, pbrt_stack_element> coordsys        = {};
  string                                    cur_object      = "";
  vec2i                                     film_resolution = {512, 512};
};

// load pbrt
inline void load_pbrt(const string& filename, pbrt_model& pbrt,
    pbrt_context& ctx, unordered_map<string, int>& material_map,
    unordered_map<string, int>&           texture_map,
    unordered_map<string, pbrt_material>& named_materials,
    unordered_map<string, pbrt_texture>&  named_textures,
    unordered_map<string, pbrt_medium>&   named_mediums,
    unordered_map<string, vector<int>>&   named_objects,
    const string& ply_dirname, bool ply_meshes) {
  // open file
  auto fs = open_file(filename, "rt");

  // helpers
  auto set_transform = [](pbrt_stack_element& ctx, const frame3f& xform) {
    if (ctx.active_transform_start) ctx.transform_start = xform;
    if (ctx.active_transform_end) ctx.transform_end = xform;
  };
  auto concat_transform = [](pbrt_stack_element& ctx, const frame3f& xform) {
    if (ctx.active_transform_start) ctx.transform_start *= xform;
    if (ctx.active_transform_end) ctx.transform_end *= xform;
  };

  // init stack
  if (ctx.stack.empty()) ctx.stack.emplace_back();

  // parse command by command
  auto line = ""s;
  while (read_pbrt_cmdline(fs, line)) {
    try {
      auto str = string_view{line};
      // get command
      auto cmd = ""s;
      parse_command(str, cmd);
      if (cmd == "WorldBegin") {
        ctx.stack.push_back({});
      } else if (cmd == "WorldEnd") {
        if (ctx.stack.empty()) throw std::out_of_range{"invalid stack"};
        ctx.stack.pop_back();
        if (ctx.stack.size() != 1) throw std::out_of_range{"invalid stack"};
      } else if (cmd == "AttributeBegin") {
        ctx.stack.push_back(ctx.stack.back());
      } else if (cmd == "AttributeEnd") {
        if (ctx.stack.empty()) throw std::out_of_range{"invalid stack"};
        ctx.stack.pop_back();
      } else if (cmd == "TransformBegin") {
        ctx.stack.push_back(ctx.stack.back());
      } else if (cmd == "TransformEnd") {
        if (ctx.stack.empty()) throw std::out_of_range{"invalid stack"};
        ctx.stack.pop_back();
      } else if (cmd == "ObjectBegin") {
        ctx.stack.push_back(ctx.stack.back());
        parse_param(str, ctx.cur_object);
        named_objects[ctx.cur_object] = {};
      } else if (cmd == "ObjectEnd") {
        ctx.stack.pop_back();
        ctx.cur_object = "";
      } else if (cmd == "ObjectInstance") {
        auto object = ""s;
        parse_param(str, object);
        if (named_objects.find(object) == named_objects.end())
          throw io_error::object_error(filename, object);
        auto& named_object = named_objects.at(object);
        for (auto& shape_id : named_object) {
          pbrt.shapes[shape_id].instances.push_back(
              ctx.stack.back().transform_start);
          pbrt.shapes[shape_id].instaends.push_back(
              ctx.stack.back().transform_end);
        }
      } else if (cmd == "ActiveTransform") {
        auto name = ""s;
        parse_command(str, name);
        if (name == "StartTime") {
          ctx.stack.back().active_transform_start = true;
          ctx.stack.back().active_transform_end   = false;
        } else if (name == "EndTime") {
          ctx.stack.back().active_transform_start = false;
          ctx.stack.back().active_transform_end   = true;
        } else if (name == "All") {
          ctx.stack.back().active_transform_start = true;
          ctx.stack.back().active_transform_end   = true;
        } else {
          std::out_of_range{"invalid command"};
        }
      } else if (cmd == "Transform") {
        auto xf = identity4x4f;
        parse_param(str, xf);
        set_transform(ctx.stack.back(), mat_to_frame(xf));
      } else if (cmd == "ConcatTransform") {
        auto xf = identity4x4f;
        parse_param(str, xf);
        concat_transform(ctx.stack.back(), mat_to_frame(xf));
      } else if (cmd == "Scale") {
        auto v = zero3f;
        parse_param(str, v);
        concat_transform(ctx.stack.back(), scaling_frame(v));
      } else if (cmd == "Translate") {
        auto v = zero3f;
        parse_param(str, v);
        concat_transform(ctx.stack.back(), translation_frame(v));
      } else if (cmd == "Rotate") {
        auto v = zero4f;
        parse_param(str, v);
        concat_transform(ctx.stack.back(),
            rotation_frame(vec3f{v.y, v.z, v.w}, radians(v.x)));
      } else if (cmd == "LookAt") {
        auto from = zero3f, to = zero3f, up = zero3f;
        parse_param(str, from);
        parse_param(str, to);
        parse_param(str, up);
        auto frame = lookat_frame(from, to, up, true);
        concat_transform(ctx.stack.back(), inverse(frame));
      } else if (cmd == "ReverseOrientation") {
        ctx.stack.back().reverse = !ctx.stack.back().reverse;
      } else if (cmd == "CoordinateSystem") {
        auto name = ""s;
        parse_param(str, name);
        ctx.coordsys[name].transform_start = ctx.stack.back().transform_start;
        ctx.coordsys[name].transform_end   = ctx.stack.back().transform_end;
      } else if (cmd == "CoordSysTransform") {
        auto name = ""s;
        parse_param(str, name);
        if (ctx.coordsys.find(name) != ctx.coordsys.end()) {
          ctx.stack.back().transform_start =
              ctx.coordsys.at(name).transform_start;
          ctx.stack.back().transform_end = ctx.coordsys.at(name).transform_end;
        }
      } else if (cmd == "Integrator") {
        auto command = pbrt_command{};
        parse_param(str, command.type);
        parse_params(str, command.values);
      } else if (cmd == "Sampler") {
        auto command = pbrt_command{};
        parse_param(str, command.type);
        parse_params(str, command.values);
      } else if (cmd == "PixelFilter") {
        auto command = pbrt_command{};
        parse_param(str, command.type);
        parse_params(str, command.values);
      } else if (cmd == "Film") {
        auto command = pbrt_command{};
        parse_param(str, command.type);
        parse_params(str, command.values);
        auto film = pbrt_film{};
        convert_film(film, command, filename);
        ctx.film_resolution = film.resolution;
      } else if (cmd == "Accelerator") {
        auto command = pbrt_command{};
        parse_param(str, command.type);
        parse_params(str, command.values);
      } else if (cmd == "Camera") {
        auto command = pbrt_command{};
        parse_param(str, command.type);
        parse_params(str, command.values);
        command.frame = ctx.stack.back().transform_start;
        command.frend = ctx.stack.back().transform_end;
        auto& camera  = pbrt.cameras.emplace_back();
        convert_camera(camera, command, ctx.film_resolution, filename);
      } else if (cmd == "Texture") {
        auto command  = pbrt_command{};
        auto comptype = ""s;
        auto str_     = string{str};
        parse_param(str, command.name);
        parse_param(str, comptype);
        parse_param(str, command.type);
        parse_params(str, command.values);
        convert_texture(
            named_textures[command.name], command, named_textures, filename);
      } else if (cmd == "Material") {
        static auto material_id = 0;
        auto        command     = pbrt_command{};
        command.name = "__unnamed__material__" + std::to_string(material_id++);
        parse_param(str, command.type);
        parse_params(str, command.values);
        if (command.type.empty()) {
          ctx.stack.back().material = {};
        } else {
          ctx.stack.back().material = {};
          convert_material(ctx.stack.back().material, command, texture_map,
              named_materials, named_textures, filename);
        }
      } else if (cmd == "MakeNamedMaterial") {
        auto command = pbrt_command{};
        parse_param(str, command.name);
        parse_params(str, command.values);
        command.type = "";
        for (auto& value : command.values)
          if (value.name == "type") command.type = value.value1s;
        convert_material(named_materials[command.name], command, texture_map,
            named_materials, named_textures, filename);
      } else if (cmd == "NamedMaterial") {
        auto name = ""s;
        parse_param(str, name);
        if (named_materials.find(name) == named_materials.end())
          throw io_error::material_error(filename, name);
        ctx.stack.back().material = named_materials.at(name);
      } else if (cmd == "Shape") {
        auto command = pbrt_command{};
        parse_param(str, command.type);
        parse_params(str, command.values);
        command.frame  = ctx.stack.back().transform_start;
        command.frend  = ctx.stack.back().transform_end;
        auto& shape    = pbrt.shapes.emplace_back();
        auto  alphamap = ""s;
        convert_shape(shape, command, alphamap, named_textures, ply_dirname,
            ply_meshes, filename);
        auto matkey = "?!!!?" + ctx.stack.back().material.name + "?!!!?" +
                      ctx.stack.back().arealight.name + "?!!!?" + alphamap;
        if (material_map.find(matkey) == material_map.end()) {
          auto& material = pbrt.materials.emplace_back();
          material       = ctx.stack.back().material;
          material.name  = "material" + std::to_string(pbrt.materials.size());
          material.emission = ctx.stack.back().arealight.emission;
          // material.alpha_tex = alphamap;
          material_map[matkey] = (int)pbrt.materials.size() - 1;
        }
        shape.material = material_map.at(matkey);
        if (!ctx.cur_object.empty()) {
          named_objects[ctx.cur_object].push_back((int)pbrt.shapes.size() - 1);
          shape.instanced = true;
        }
      } else if (cmd == "AreaLightSource") {
        static auto arealight_id = 0;
        auto        command      = pbrt_command{};
        command.name             = "__unnamed__arealight__" +
                       std::to_string(arealight_id++);
        parse_param(str, command.type);
        parse_params(str, command.values);
        command.frame = ctx.stack.back().transform_start;
        command.frend = ctx.stack.back().transform_end;
        convert_arealight(ctx.stack.back().arealight, command, filename);
      } else if (cmd == "LightSource") {
        auto command = pbrt_command{};
        parse_param(str, command.type);
        parse_params(str, command.values);
        command.frame = ctx.stack.back().transform_start;
        command.frend = ctx.stack.back().transform_end;
        if (command.type == "infinite") {
          auto& environment = pbrt.environments.emplace_back();
          convert_environment(environment, command, texture_map, filename);
        } else {
          auto& light = pbrt.lights.emplace_back();
          convert_light(light, command, filename);
        }
      } else if (cmd == "MakeNamedMedium") {
        auto command = pbrt_command{};
        parse_param(str, command.name);
        parse_params(str, command.values);
        command.type = "";
        for (auto& value : command.values)
          if (command.name == "type") command.type = value.value1s;
        auto medium                 = pbrt_medium{};
        named_mediums[command.name] = medium;
      } else if (cmd == "MediumInterface") {
        auto interior = ""s, exterior = ""s;
        parse_param(str, interior);
        parse_param(str, exterior);
        ctx.stack.back().interior = named_mediums.at(interior);
        ctx.stack.back().exterior = named_mediums.at(exterior);
      } else if (cmd == "Include") {
        auto includename = ""s;
        parse_param(str, includename);
        try {
          load_pbrt(path_join(path_dirname(filename), includename), pbrt, ctx,
              material_map, texture_map, named_materials, named_textures,
              named_mediums, named_objects, ply_dirname, ply_meshes);
        } catch (const io_error& error) {
          throw io_error::dependent_error(filename, error);
        }
      } else {
        throw io_error::command_error(filename, cmd);
      }
    } catch (const io_error& error) {
      throw;
    } catch (...) {
      throw io_error::parse_error(filename);
    }
  }
}

// Load pbrt
pbrt_model load_pbrt(const string& filename, bool ply_meshes) {
  auto pbrt = pbrt_model{};
  load_pbrt(filename, pbrt, ply_meshes);
  return pbrt;
}

// load pbrt
void load_pbrt(const string& filename, pbrt_model& pbrt, bool ply_meshes) {
  auto ctx             = pbrt_context{};
  auto material_map    = unordered_map<string, int>{};
  auto texture_map     = unordered_map<string, int>{};
  auto named_materials = unordered_map<string, pbrt_material>{{"", {}}};
  auto named_mediums   = unordered_map<string, pbrt_medium>{{"", {}}};
  auto named_textures  = unordered_map<string, pbrt_texture>{{"", {}}};
  auto named_objects   = unordered_map<string, vector<int>>{};
  load_pbrt(filename, pbrt, ctx, material_map, texture_map, named_materials,
      named_textures, named_mediums, named_objects, path_dirname(filename),
      ply_meshes);
  pbrt.textures.resize(texture_map.size());
  for (auto& [path, texture_id] : texture_map) {
    pbrt.textures[texture_id].filename = path;
  }
}

inline void format_value(string& str, const pbrt_value& value) {
  static auto type_labels = unordered_map<pbrt_type, string>{
      {pbrt_type::real, "float"},
      {pbrt_type::integer, "integer"},
      {pbrt_type::boolean, "bool"},
      {pbrt_type::string, "string"},
      {pbrt_type::point, "point"},
      {pbrt_type::normal, "normal"},
      {pbrt_type::vector, "vector"},
      {pbrt_type::texture, "texture"},
      {pbrt_type::color, "rgb"},
      {pbrt_type::point2, "point2"},
      {pbrt_type::vector2, "vector2"},
      {pbrt_type::spectrum, "spectrum"},
  };

  auto format_vector = [](string& str, auto& values) {
    str += "[ ";
    for (auto& value : values) {
      str += " ";
      format_value(str, value);
    }
    str += " ]";
  };

  format_values(str, "\"{} {}\" ", type_labels.at(value.type), value.name);
  switch (value.type) {
    case pbrt_type::real:
      if (!value.vector1f.empty()) {
        format_vector(str, value.vector1f);
      } else {
        format_value(str, value.value1f);
      }
      break;
    case pbrt_type::integer:
      if (!value.vector1f.empty()) {
        format_vector(str, value.vector1i);
      } else {
        format_value(str, value.value1i);
      }
      break;
    case pbrt_type::boolean:
      format_values(str, "\"{}\"", value.value1b ? "true" : "false");
      break;
    case pbrt_type::string:
    case pbrt_type::texture: format_values(str, "\"{}\"", value.value1s); break;
    case pbrt_type::point:
    case pbrt_type::vector:
    case pbrt_type::normal:
    case pbrt_type::color:
      if (!value.vector3f.empty()) {
        format_vector(str, value.vector3f);
      } else {
        format_values(str, "[ {} ]", value.value3f);
      }
      break;
    case pbrt_type::spectrum: format_vector(str, value.vector1f); break;
    case pbrt_type::point2:
    case pbrt_type::vector2:
      if (!value.vector2f.empty()) {
        format_vector(str, value.vector2f);
      } else {
        format_values(str, "[ {} ]", value.value2f);
      }
      break;
  }
}

inline void format_value(string& str, const vector<pbrt_value>& values) {
  for (auto& value : values) {
    str += " ";
    format_value(str, value);
  }
}

void save_pbrt(
    const string& filename, const pbrt_model& pbrt, bool ply_meshes) {
  // open file
  auto fs = open_file(filename, "wt");

  // save comments
  format_values(fs, "#\n");
  format_values(fs, "# Written by Yocto/GL\n");
  format_values(fs, "# https://github.com/xelatihy/yocto-gl\n");
  format_values(fs, "#\n\n");
  for (auto& comment : pbrt.comments) {
    format_values(fs, "# {}\n", comment);
  }
  format_values(fs, "\n");

  for (auto& camera : pbrt.cameras) {
    auto command = pbrt_command{};
    command.type = "image";
    command.values.push_back(
        make_pbrt_value("xresolution", camera.resolution.x));
    command.values.push_back(
        make_pbrt_value("yresolution", camera.resolution.y));
    command.values.push_back(make_pbrt_value("filename", "image.exr"s));
    format_values(fs, "Film \"{}\" {}\n", command.type, command.values);
  }

  for (auto& camera : pbrt.cameras) {
    auto command  = pbrt_command{};
    command.type  = "perspective";
    command.frame = camera.frame;
    command.values.push_back(make_pbrt_value(
        "fov", 2 * tan(0.036f / (2 * camera.lens)) * 180 / pif));
    format_values(fs, "LookAt {} {} {}\n", command.frame.o,
        command.frame.o - command.frame.z, command.frame.y);
    format_values(fs, "Camera \"{}\" {}\n", command.type, command.values);
  }

  format_values(fs, "\nWorldBegin\n\n");

  for (auto& light : pbrt.lights) {
    auto command  = pbrt_command{};
    command.frame = light.frame;
    if (light.distant) {
      command.type = "distance";
      command.values.push_back(make_pbrt_value("L", light.emission));
    } else {
      command.type = "point";
      command.values.push_back(make_pbrt_value("I", light.emission));
    }
    format_values(fs, "AttributeBegin\n");
    format_values(fs, "Transform {}\n", frame_to_mat(command.frame));
    format_values(fs, "LightSource \"{}\" {}\n", command.type, command.values);
    format_values(fs, "AttributeEnd\n");
  }

  for (auto& environment : pbrt.environments) {
    auto command  = pbrt_command{};
    command.frame = environment.frame;
    command.type  = "infinite";
    command.values.push_back(make_pbrt_value("L", environment.emission));
    command.values.push_back(
        make_pbrt_value("mapname", environment.emission_tex));
    format_values(fs, "AttributeBegin\n");
    format_values(fs, "Transform {}\n", frame_to_mat(command.frame));
    format_values(fs, "LightSource \"{}\" {}\n", command.type, command.values);
    format_values(fs, "AttributeEnd\n");
  }

  auto reflectivity_to_eta = [](const vec3f& reflectivity) {
    return (1 + sqrt(reflectivity)) / (1 - sqrt(reflectivity));
  };

  for (auto& material : pbrt.materials) {
    auto command = pbrt_command{};
    switch (material.type) {
      case pbrt_mtype::matte: {
        command.type = "matte";
        command.values.push_back(make_pbrt_value("Kd", material.color));
      } break;
      case pbrt_mtype::plastic: {
        command.type = "matte";
        command.values.push_back(make_pbrt_value("Kd", material.color));
        command.values.push_back(make_pbrt_value("Ks", vec3f{1, 1, 1}));
        command.values.push_back(
            make_pbrt_value("roughness", pow(material.roughness, 2)));
        command.values.push_back(
            make_pbrt_value("eta", reflectivity_to_eta(material.color)));
        command.values.push_back(make_pbrt_value("remaproughness", false));
      } break;
      case pbrt_mtype::metal: {
        command.type = "metal";
        command.values.push_back(make_pbrt_value("Kr", vec3f{1, 1, 1}));
        command.values.push_back(
            make_pbrt_value("roughness", pow(material.roughness, 2)));
        command.values.push_back(
            make_pbrt_value("eta", reflectivity_to_eta(material.color)));
        command.values.push_back(make_pbrt_value("remaproughness", false));
      } break;
      case pbrt_mtype::thinglass: {
        command.type = "uber";
        command.values.push_back(make_pbrt_value("Ks", vec3f{1, 1, 1}));
        command.values.push_back(make_pbrt_value("Kt", material.color));
        command.values.push_back(
            make_pbrt_value("roughness", pow(material.roughness, 2)));
        command.values.push_back(
            make_pbrt_value("eta", reflectivity_to_eta(material.color)));
        command.values.push_back(make_pbrt_value("remaproughness", false));
      } break;
      case pbrt_mtype::glass: {
        command.type = "glass";
        command.values.push_back(make_pbrt_value("Kr", vec3f{1, 1, 1}));
        command.values.push_back(make_pbrt_value("Kt", vec3f{1, 1, 1}));
        command.values.push_back(
            make_pbrt_value("roughness", pow(material.roughness, 2)));
        command.values.push_back(make_pbrt_value("eta", material.ior));
        command.values.push_back(make_pbrt_value("remaproughness", false));
      } break;
      case pbrt_mtype::subsurface: {
        command.type = "matte";
        command.values.push_back(make_pbrt_value("Kd", material.color));
      } break;
    }

    format_values(fs, "MakeNamedMaterial \"{}\" \"string type\" \"{}\" {}\n",
        material.name, command.type, command.values);
  }

  auto object_id = 0;
  for (auto& shape : pbrt.shapes) {
    auto& material = pbrt.materials.at(shape.material);
    auto  command  = pbrt_command{};
    command.frame  = shape.frame;
    if (ply_meshes) {
      command.type = "plymesh";
      command.values.push_back(make_pbrt_value("filename", shape.filename_));
    } else {
      command.type = "trianglemesh";
      command.values.push_back(make_pbrt_value("indices", shape.triangles));
      command.values.push_back(
          make_pbrt_value("P", shape.positions, pbrt_type::point));
      if (!shape.normals.empty())
        command.values.push_back(
            make_pbrt_value("N", shape.triangles, pbrt_type::normal));
      if (!shape.texcoords.empty())
        command.values.push_back(make_pbrt_value("uv", shape.texcoords));
    }
    if (ply_meshes) {
      auto ply = ply_model{};
      add_positions(ply, shape.positions);
      add_normals(ply, shape.normals);
      add_texcoords(ply, shape.texcoords);
      add_triangles(ply, shape.triangles);
      try {
        save_ply(path_dirname(filename) + "/" + shape.filename_, ply);
      } catch (const io_error& error) {
        throw io_error::dependent_error(filename, error);
      }
    }
    auto object = "object" + std::to_string(object_id++);
    if (!shape.instances.empty())
      format_values(fs, "ObjectBegin \"{}\"\n", object);
    format_values(fs, "AttributeBegin\n");
    format_values(fs, "Transform {}\n", frame_to_mat(shape.frame));
    if (material.emission != zero3f) {
      auto acommand = pbrt_command{};
      acommand.type = "diffuse";
      acommand.values.push_back(make_pbrt_value("L", material.emission));
      format_values(
          fs, "AreaLightSource \"{}\" {}\n", acommand.type, acommand.values);
    }
    format_values(fs, "NamedMaterial \"{}\"\n", material.name);
    format_values(fs, "Shape \"{}\" {}\n", command.type, command.values);
    format_values(fs, "AttributeEnd\n");
    if (!shape.instances.empty()) format_values(fs, "ObjectEnd\n");
    for (auto& iframe : shape.instances) {
      format_values(fs, "AttributeBegin\n");
      format_values(fs, "Transform {}\n", frame_to_mat(iframe));
      format_values(fs, "ObjectInstance \"{}\"\n", object);
      format_values(fs, "AttributeEnd\n");
    }
  }

  format_values(fs, "\nWorldEnd\n\n");
}

// Load/save pbrt
bool load_pbrt(
    const string& filename, pbrt_model& pbrt, string& error, bool ply_meshes) {
  try {
    load_pbrt(filename, pbrt, ply_meshes);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}
bool save_pbrt(const string& filename, const pbrt_model& pbrt, string& error,
    bool ply_meshes) {
  try {
    save_pbrt(filename, pbrt, ply_meshes);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

}  // namespace yocto
