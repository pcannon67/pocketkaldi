// Created at 2017-4-19

#include "configuration.h"

#include <string>
#include <vector>
#include "util.h"
#include "status.h"

namespace pocketkaldi {

const std::string Configuration::kDefaultString = "~_$$NOT_Ext_Val$$>-^-<@@_";

Status Configuration::Read(const std::string &filename) {
  util::ReadableFile fd;
  std::string line;
  std::vector<std::string> fields;
  std::string key;
  std::string value;
  Status status;

  filename_ = filename;
  status = fd.Open(filename);
  if (!status.ok()) return status;

  // Read line by line
  while (!fd.Eof() && fd.ReadLine(&line, &status)) {
    line = util::Trim(line);

    // Empty line and comment line
    if (line.empty() || line.front() == '#') continue;

    // Parse the line with key=value format
    fields = util::Split(line, "=");
    if (fields.size() != 2) {
      return Status::Corruption(util::Format(
          "Unexpected line in {}: {}",
          filename_,
          line));
    }
    key = util::Tolower(util::Trim(fields[0]));
    value = util::Trim(fields[1]);
    if (value.empty()) {
      return Status::Corruption(util::Format(
          "Value cound not be empty: {}",
          filename_,
          line));
    }

    table_[key] = value;
  }

  return status;
}

std::string Configuration::GetPathOrElse(
    const std::string &key,
    const std::string &default_val) const {
  // Currently it is for *nix only
  std::string path = GetStringOrElse(key, kDefaultString);
  if (path == kDefaultString) return default_val;

  // Return the path directly if it is a absolute path
  if (path.front() == '/') return path; 

  // Get the directory name of filename_
  int pos = filename_.rfind('/');
  if (pos == std::string::npos) return path; 
  std::string directory(filename_.cbegin(), filename_.cbegin() + pos + 1);
  return directory + path;
}

std::string Configuration::GetStringOrElse(
    const std::string &key,
    const std::string &default_val) const {
  std::string lower_key = util::Tolower(key);
  std::unordered_map<std::string, std::string>::const_iterator
  it = table_.find(lower_key);
  if (it == table_.end()) return default_val;

  return it->second;
}

int Configuration::GetIntegerOrElse(const std::string &key, int default_val) const {
  std::string int_string = GetStringOrElse(key, kDefaultString);
  if (int_string == kDefaultString) return default_val;
  return std::stoi(int_string);
}

}  // namespace pocketkaldi
