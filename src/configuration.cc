// Created at 2017-4-19

#include "configuration.h"

#include <string>
#include <vector>
#include "util.h"
#include "status.h"

namespace pocketkaldi {

void Configuration::Read(const std::string &filename, Status *status) {
  util::ReadableFile fd;
  std::string line;
  std::vector<std::string> fields;
  std::string key;
  std::string value;

  filename_ = filename;
  fd.Open(filename, status);
  if (!status->ok()) goto Configuration_Read_Failed;

  // Read line by line
  while (!fd.Eof() && fd.ReadLine(&line, status)) {
    line = util::Trim(line);

    // Empty line and comment line
    if (line.empty() || line.front() == '#') continue;

    // Parse the line with key=value format
    fields = util::Split(line, "=");
    if (fields.size() != 2) {
      *status = Status::Corruption(util::Format(
          "Unexpected line in {}: {}",
          filename_,
          line));
      goto Configuration_Read_Failed;
    }
    key = util::Tolower(util::Trim(fields[0]));
    value = util::Trim(fields[1]);
    if (value.empty()) {
      *status = Status::Corruption(util::Format(
          "Value cound not be empty: {}",
          filename_,
          line));
      goto Configuration_Read_Failed;
    }

    table_[key] = value;
  }

Configuration_Read_Failed:
  return;
}

std::string Configuration::GetPath(
    const std::string &key,
    const std::string &default_val) const {
  std::string lower_key = util::Tolower(key);
  std::unordered_map<std::string, std::string>::const_iterator
  it = table_.find(lower_key);
  if (it == table_.end()) return default_val;

  // Currently it is for *nix only
  std::string path = it->second;

  // Return the path directly if it is a absolute path
  if (path.front() == '/') return path; 

  // Get the directory name of filename_
  int pos = filename_.rfind('/');
  if (pos == std::string::npos) return path; 
  std::string directory(filename_.cbegin(), filename_.cbegin() + pos + 1);
  return directory + path;
}

}  // namespace pocketkaldi
