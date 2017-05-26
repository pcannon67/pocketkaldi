// Created at 2017-4-19

#ifndef POCKETKALDI_CONFIGURATION_H_
#define POCKETKALDI_CONFIGURATION_H_

#include <string>
#include <unordered_map>
#include <limits.h>
#include "util.h"

namespace pocketkaldi {

// Configuration is a class to read configuration files for pocketkaldi.
class Configuration {
 public:
  // Read configuration from filename. If success, `status->ok` will be true
  // Otherwise, status->message will contain error message
  Status Read(const std::string &filename);

  // Get a path from configuration file. If the path is a relative path, return
  // a path combines path to configuatin file and path itself. If the key not
  // exist, returns `default`
  std::string GetPathOrElse(
      const std::string &key,
      const std::string &default_val) const;

  // Get a different type of value from configuration file. If the key not
  // exist, returns `default_val`
  std::string GetStringOrElse(
      const std::string &key,
      const std::string &default_val) const;
  int GetIntegerOrElse(const std::string &key, int default_val) const;

  // Return Status version
  Status GetPath(const std::string &key, std::string *val) const {
    std::string path = GetPathOrElse(key, kDefaultString);
    if (path == kDefaultString) return KeyNotFound(key);
    *val = path;
    return Status::OK();
  }
  Status GetString(const std::string &key, std::string *val) const {
    std::string path = GetStringOrElse(key, kDefaultString);
    if (path == kDefaultString) return KeyNotFound(key);
    *val = path;
    return Status::OK();
  }
  Status GetInteger(const std::string &key, int *val) const {
    int int_val = GetIntegerOrElse(key, INT_MIN);
    if (int_val == LONG_MIN) return KeyNotFound(key);
    *val = int_val;
    return Status::OK();
  }

  // Get filename
  const std::string &filename() const {
    return filename_;
  }

 private:
  static const std::string kDefaultString;
  
  std::string filename_;
  std::unordered_map<std::string, std::string> table_;

  Status KeyNotFound(const std::string &key) const {
    return Status::Corruption(util::Format(
        "Unable to find key '{}' in '{}'",
        key,
        filename()));
  }
};

}  // namespace pocketkaldi

#endif  // POCKETKALDI_CONFIGURATION_H_