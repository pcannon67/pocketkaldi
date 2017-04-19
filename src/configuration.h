// Created at 2017-4-19

#ifndef POCKETKALDI_CONFIGURATION_H_
#define POCKETKALDI_CONFIGURATION_H_

#include <string>
#include <unordered_map>
#include "util.h"

namespace pocketkaldi {

// Configuration is a class to read configuration files for pocketkaldi.
class Configuration {
 public:
  // Read configuration from filename. If success, `status->ok` will be true
  // Otherwise, status->message will contain error message
  void Read(const std::string &filename, Status *status);

  // Get a path from configuration file. If the path is a relative path, return
  // a path combines path to configuatin file and path itself. If the key not
  // exist, returns `default`
  std::string GetPath(
  	  const std::string &key,
  	  const std::string &default_val) const;

 private:
  std::string filename_;
  std::unordered_map<std::string, std::string> table_;
};

}  // namespace pocketkaldi

#endif  // POCKETKALDI_CONFIGURATION_H_