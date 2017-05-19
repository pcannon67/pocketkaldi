// Created at 2014-04-19

#include "configuration.h"

#include <stdio.h>
#include <assert.h>
#include <string>
#include "status.h"
#include "util.h"

void TestConf() {
  std::string conf_file = TESTDIR "data/test.conf";
  pocketkaldi::Status status;
  pocketkaldi::Configuration conf;

  // Read
  status = conf.Read(conf_file);
  puts(status.what().c_str());
  assert(status.ok());

  // Get path of test_conf.txt
  std::string path = conf.GetPath("TestConf", "");
  assert(path != "");

  // Check the path
  pocketkaldi::util::ReadableFile fd;
  std::string line;
  status = fd.Open(path);
  assert(status.ok());
  assert(fd.ReadLine(&line, &status));
  assert(line == "Success!");
}

int main() {
  TestConf();
  return 0;
}
