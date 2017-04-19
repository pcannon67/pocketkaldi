//
// The MIT License (MIT)
//
// Copyright 2013-2014 The MilkCat Project Developers
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
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
// status.h --- Created at 2013-12-06
//


#ifndef POCKETKALDI_STATUS_H_
#define POCKETKALDI_STATUS_H_

#include <stdio.h>
#include <assert.h>
#include <string>

namespace pocketkaldi {

class Status {
 public:
  Status(): code_(0) {}

  static Status OK() { return Status(); }
  static Status IOError(const std::string &message) {
    return Status(kIOError, message);
  }
  static Status Corruption(const std::string &message) {
    return Status(kCorruption, message);
  }
  static Status NotImplemented(const std::string &message) {
    return Status(kNotImplemented, message);
  }
  static Status RuntimeError(const std::string &message) {
    return Status(kRuntimeError, message);
  }
  static Status Info(const std::string &message) {
    return Status(kInfo, message);
  }

  // Return true if the state is success
  bool ok() { return code_ == 0; }

  // Return a string representation of what has happened. If the status
  // is success return a string with length 0
  const std::string &what() const { return message_; }

 private:
  int code_;
  std::string message_;

  enum {
    kIOError = 1,
    kCorruption = 2,
    kRuntimeError = 3,
    kNotImplemented = 4,
    kInfo = 5
  };

  Status(int code, const std::string &error_message): code_(code) {
    std::string message;
    switch (code) {
     case kIOError:
      message = "IOError: ";
      break;
     case kCorruption:
      message = "Corruption: ";
      break;
     case kRuntimeError:
      message = "RuntimeError: ";
      break;
     case kNotImplemented:
      message = "NotImplemented: ";
      break;
     case kInfo:
      message = "";
      break;
    }

    message += error_message;
    message_ = message;
  }
};

}  // namespace pocketkaldi

#endif  // POCKETKALDI_STATUS_H_
