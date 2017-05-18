// Created at 2016-11-24

#ifndef POCKETKALDI_UTIL_H_
#define POCKETKALDI_UTIL_H_

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <iostream>
#include <utility>
#include <vector>
#include "pocketkaldi.h"
#include "status.h"

// Error types for status
#define PK_STATUS_STSUCCESS 0
#define PK_STATUS_STIOERROR 1
#define PK_STATUS_STCORRUPTED 2

#define PK_WARN fprintf
#define PK_UNUSED(x) (void)(x)
#define PK_MIN(a, b) ((a) < (b) ? (a) : (b))

#define PK_PATHMAX 1024

#define PK_INFO(msg) std::cout << __FILE__ << ": " << (msg) << std::endl;
// #define PK_DEBUG(msg) std::cout << __FILE__ << ": " << (msg) << std::endl;
#define PK_DEBUG(msg)

typedef struct pk_readable_t {
  char filename[PK_PATHMAX];
  int64_t filesize;
  FILE *fd;
} pk_readable_t;

// Stores a byte array as buffer and supports read all kinds of data (like
// string, uint8, int16, ...)
typedef struct pk_bytebuffer_t {
  char *buffer;
  int64_t size; 
  int64_t current_position;
} pk_bytebuffer_t;

// Initialize the status set to success (ok) state
POCKETKALDI_EXPORT
void pk_status_init(pk_status_t *status);

// Set status to failed state with message
POCKETKALDI_EXPORT
void pk_status_fail(pk_status_t *status, int errcode, const char *fmsg, ...);

#define PK_STATUS_IOERROR(status, fmsg, ...) \
    pk_status_fail(status, PK_STATUS_STIOERROR, fmsg, __VA_ARGS__)

#define PK_STATUS_CORRUPTED(status, fmsg, ...) \
    pk_status_fail(status, PK_STATUS_STCORRUPTED, fmsg, __VA_ARGS__)

// Allocate a memory block contains at least `size` bytes
POCKETKALDI_EXPORT
void *pk_alloc(size_t size);

// Allocate a memory block contains at least `size` bytes
POCKETKALDI_EXPORT
void *pk_realloc(void *ptr, size_t size);

// Remove a memory block pointed by `pointer`
POCKETKALDI_EXPORT
void pk_free(void *pointer);

// The same as strlcpy in FreeBSD
size_t pk_strlcpy(char *dst, const char *src, size_t siz);

// Open the readable file and returns an instance of it. If failed, return NULL
// and status will be set to failed state
POCKETKALDI_EXPORT
pk_readable_t *pk_readable_open(const char *filename, pk_status_t *status);

// Close the readable file
POCKETKALDI_EXPORT
void pk_readable_close(pk_readable_t *self);

// Reads n bytes into buffer from file and store into buffer. buffer should have
// at least n bytes.
POCKETKALDI_EXPORT
void pk_readable_read(
    pk_readable_t *self,
    char *buffer,
    int n,
    pk_status_t *status);

// Reads a line from file and store it into buffer. buffer should have at least
// buffer_size bytes. Returns true if read successfully, otherwise returns
// false. If failure is caused by reaching end-of-file, status->ok == true,
// otherwise status->ok == false 
POCKETKALDI_EXPORT
bool pk_readable_readline(
    pk_readable_t *self,
    char *buffer,
    int buffer_size,
    pk_status_t *status);

// Reads a int32 value from file. 
POCKETKALDI_EXPORT
int32_t pk_readable_readint32(pk_readable_t *self, pk_status_t *status);

// Reads a float value from file. 
POCKETKALDI_EXPORT
float pk_readable_readfloat(pk_readable_t *self, pk_status_t *status);

// Reads bytes from file into byte buffer. The number of bytes to read is
// specified by the size of buffer. And it also reset the current_position of 
// the byte buffer
POCKETKALDI_EXPORT
void pk_readable_readbuffer(
    pk_readable_t *self,
    pk_bytebuffer_t *buffer,
    pk_status_t *status);

// Read a pocketkaldi section head from file, check if it is the samew as
// expected_section. Then return the size of this section. If failed or
// different with expected_section status will be set to fail
POCKETKALDI_EXPORT
int pk_readable_readsectionhead(
    pk_readable_t *self,
    const char *expected_section,
    pk_status_t *status);

// Initialize the bytebuffer 
POCKETKALDI_EXPORT
void pk_bytebuffer_init(pk_bytebuffer_t *self);

// Reset the bytebuffer, allocate `size` bytes
POCKETKALDI_EXPORT
void pk_bytebuffer_reset(pk_bytebuffer_t *self, int64_t size);

// Destroy the bytebuffer  
POCKETKALDI_EXPORT
void pk_bytebuffer_destroy(pk_bytebuffer_t *self);

// Read n bytes from array buffer and store them into `buffer`. And `buffer` 
// should have at least n bytes
POCKETKALDI_EXPORT
void pk_bytebuffer_readbytes(pk_bytebuffer_t *self, char *buffer, int64_t n);

// Read an int32 from array buffer and return it
POCKETKALDI_EXPORT
int32_t pk_bytebuffer_readint32(pk_bytebuffer_t *self);

// Read an int32 from array buffer and return it
POCKETKALDI_EXPORT
int64_t pk_bytebuffer_readint64(pk_bytebuffer_t *self);

// Read a float from array buffer and return it
POCKETKALDI_EXPORT
float pk_bytebuffer_readfloat(pk_bytebuffer_t *self);

// The same as std::nth_elemnt: get the n-th element in array. And the element
// at the nth position is the element that would be in that position in a
// sorted sequence
// base_ptr: pointer to the array
// nmemb: number of elements in array
// size: size of each element
// nth: to get the n-th element
// compar: comparator
// thunk: thunk used in comparasion, usually NULL
// Thanks wengxt :)
POCKETKALDI_EXPORT
void pk_introselect_r(
    void *base_ptr,
    size_t nmemb,
    size_t size,
    size_t nth,
    int (*compar)(const void *, const void *, void *),
    void *thunk);

// 
namespace pocketkaldi {
namespace util {

template<typename T,
         typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
inline std::string ToString(const T &val) {
  return std::to_string(val);
}
template<typename T>
inline std::string ToString(const T * const val) {
  return std::to_string(reinterpret_cast<uint64_t>(val));
}
inline std::string ToString(const std::string &val) {
  return val;
}
inline std::string ToString(const char *const &val) {
  return std::string(val);
}

namespace {

inline std::string FormatImpl(const std::string &formatted) {
  return formatted;
}
template<typename T, typename... Args>
inline std::string FormatImpl(
    const std::string &formatted,
    T &&item,
    Args &&...args) {
  size_t pos = formatted.find("{}");
  std::string repl = ToString(item);
  std::string next_formatted = formatted;
  if (pos != std::string::npos) {
    next_formatted.replace(pos, 2, repl);
  }
  return FormatImpl(next_formatted, std::forward<Args>(args)...);
}

}  // namespace

// Format string function, just like Python, it uses '{}' for replacement. For
// example:
//   util::Format("Hello, {}, {}!", "World", "2233");
template<typename... Args>
inline std::string Format(const std::string &fmt, Args &&...args) {
  return FormatImpl(fmt, std::forward<Args>(args)...);
}

// Trim string 
std::string Trim(const std::string &str);

// Split string by delim and returns as a vector of strings
std::vector<std::string> Split(
    const std::string &str,
    const std::string &delim);

// String tolower
std::string Tolower(const std::string &str);

// A wrapper of FILE in stdio.h
class ReadableFile {
 public:
  ReadableFile();
  ~ReadableFile();

  // Return true if end-of-file reached
  bool Eof() const;

  // Open a file for read. If success, status->ok() will be true. Otherwise,
  // status->ok() == false
  void Open(const std::string &filename, Status *status);

  // Read a line from file. If success, the line will be stored in `line` and
  // status->ok() will be true and return true. Otherwise, if EOF reached
  // first time, return false but status->ok() will still be true.
  // If other error occured, return false anf status->ok() will be false
  bool ReadLine(std::string *line, Status *status);

 private:
  std::string filename_;
  FILE *fd_;
};

// To check if a class have 'previous' field
template <typename T>
struct has_previous {
  template<typename C> static int8_t check(decltype(&C::previous)) ;
  template<typename C> static int16_t check(...);    
  enum {
    value = (sizeof(check<T>(0)) == sizeof(int8_t))
  };
};

}  // namespace util
}  // namespace pocketkaldi

#endif  // POCKETKALDI_UTIL_H_
