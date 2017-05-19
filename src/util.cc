// Created at 2016-11-24

#include "util.h"

#include <assert.h>
#include <ctype.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <array>
#include <algorithm>

void pk_status_init(pk_status_t *status) {
  status->ok = true;
  status->errcode = 0;
  status->message[0] = '\0';
}

void pk_status_fail(pk_status_t *status, int errcode, const char *fmsg, ...) {
  status->ok = false;
  status->errcode = errcode;
  const char *error_prefix = NULL;
  char unknown_error[128];
  char msg[128];

  // Gets error message from fmsg
  va_list args;
  va_start(args, fmsg);
  vsnprintf(msg, sizeof(msg), fmsg, args);
  va_end(args);
  
  // Gets string representation of error code
  switch (errcode) {
    case PK_STATUS_STIOERROR:
      error_prefix = "IOError";
      break;
    case PK_STATUS_STCORRUPTED:
      error_prefix = "Corrupted";
      break;
    default:
      snprintf(
          unknown_error,
          sizeof(unknown_error),
          "UnknownError(%d)",
          errcode);
      error_prefix = unknown_error;
  }

  snprintf(
      status->message,
      PK_STATUS_MSGMAX,
      "%s: %s",
      error_prefix,
      msg);
}

void *pk_alloc(size_t size) {
  return malloc(size);
}

void *pk_realloc(void *ptr, size_t size) {
  return realloc(ptr, size);
}

void pk_free(void *pointer) {
  free(pointer);
}

pk_readable_t *pk_readable_open(const char *filename, pk_status_t *status) {
  pk_readable_t *self = (pk_readable_t *)malloc(sizeof(pk_readable_t));
  pk_strlcpy(self->filename, filename, PK_PATHMAX);
  self->fd = fopen(filename, "rb");
  if (self->fd == NULL) {
    PK_STATUS_IOERROR(status, "unable to open: %s", filename);
    free(self);
    self = NULL;
  } else {
    // Get file size
    fseek(self->fd, 0, SEEK_END);
    self->filesize = ftell(self->fd);
    fseek(self->fd, 0, SEEK_SET);
  }

  return self;
}

void pk_readable_close(pk_readable_t *self) {
  if (self->fd != NULL) fclose(self->fd);
  self->fd = NULL;
  self->filesize = 0;
  self->filename[0] = '\0';
  free(self);
}

void pk_readable_read(
    pk_readable_t *self,
    char *buffer,
    int n,
    pk_status_t *status) {
  int bytes_read = fread(buffer, 1, n, self->fd);
  if (bytes_read != n) {
    PK_STATUS_IOERROR(status, "read failed: %s", self->filename);
  }
}

int32_t pk_readable_readint32(pk_readable_t *self, pk_status_t *status) {
  int32_t val = 0;;
  pk_readable_read(self, (char *)&val, sizeof(int32_t), status);
  return val;
}

float pk_readable_readfloat(pk_readable_t *self, pk_status_t *status) {
  float val = 0.0f;
  pk_readable_read(self, (char *)&val, sizeof(float), status);
  return val;
}

void pk_readable_readbuffer(
    pk_readable_t *self,
    pk_bytebuffer_t *buffer,
    pk_status_t *status) {
  int bytes_to_read = buffer->size;
  assert(bytes_to_read != 0 && "read into a empty buffer");

  buffer->current_position = 0;
  pk_readable_read(self, buffer->buffer, bytes_to_read, status);
}

bool pk_readable_readline(
    pk_readable_t *self,
    char *buffer,
    int buffer_size,
    pk_status_t *status) {
  // Failed if it already reached EOF
  if (feof(self->fd)) {
    PK_STATUS_IOERROR(status, "EOF already reached: %s", self->filename);
    return false;
  }

  // Readline
  char *s = fgets(buffer, buffer_size, self->fd);
  if (s == NULL) {
    if (feof(self->fd)) {
      // First time that reached EOF
      return false;
    } else {
      PK_STATUS_IOERROR(status, "%s", self->filename);
      return false;
    }
  }

  // Trim the tailing '\r' or '\n'
  char *tail = buffer + strlen(buffer) - 1;
  while ((*tail == '\r' || *tail == '\n') && tail >= buffer) {
    *tail = '\0';
    --tail;
  }
}

int pk_readable_readsectionhead(
    pk_readable_t *self,
    const char *expected_section,
    pk_status_t *status) {
  assert(strlen(expected_section) == 4 && "section_name should be 4 chars");

  // Read section header
  pk_bytebuffer_t section_header;
  pk_bytebuffer_init(&section_header);
  pk_bytebuffer_reset(&section_header, 8);
  pk_readable_readbuffer(self, &section_header, status);

  // Check section name
  if (status->ok) {
    char section_name[5];
    pk_bytebuffer_readbytes(&section_header, section_name, 4);
    section_name[4] = '\0';
    if (strcmp(section_name, expected_section) != 0) {
      PK_STATUS_CORRUPTED(
          status,
          "section_name == %s expected, but %s found (%s)",
          expected_section,
          section_name,
          self->filename);
    }
  }

  int32_t section_size = 0;
  if (status->ok) {
    section_size = pk_bytebuffer_readint32(&section_header);
  }

  pk_bytebuffer_destroy(&section_header);
  return section_size;
}


void pk_bytebuffer_init(pk_bytebuffer_t *self) {
  self->buffer = NULL;
  self->size = 0;
  self->current_position = 0;
}

void pk_bytebuffer_reset(pk_bytebuffer_t *self, int64_t size) {
  free(self->buffer);
  self->buffer = (char *)malloc(sizeof(char) * size);
  self->size = size;
  self->current_position = 0;
}

void pk_bytebuffer_destroy(pk_bytebuffer_t *self) {
  free(self->buffer);
  self->buffer = NULL;
  self->size = 0;
  self->current_position = 0;
}

void pk_bytebuffer_readbytes(
    pk_bytebuffer_t *self,
    char *buffer,
    int64_t n) {
  assert(self->current_position + n <= self->size &&
         "pk_bytebuffer: out of boundary");
  for (int i = 0; i < n; ++i) {
    buffer[i] = self->buffer[self->current_position + i];
  }
  self->current_position += n;
}

int32_t pk_bytebuffer_readint32(pk_bytebuffer_t *self) {
  int32_t val;
  pk_bytebuffer_readbytes(self, (char *)&val, sizeof(int32_t));
  return val;
}

int64_t pk_bytebuffer_readint64(pk_bytebuffer_t *self) {
  int64_t val;
  pk_bytebuffer_readbytes(self, (char *)&val, sizeof(int64_t));
  return val;
}

float pk_bytebuffer_readfloat(pk_bytebuffer_t *self) {
  float val;
  pk_bytebuffer_readbytes(self, (char *)&val, sizeof(float));
  return val;
}

namespace pocketkaldi {
namespace util {

std::string Trim(const std::string &str) {
  std::string::const_iterator begin = str.cbegin();
  std::string::const_iterator end = str.cend() - 1;

  while (begin < str.cend() && isspace(*begin)) ++begin;
  while (end > begin && isspace(*end)) --end;
  return std::string(begin, end + 1);
}

std::vector<std::string> Split(
    const std::string &str,
    const std::string &delim) {
  std::vector<std::string> fields;
  int start = 0;
  int pos = 0;
  while ((pos = str.find(delim, start)) != std::string::npos) {
    fields.emplace_back(str.cbegin() + start, str.cbegin() + pos);
    start = pos + delim.size();
  }
  if (str.cbegin() + start < str.cend()) {
    fields.emplace_back(str.cbegin() + start, str.cend());
  }

  return fields;
}

std::string Tolower(const std::string &str) {
  std::string lower(str.begin(), str.end());
  std::transform(lower.begin(), lower.end(), lower.begin(), tolower);
  return lower;
}


ReadableFile::ReadableFile(): fd_(nullptr) {
}

ReadableFile::~ReadableFile() {
  if (fd_ != nullptr) fclose(fd_);
  fd_ = NULL;
}

Status ReadableFile::Open(const std::string &filename) {
  filename_ = filename;
  fd_ = fopen(filename.c_str(), "rb");
  if (fd_ == NULL) {
    return Status::IOError(util::Format("Unable to open {}", filename));
  } else {
    return Status::OK();
  }
}

bool ReadableFile::ReadLine(std::string *line, Status *status) {
  // Failed if it already reached EOF
  if (feof(fd_)) {
    *status = Status::IOError(
        util::Format("EOF already reached: {}", filename_));
    return false;
  }

  // Readline
  std::array<char, 4096> chunk;
  char *s = fgets(chunk.data(), chunk.size(), fd_);
  if (s == NULL) {
    if (feof(fd_)) {
      // First time that reached EOF
      return false;
    } else {
      *status = Status::IOError(filename_);
      return false;
    }
  }

  // Trim the tailing '\r' or '\n'
  *line = s;
  while (line->empty() == false &&
         (line->back() == '\r' || line->back() == '\n')) {
    line->pop_back();
  }
  return true;
}

Status ReadableFile::Read(void *ptr, int size) {
  if (1 != fread(ptr, size, 1, fd_)) {
    return Status::IOError(util::Format("failed to read: {}", filename_));
  } else {
    return Status::OK();
  }
}

bool ReadableFile::Eof() const {
  return feof(fd_) != 0;
}

}  // namespace util
}  // namespace pocketkaldi
