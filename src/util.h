// Created at 2016-11-24

#ifndef POCKETKALDI_UTIL_H_
#define POCKETKALDI_UTIL_H_

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include "pocketkaldi.h"

#define PK_STATUS_MSGMAX 256

// Error types for status
#define PK_STATUS_STSUCCESS 0
#define PK_STATUS_STIOERROR 1
#define PK_STATUS_STCORRUPTED 2

#define PK_WARN fprintf
#define PK_UNUSED(x) (void)(x)

#define PK_PATHMAX 1024

typedef struct pk_status_t {
  bool ok;
  int errcode;
  char message[PK_STATUS_MSGMAX];
} pk_status_t;

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



#endif  // POCKETKALDI_UTIL_H_
