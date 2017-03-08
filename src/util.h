// Created at 2016-11-24

#ifndef POCKETKALDI_UTIL_H_
#define POCKETKALDI_UTIL_H_

#include <stdbool.h>
#include "pocketkaldi.h"

#define PK_STATUS_MSGMAX 256

// Error types for status
#define PK_STATUS_SUCCESS 0
#define PK_STATUS_IOERROR 1
#define PK_STATUS_CORRUPTED 2

typedef struct {
  bool ok;
  int errcode;
  char message[PK_STATUS_MSGMAX];
} pk_status_t;

// Initialize the status set to success (ok) state
POCKETKALDI_EXPORT
void pk_status_init(pk_status_t *status);

// Set status to failed state with message
POCKETKALDI_EXPORT
void pk_status_fail(pk_status_t *status, int errcode, const char *msg);

#endif  // POCKETKALDI_UTIL_H_
