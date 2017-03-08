// Created at 2016-11-24

#include "util.h"

#include <stdarg.h>
#include <string.h>
#include <stdio.h>

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
    case PK_STATUS_IOERROR:
      error_prefix = "IOError";
      break;
    case PK_STATUS_CORRUPTED:
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
