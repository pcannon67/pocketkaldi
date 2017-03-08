// Created at 2016-11-24

#include "util.h"

#include <string.h>
#include <stdio.h>

void pk_status_init(pk_status_t *status) {
  status->ok = true;
  status->errcode = 0;
  status->message[0] = '\0';
}

void pk_status_fail(pk_status_t *status, int errcode, const char *msg) {
  status->ok = false;
  status->errcode = errcode;
  const char *error_prefix = NULL;
  char unknown_error[128];
  
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
