// Created at 2016-11-24

#include "fst.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

void pk_fst_init(pk_fst_t *self) {
  self->start_state = 0;
  self->state_number = 0;
  self->arc_number = 0;

  self->arc = NULL;
  self->arc_index = NULL;
  self->final = NULL;
}

void pk_fst_destroy(pk_fst_t *self) {
  self->start_state = 0;
  self->state_number = 0;
  self->arc_number = 0;

  free(self->arc);
  self->arc = NULL;
  free(self->arc_index);
  self->arc_index = NULL;
  free(self->final);
  self->final = NULL;
}

void pk_fst_read(
    pk_fst_t *self,
    pk_readable_t *fd,
    pk_status_t *status) {
  // Check magic number
  int32_t magic_number = pk_readable_readint32(fd, status);
  if (!status->ok) goto pk_fst_read_failed;
  if (magic_number != 0x3323) {
    PK_STATUS_CORRUPTED(status, "%s", fd->filename);
    goto pk_fst_read_failed;
  }

  // Read metadata and check filesize
  self->state_number = pk_readable_readint32(fd, status);
  if (!status->ok) goto pk_fst_read_failed;
  self->arc_number = pk_readable_readint32(fd, status);
  if (!status->ok) goto pk_fst_read_failed;
  self->start_state = pk_readable_readint32(fd, status);
  if (!status->ok) goto pk_fst_read_failed;

  // Final weight
  self->final = (float *)malloc(self->state_number * sizeof(float));
  for (int idx = 0; idx < self->state_number; ++idx) {
    self->final[idx] = pk_readable_readfloat(fd, status);
    if (!status->ok) goto pk_fst_read_failed;
  }

  self->arc_index = (int32_t *)malloc(self->state_number * sizeof(int32_t));
  for (int idx = 0; idx < self->state_number; ++idx) {
    self->arc_index[idx] = pk_readable_readint32(fd, status);
    if (!status->ok) goto pk_fst_read_failed;
  }

  self->arc = (pk_fst_arc_t *)malloc(self->arc_number * sizeof(pk_fst_arc_t));
  for (int idx = 0; idx < self->arc_number; ++idx) {
    self->arc[idx].next_state = pk_readable_readint32(fd, status);
    if (!status->ok) goto pk_fst_read_failed;
    self->arc[idx].input_label = pk_readable_readint32(fd, status);
    if (!status->ok) goto pk_fst_read_failed;
    self->arc[idx].output_label = pk_readable_readint32(fd, status);
    if (!status->ok) goto pk_fst_read_failed;
    self->arc[idx].weight = pk_readable_readfloat(fd, status);
    if (!status->ok) goto pk_fst_read_failed;
  }

  if (false) {
pk_fst_read_failed:
    pk_fst_destroy(self);
  }
}

// Calcuate the number of outcoming arcs for state
static int calc_state_outcomming_arc_number(const pk_fst_t *fst, int state) {
  int32_t *arcidx = fst->arc_index;
  int state_number = fst->state_number;
  int state_idx = arcidx[state];

  if (state_idx < 0) return 0;
  int count;
  int next_state = -1;

  // Find the next state that have outcoming arcs
  for (int check_state = state + 1; check_state < state_number; ++check_state) {
    if (arcidx[check_state] > 0) {
      next_state = check_state;
      break;
    }
  }
  int next_idx = next_state >= 0 ? arcidx[next_state] : fst->arc_number;
  return next_idx - state_idx;
}

const pk_fst_arc_t *pk_fst_iter_next(pk_fst_iter_t *self) {
  if (self->offset < self->total) {
    pk_fst_arc_t *arc = &(self->fst->arc[self->base_idx + self->offset]);
    ++self->offset;
    return arc;
  } else {
    return NULL;
  }
}

void pk_fst_iterate_arc(
    const pk_fst_t *self,
    int state,
    pk_fst_iter_t *iter) {
  assert(state < self->state_number && "State index out-of-boundary");
  iter->fst = self;
  iter->total = calc_state_outcomming_arc_number(self, state);
  assert(iter->total >= 0);
  iter->offset = 0;
  iter->base_idx = self->arc_index[state];
}

int pk_fst_startstate(const pk_fst_t *self) {
  return self->start_state;
}

float pk_fst_final(const pk_fst_t *self, int state) {
  assert(state < self->state_number && "State index out-of-boundary");
  return self->final[state];
}
