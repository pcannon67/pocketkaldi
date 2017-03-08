// Created at 2016-11-24

#include "fst.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

static void pk_fst_init(pk_fst_t *self) {
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
    const char *filename,
    pk_status_t *status) {
  pk_status_init(status);
  pk_fst_init(self);

  FILE *fd = fopen(filename, "rb");
  if (fd == NULL) pk_status_fail(status, PK_STATUS_IOERROR, filename);

  int64_t filesize = 0;
  if (status->ok) {
    fseek(fd, 0L, SEEK_END);
    filesize = ftell(fd);
    fseek(fd, 0L, SEEK_SET);
    
    // Fst file should have at least 16 bytes (4 int32_t)
    if (filesize <= 16) pk_status_fail(status, PK_STATUS_CORRUPTED, filename);
  }

  char *data = NULL;
  if (status->ok) {
    data = (char *)malloc(filesize);
    int blocks_read = fread(data, filesize, 1, fd);
    if (blocks_read != 1) {
      pk_status_fail(status, PK_STATUS_IOERROR, filename);
    }
  }
  char *buffer_ptr = data;

  // Check magic number
  if (status->ok) {
    int32_t magic_number = *((int32_t *)buffer_ptr);
    buffer_ptr += 4;
    if (magic_number != 0x3323) {
      pk_status_fail(status, PK_STATUS_CORRUPTED, filename);
    }
  }

  // Read metadata and check filesize
  int state_number = 0;
  int arc_number = 0;
  int start_state = 0;
  if (status->ok) {
    state_number = *((int32_t *)buffer_ptr);
    buffer_ptr += sizeof(int32_t);
    arc_number = *((int32_t *)buffer_ptr);
    buffer_ptr += sizeof(int32_t);
    start_state = *((int32_t *)buffer_ptr);
    buffer_ptr += sizeof(int32_t);

    // Now we can calculate the expected filesize
    int64_t expected_size = 16;  // Base bytes
    expected_size += state_number * sizeof(float);  // Final
    expected_size += state_number * sizeof(int32_t);  // Arc index
    expected_size += arc_number * sizeof(pk_fst_arc_t);  // Arc

    self->state_number = state_number;
    self->arc_number = arc_number;
    self->start_state = start_state;

    if (expected_size != filesize) {
      pk_status_fail(status, PK_STATUS_CORRUPTED, filename);
    }
  }

  // Read other parts
  if (status->ok) {
    self->final = (float *)malloc(state_number * sizeof(float));
    for (int idx = 0; idx < state_number; ++idx) {
      self->final[idx] = *((float *)buffer_ptr);
      buffer_ptr += sizeof(float);
    }

    self->arc_index = (int32_t *)malloc(state_number * sizeof(int32_t));
    for (int idx = 0; idx < state_number; ++idx) {
      self->arc_index[idx] = *((int32_t *)buffer_ptr);
      buffer_ptr += sizeof(int32_t);
    }

    self->arc = (pk_fst_arc_t *)malloc(arc_number * sizeof(pk_fst_arc_t));
    for (int idx = 0; idx < arc_number; ++idx) {
      self->arc[idx].next_state = *((int32_t *)buffer_ptr);
      buffer_ptr += sizeof(int32_t);
      self->arc[idx].input_label = *((int32_t *)buffer_ptr);
      buffer_ptr += sizeof(int32_t);
      self->arc[idx].output_label = *((int32_t *)buffer_ptr);
      buffer_ptr += sizeof(int32_t);
      self->arc[idx].weight = *((float *)buffer_ptr);
      buffer_ptr += sizeof(float);
    }

    if (buffer_ptr - data != filesize) {
      pk_status_fail(status, PK_STATUS_CORRUPTED, filename);
    }
  }

  if (fd != NULL) fclose(fd);
  if (!status->ok) {
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
