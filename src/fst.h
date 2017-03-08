// Created at 2016-11-24

#ifndef POCKETKALDI_FST_H_
#define POCKETKALDI_FST_H_

#include <stdint.h>
#include "util.h"
#include "pocketkaldi.h"

typedef struct {
  int32_t next_state;
  int32_t input_label;
  int32_t output_label;
  float weight;
} pk_fst_arc_t;

typedef struct {
  int32_t start_state;
  int32_t state_number;
  int32_t arc_number;

  pk_fst_arc_t *arc;
  float *final;
  int32_t *arc_index;
} pk_fst_t;

typedef struct {
  int base_idx;
  int offset;
  int total;
  const pk_fst_t *fst;
} pk_fst_iter_t;

// Destroy it and free all memory
POCKETKALDI_EXPORT
void pk_fst_destroy(pk_fst_t *self);

// Read fst from file `filename`. On success, status->ok is true. On failed,
// status->ok is false
POCKETKALDI_EXPORT
void pk_fst_read(
    pk_fst_t *self,
    const char *filename,
    pk_status_t *status);

// Itereate all outcoming arcs of state. It will initialize an `pk_fst_iter_t`
// object `iter`. And then arcs could be iterated by following code:
// 
// pk_fst_arc_t *arc;
// while ((arc = pk_fst_iter_next(iter)) != NULL) {
//   ...
// }
POCKETKALDI_EXPORT
void pk_fst_iterate_arc(
    const pk_fst_t *self,
    int state,
    pk_fst_iter_t *iter);

// Iterate arcs from an iterator. If next arc exists, retrun it and move the
// iterator forward, else return NULL
POCKETKALDI_EXPORT
const pk_fst_arc_t *pk_fst_iter_next(pk_fst_iter_t *self);

// Gets the start state of fst
POCKETKALDI_EXPORT
int pk_fst_startstate(const pk_fst_t *self);

// Get the final score of state. If the state is non-terminal, returns 0
POCKETKALDI_EXPORT
float pk_fst_final(const pk_fst_t *self, int state);

#endif  // POCKETKALDI_FST_H_