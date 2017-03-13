// Created at 2017-03-09

#ifndef POCKETKALDI_TRANSMODEL_H_
#define POCKETKALDI_TRANSMODEL_H_

#include <stdint.h>
#include "util.h"

#define PK_TRANSITION_SECTION "TRM0"

// The transtion model used in Kaldi. And it only implements the functions
// useful for decoding. Full implementation could be found in
// kaldi::TransitionModel
typedef struct pk_transition_t {
  int32_t num_transition_ids;
  int32_t num_pdfs;
  int32_t *id2pdf;
} pk_transition_t;

// Read the transition model from fd. And when failed, set status to failed
// state
POCKETKALDI_EXPORT
void pk_transition_init(
    pk_transition_t *self,
    pk_readable_t *fd,
    pk_status_t *status);

// Map the transtion-id to pdf-id
POCKETKALDI_EXPORT
int pk_transition_tid2pdf(const pk_transition_t *self, int transtion_id);

// Destroy the transtion model
POCKETKALDI_EXPORT
void pk_transition_destroy(pk_transition_t *self);

#endif