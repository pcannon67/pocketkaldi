// Created at 2017-03-22

#ifndef POCKETKALDI_AM_H_
#define POCKETKALDI_AM_H_

#include "nnet.h"
#include "util.h"

#define PK_AM_SECTION "AM~0"

// AM struct
typedef struct pk_am_t {
  pk_nnet_t nnet;
  pk_vector_t log_prior;
  int left_context;
  int right_context;
} pk_am_t;

// Initialize the AM
POCKETKALDI_EXPORT
void pk_am_init(pk_am_t *self);

// Read am from fd. When any error occured, status->ok == false
POCKETKALDI_EXPORT
void pk_am_read(pk_am_t *self, pk_readable_t *fd, pk_status_t *status);

// Destpry the AM
POCKETKALDI_EXPORT
void pk_am_destroy(pk_am_t *self);

// Compute the log-likelihood of the frame indexed by frame_idx in frames
// And output as a vector loglikelihood
POCKETKALDI_EXPORT
void pk_am_compute(
    const pk_am_t *self,
    const pk_matrix_t *frames,
    int frame_idx,
    pk_vector_t *loglikelihood);

#endif
