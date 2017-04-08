// Create at 2017-02-23

#ifndef POCKETKALDI_DECODABLE_H_
#define POCKETKALDI_DECODABLE_H_

#include "decodable.h"

#include <assert.h>
#include "am.h"
#include "transition.h"
#include "matrix.h"

// Used in decoder
typedef struct pk_decodable_t {
  pk_transition_t *trans_model;
  pk_matrix_t log_prob;
} pk_decodable_t;

// Initialize the decodable struct
POCKETKALDI_EXPORT
void pk_decodable_init(
    pk_decodable_t *self,
    pk_am_t *am,
    pk_transition_t *trans_model,
    float prob_scale,
    const pk_matrix_t *feats);

// Destroy the decodable struct
POCKETKALDI_EXPORT
void pk_decodable_destroy(pk_decodable_t *self);

// Gets the log-likelihood of frame and transition-id
POCKETKALDI_EXPORT
float pk_decodable_loglikelihood(
    pk_decodable_t *self,
    int frame,
    int trans_id);

// Return true if frame is the last frame 
POCKETKALDI_EXPORT
bool pk_decodable_islastframe(pk_decodable_t *self, int frame);

#endif  // POCKETKALDI_DECODABLE_H_
