// Create at 2017-02-23

#ifndef POCKETKALDI_DECODABLE_H_
#define POCKETKALDI_DECODABLE_H_

#include "decodable.h"

#include <assert.h>
#include "am.h"
#include "transition.h"
#include "matrix.h"

using pocketkaldi::AcousticModel;

// Used in decoder
typedef struct pk_decodable_t {
  pk_matrix_t log_prob;
  AcousticModel *am;
} pk_decodable_t;

// Initialize the decodable struct
POCKETKALDI_EXPORT
void pk_decodable_init(
    pk_decodable_t *self,
    AcousticModel *am,
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
