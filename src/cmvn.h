// Created at 2017-03-06

#ifndef POCKETKALDI_ONLINE_CMVN_H_
#define POCKETKALDI_ONLINE_CMVN_H_

#include <stdint.h>
#include "matrix.h"
#include "util.h"

#define PK_ONLINECMVN_WINDOW 600
#define PK_ONLINECMVN_GLOBALFRAMES 200

// Stores the intermediate state in computing CMVN of an utterance
// Currently only mean computation is supported
typedef struct pk_cmvn_t {
  const pk_matrix_t *raw_feats;
  const pk_vector_t *global_stats;
  pk_vector_t cached_stats;
  int cached_frame;
} pk_cmvn_t;

// Apply CMVN to frame
POCKETKALDI_EXPORT
void pk_cmvn_apply(const pk_vector_t *stats, pk_vector_t *feats);

// Initialize online cmvn. global_stats is the accumulated CMVN stats from
// training data. raw_feats is raw feature matrix. (global_stats and raw_feats
// are borrowed)
POCKETKALDI_EXPORT
void pk_cmvn_init(
    pk_cmvn_t *self,
    const pk_vector_t *global_stats,
    const pk_matrix_t *raw_feats);

// Destroy online cmvn
POCKETKALDI_EXPORT
void pk_cmvn_destroy(pk_cmvn_t *self);

// Gets the fests of frame after appling online CMVN.
POCKETKALDI_EXPORT
void pk_cmvn_getframe(pk_cmvn_t *self, int frame, pk_vector_t *feats);

#endif  // POCKETKALDI_ONLINE_CMVN_H_

