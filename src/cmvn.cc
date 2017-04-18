// online_cmvn.c (Original file in Kaldi: feat/online-feature.cc)

// Copyright    2013  Johns Hopkins University (author: Daniel Povey)
//              2014  Yanqing Sun, Junjie Wang,
//                    Daniel Povey, Korbinian Riedhammer
//              2017  Xiaoyang Chen

// See ../../COPYING for clarification regarding multiple authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
// THIS CODE IS PROVIDED *AS IS* BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY IMPLIED
// WARRANTIES OR CONDITIONS OF TITLE, FITNESS FOR A PARTICULAR PURPOSE,
// MERCHANTABLITY OR NON-INFRINGEMENT.
// See the Apache 2 License for the specific language governing permissions and
// limitations under the License.

// Created at 2017-03-06

#include "cmvn.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "fbank.h"

// Computes the raw CMVN stats for frame
static void compute_stats(pk_cmvn_t *self, int frame, pk_vector_t *stats) {
  // Currently we only support compute frame by frame. Random frame computation
  // will be implemented when needed (like ring buffer based cache in Kaldi)
  assert(self->cached_frame == frame - 1);

  pk_vector_t feats = pk_matrix_getcol(self->raw_feats, frame);

  // We do the computation in double
  double stats_dbl[PK_FBANK_DIM + 1];
  if (self->cached_frame < 0) {
    for (int i = 0; i < PK_FBANK_DIM + 1; ++i) stats_dbl[i] = 0.0;
  } else {
    assert(self->cached_stats.dim == PK_FBANK_DIM + 1);
    for (int i = 0; i < PK_FBANK_DIM + 1; ++i) {
      stats_dbl[i] = self->cached_stats.data[i];
    }
  }

  // Add the stats of current frame
  assert(feats.dim == PK_FBANK_DIM);
  for (int i = 0; i < feats.dim; ++i) {
    stats_dbl[i] += feats.data[i];
  }
  stats_dbl[PK_FBANK_DIM] += 1.0;

  // Subtract the stats of previous (slideing window)
  int prev_frame = frame - PK_ONLINECMVN_WINDOW;
  if (prev_frame >= 0) {
    pk_vector_t prev_feats = pk_matrix_getcol(self->raw_feats, prev_frame);
    for (int i = 0; i < prev_feats.dim; ++i) {
      stats_dbl[i] -= prev_feats.data[i];
    }
    stats_dbl[PK_FBANK_DIM] -= 1.0;
  }

  // Store to stats and stats cache
  pk_vector_resize(stats, PK_FBANK_DIM + 1);
  for (int i = 0; i < stats->dim; ++i) {
    stats->data[i] = (float)stats_dbl[i];
  }
  self->cached_frame = frame;
  pk_vector_copy(&self->cached_stats, stats);
}

// Smooth the CMVN stats "stats", by possibly adding some stats from
// global_stats
static void smooth_stats(pk_cmvn_t *self, pk_vector_t *stats) {
  assert(stats->dim == PK_FBANK_DIM + 1);
  assert(self->global_stats->dim == PK_FBANK_DIM + 1);
  double count = stats->data[PK_FBANK_DIM];
  assert(count <= PK_ONLINECMVN_WINDOW &&
         "something went wrong in compute_stats_nextframe()");

  // Smoothing is not necessary when frame count is enough
  if (count >= PK_ONLINECMVN_WINDOW) return;

  double count_from_global = PK_ONLINECMVN_WINDOW - count;
  double global_count = self->global_stats->data[PK_FBANK_DIM];
  assert(global_count > 0);
  if (count_from_global > PK_ONLINECMVN_GLOBALFRAMES) {
    count_from_global = PK_ONLINECMVN_GLOBALFRAMES;
  }

  double scalar = count_from_global / global_count;
  for (int i = 0; i < stats->dim; ++i) {
    double stats_val = stats->data[i];
    double global_stats_val = self->global_stats->data[i];
    stats->data[i] = (float)(stats_val + scalar * global_stats_val);
  }
}

void pk_cmvn_apply(const pk_vector_t *stats, pk_vector_t *feats) {
  assert(feats->dim == stats->dim - 1 && feats->dim == PK_FBANK_DIM);
  double count = stats->data[PK_FBANK_DIM];
  assert(count > 0);
  
  for (int d = 0; d < feats->dim; ++d) {
    double mean = stats->data[d] / count;
    feats->data[d] -= mean;
  }
}

void pk_cmvn_getframe(pk_cmvn_t *self, int frame, pk_vector_t *feats) {
  assert(feats->dim == PK_FBANK_DIM && "pk_cmvn_getframe: feats dim mismatch");

  pk_vector_t stats;
  pk_vector_init(&stats, 0, NAN);

  // Calculate and smooth stats
  compute_stats(self, frame, &stats);
  smooth_stats(self, &stats);

  // Apply CMVN with stats
  pk_vector_t raw_frame = pk_matrix_getcol(self->raw_feats, frame);
  pk_vector_copy(feats, &raw_frame);
  pk_cmvn_apply(&stats, feats);

  pk_vector_destroy(&stats);
}

void pk_cmvn_init(
    pk_cmvn_t *self,
    const pk_vector_t *global_stats,
    const pk_matrix_t *raw_feats) {
  pk_vector_init(&self->cached_stats, 0, NAN);
  self->global_stats = global_stats;
  self->raw_feats = raw_feats;
  self->cached_frame = -1;
}

void pk_cmvn_destroy(pk_cmvn_t *self) {
  pk_vector_destroy(&self->cached_stats);
  self->global_stats = NULL;
  self->raw_feats = NULL;
  self->cached_frame = 0;
}
