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

namespace pocketkaldi {

// Computes the raw CMVN stats for frame
void CMVN::ComputeStats(int frame, Vector<float> *stats) {
  // Currently we only support compute frame by frame. Random frame computation
  // will be implemented when needed (like ring buffer based cache in Kaldi)
  assert(cached_frame_ == frame - 1);

  pk_vector_t c_feats = pk_matrix_getcol(raw_feats_, frame);
  SubVector<float> feats(c_feats.data, c_feats.dim);

  // We do the computation in double
  Vector<double> stats_dbl(PK_FBANK_DIM + 1);
  if (cached_frame_ >= 0) {
    // if cached_frame_ < 0, stats_dbl is filled with 0
    assert(cached_stats_.Dim() == PK_FBANK_DIM + 1);
    stats_dbl.CopyFromVec(cached_stats_);
  }

  // Add the stats of current frame
  assert(feats.Dim() == PK_FBANK_DIM);
  stats_dbl.Range(0, PK_FBANK_DIM).AddVec(1.0f, feats);
  stats_dbl(PK_FBANK_DIM) += 1.0;

  // Subtract the stats of previous (slideing window)
  int prev_frame = frame - PK_ONLINECMVN_WINDOW;
  if (prev_frame >= 0) {
    pk_vector_t c_prev_feats = pk_matrix_getcol(raw_feats_, prev_frame);
    SubVector<float> prev_feats(c_prev_feats.data, c_prev_feats.dim);
    stats_dbl.Range(0, PK_FBANK_DIM).AddVec(-1.0f, prev_feats);
    stats_dbl(PK_FBANK_DIM) -= 1.0;
  }

  // Store to stats and stats cache
  stats->Resize(PK_FBANK_DIM + 1, Vector<float>::kUndefined);
  stats->CopyFromVec(stats_dbl);
  cached_frame_ = frame;
  cached_stats_.Resize(stats->Dim());
  cached_stats_.CopyFromVec(*stats);
}

void CMVN::SmoothStats(Vector<float> *stats) {
  assert(stats->Dim() == PK_FBANK_DIM + 1);
  assert(global_stats_.Dim() == PK_FBANK_DIM + 1);
  double count = (*stats)(PK_FBANK_DIM);
  assert(count <= PK_ONLINECMVN_WINDOW &&
         "something went wrong in compute_stats_nextframe()");

  // Smoothing is not necessary when frame count is enough
  if (count >= PK_ONLINECMVN_WINDOW) return;

  double count_from_global = PK_ONLINECMVN_WINDOW - count;
  double global_count = global_stats_(PK_FBANK_DIM);
  assert(global_count > 0);
  if (count_from_global > PK_ONLINECMVN_GLOBALFRAMES) {
    count_from_global = PK_ONLINECMVN_GLOBALFRAMES;
  }

  double scalar = count_from_global / global_count;
  stats->AddVec(scalar, global_stats_);
}

void CMVN::Apply(const VectorBase<float> &stats, VectorBase<float> *feats) {
  assert(feats->Dim() == stats.Dim() - 1 && feats->Dim() == PK_FBANK_DIM);
  double count = stats(PK_FBANK_DIM);
  assert(count > 0);

  float scale = 1 / count;
  feats->AddVec(-scale, stats.Range(0, feats->Dim()));
}

void CMVN::GetFrame(int frame, pk_vector_t *feats) {
  Vector<float> stats;

  // Calculate and smooth stats
  ComputeStats(frame, &stats);
  SmoothStats(&stats);

  // Apply CMVN with stats
  pk_vector_t raw_frame = pk_matrix_getcol(raw_feats_, frame);
  pk_vector_copy(feats, &raw_frame);
  SubVector<float> frame_feats(feats->data, feats->dim);
  Apply(stats, &frame_feats);
}


CMVN::CMVN(const pk_vector_t *global_stats, const pk_matrix_t *raw_feats) {
  global_stats_.Resize(global_stats->dim);
  global_stats_.CopyFromVec(SubVector<float>(
      global_stats->data,
      global_stats->dim));
  raw_feats_ = raw_feats;
  cached_frame_ = -1;
}

CMVN::~CMVN() {
  raw_feats_ = nullptr;
  cached_frame_ = 0;
}

}  // namespace pocketkaldi
