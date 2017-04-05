// decoder/simple-decoder.h

// Copyright 2009-2013  Microsoft Corporation;  Lukas Burget;
//                      Saarland University (author: Arnab Ghoshal);
//                      Johns Hopkins University (author: Daniel Povey)

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

#ifndef POCKETKALDI_DECODER_H_
#define POCKETKALDI_DECODER_H_

#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include "hashlist.h"
#include "transition.h"
#include "decodable.h"
#include "vector.h"
#include "fst.h"
#include "am.h"

#define PK_DECODER_BEAMSIZE 30000
#define PK_DECODER_BEAMDELTA 0.5
#define PK_DECODER_MINACTIVE 200

typedef struct {
  int32_t *alignment;
  int32_t alignment_size;
  int32_t *words;
  int32_t size;
  float weight;
} pk_decoder_result_t;

// Initialize the pk_decoder_result_t
POCKETKALDI_EXPORT
void pk_decoder_result_init(pk_decoder_result_t *self);

// Destroy the pk_decoder_result_t.
POCKETKALDI_EXPORT
void pk_decoder_result_destroy(pk_decoder_result_t *best_path);

typedef struct pk_decoder_token_t {
  struct pk_decoder_token_t *previous;
  int input_label;
  int output_label;
  int next_state;
  float arc_weight;
  float cost;
  int32_t ref_count;
} pk_decoder_token_t;


typedef struct pk_decoder_t {
  pk_hashlist_t cur_toks;
  pk_hashlist_t prev_toks;
  pk_hashlist_t tmp_toks;

  // empty_toks is a linklist of empty (freed) tokens. And these tokens a linked
  // by the previous field
  pk_decoder_token_t *empty_toks;

  const pk_fst_t *fst;
  float beam;

  // Keep track of the number of frames decoded in the current file.
  int32_t num_frames_decoded;
} pk_decoder_t;

// Initialize the decoder. It just borrows the fst
POCKETKALDI_EXPORT
void pk_decoder_init(pk_decoder_t *self, const pk_fst_t *fst, float beam);

// Destroy the decoder
POCKETKALDI_EXPORT
void pk_decoder_destroy(pk_decoder_t *self);

// Decode given decoable object
POCKETKALDI_EXPORT
bool pk_decoder_decode(pk_decoder_t *self, pk_decodable_t *decodable);

// Check if it reached the final state
POCKETKALDI_EXPORT
bool pk_decoder_reachedfinal(pk_decoder_t *self);

// Outputs best_path corresponding to the single best path through the lattice.
POCKETKALDI_EXPORT
bool pk_decoder_bestpath(
    pk_decoder_t *self,
    pk_decoder_result_t *best_path,
    bool use_final_probs);

#endif  // POCKETKALDI_DECODER_H_
