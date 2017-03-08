// decoder/simple-decoder.cc

// Copyright 2009-2011 Microsoft Corporation
//           2012-2013 Johns Hopkins University (author: Daniel Povey)

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

#define HAVE_ATLAS 1

#include "decoder.h"
#include "fst.h"
#include "vector.h"
#include <algorithm>
#include <stdio.h>
#include <math.h>

void pk_decoder_result_destroy(pk_decoder_result_t *best_path) {
  best_path->alignment_size = 0;
  pk_free(best_path->alloc, best_path->alignment);
  best_path->alignment = NULL;

  best_path->size = 0;
  pk_free(best_path->alloc, best_path->words);
  best_path->words = NULL;

  best_path->weight = 0;
  best_path->alloc = NULL;
}

pk_decoder_token_t *pk_decoder_token_new(
    pk_alloc_t *alloc,
    const pk_fst_arc_t *arc,
    double acoustic_cost,
    pk_decoder_token_t *previous) {
  pk_decoder_token_t *self = (pk_decoder_token_t *)malloc(
      sizeof(pk_decoder_token_t));
  self->previous = previous;
  self->ref_count = 1;
  self->input_label = arc->input_label;
  self->output_label = arc->output_label;
  self->next_state = arc->next_state;
  self->arc_weight = arc->weight + acoustic_cost;
  if (previous) {
    ++previous->ref_count;
    self->cost = previous->cost + arc->weight + acoustic_cost;
  } else {
    self->cost = arc->weight + acoustic_cost;
  }
  return self;
}

void pk_decoder_token_delete(pk_alloc_t *alloc, pk_decoder_token_t *self) {
  pk_decoder_token_t *token = self;
  while (--(token->ref_count) == 0) {
    pk_decoder_token_t *previous = token->previous;
    token->previous = NULL;
    token->ref_count = 0;
    self->input_label = 0;
    self->output_label = 0;
    self->arc_weight = 0;
    self->next_state = 0;
    self->cost = 0;
    pk_free(alloc, token);
    if (previous == NULL) break;
    token = previous;
  }
}

PKVECTOR_DEFINE(pk_decoder_token_t *, token_vector)
PKVECTOR_DEFINE(int, int_vector)
PKVECTOR_DEFINE(int32_t, int32_vector)

namespace kaldi {

PkSimpleDecoder::PkSimpleDecoder(
    const pk_fst_t *fst,
    BaseFloat beam): fst_(fst), beam_(beam) {
  pk_alloc_init(&alloc_);
  pk_hashlist_init(&cur_toks_, &alloc_);
  pk_hashlist_init(&prev_toks_, &alloc_);
  pk_hashlist_init(&tmp_toks_, &alloc_);
}

PkSimpleDecoder::~PkSimpleDecoder() {
  ClearToks(&cur_toks_);
  ClearToks(&prev_toks_);
  pk_hashlist_destroy(&cur_toks_);
  pk_hashlist_destroy(&prev_toks_);
  pk_hashlist_destroy(&tmp_toks_);
}


bool PkSimpleDecoder::Decode(DecodableInterface *decodable) {
  InitDecoding();
  while( !decodable->IsLastFrame(num_frames_decoded_ - 1)) {
    ClearToks(&prev_toks_);
    pk_hashlist_swap(&cur_toks_, &prev_toks_);
    ProcessEmitting(decodable);
    ProcessNonemitting();
    PruneToks(beam_, &cur_toks_);
  }
  return cur_toks_.size > 0;
}

void PkSimpleDecoder::InitDecoding() {
  // clean up from last time:
  ClearToks(&cur_toks_);
  ClearToks(&prev_toks_);
  // initialize decoding:
  int start_state = pk_fst_startstate(fst_);
  KALDI_ASSERT(start_state != fst::kNoStateId);
  pk_fst_arc_t dummy_arc;
  dummy_arc.input_label = 0;
  dummy_arc.output_label = 0;
  dummy_arc.next_state = start_state;
  dummy_arc.weight = 0;

  pk_decoder_token_t *start_token = pk_decoder_token_new(
      &alloc_,
      &dummy_arc,
      0,
      NULL);
  pk_hashlist_insert(&cur_toks_, start_state, start_token);
  num_frames_decoded_ = 0;
  ProcessNonemitting();
}

void PkSimpleDecoder::AdvanceDecoding(DecodableInterface *decodable,
                                      int max_num_frames) {
  KALDI_ASSERT(num_frames_decoded_ >= 0 &&
               "You must call InitDecoding() before AdvanceDecoding()");
  int num_frames_ready = decodable->NumFramesReady();
  // num_frames_ready must be >= num_frames_decoded, or else
  // the number of frames ready must have decreased (which doesn't
  // make sense) or the decodable object changed between calls
  // (which isn't allowed).
  KALDI_ASSERT(num_frames_ready >= num_frames_decoded_);
  int target_frames_decoded = num_frames_ready;
  if (max_num_frames >= 0)
    target_frames_decoded = std::min(target_frames_decoded,
                                     num_frames_decoded_ + max_num_frames);
  while (num_frames_decoded_ < target_frames_decoded) {
    // note: ProcessEmitting() increments num_frames_decoded_
    ClearToks(&prev_toks_);
    pk_hashlist_swap(&cur_toks_, &prev_toks_);
    ProcessEmitting(decodable);
    ProcessNonemitting();
    PruneToks(beam_, &cur_toks_);
  }   
}

bool PkSimpleDecoder::ReachedFinal() const {
  pk_hashlist_elem_t *elem = cur_toks_.head;
  while (elem) {
    int state = elem->key;
    pk_decoder_token_t *token = (pk_decoder_token_t *)elem->value;
    if (token->cost != std::numeric_limits<BaseFloat>::infinity() &&
        pk_fst_final(fst_, state) != INFINITY) {
      return true;
    }
    elem = elem->next;
  }
  return false;
}

BaseFloat PkSimpleDecoder::FinalRelativeCost() const {
  // as a special case, if there are no active tokens at all (e.g. some kind of
  // pruning failure), return infinity.
  double infinity = std::numeric_limits<double>::infinity();
  if (cur_toks_.size == 0) {
    return infinity;
  }
  double best_cost = infinity;
  double best_cost_with_final = infinity;
  pk_hashlist_elem_t *elem = cur_toks_.head;
  while (elem) {
    int state = elem->key;
    pk_decoder_token_t *token = (pk_decoder_token_t *)elem->value;

    // Note: Plus is taking the minimum cost, since we're in the tropical
    // semiring.
    best_cost = std::min(best_cost, (double)token->cost);
    best_cost_with_final = std::min(
        best_cost_with_final,
        double(token->cost + pk_fst_final(fst_, state)));
    elem = elem->next;
  }

  BaseFloat extra_cost = best_cost_with_final - best_cost;
  if (extra_cost != extra_cost) { // NaN.  This shouldn't happen; it indicates some
                                  // kind of error, most likely.
    KALDI_WARN << "Found NaN (likely search failure in decoding)";
    return infinity;
  }
  // Note: extra_cost will be infinity if no states were final.
  return extra_cost;
}

// Outputs an FST corresponding to the single best path
// through the lattice.
bool PkSimpleDecoder::GetBestPath(
    pk_decoder_result_t *best_path,
    bool use_final_probs)  {
  pk_decoder_token_t *best_tok = NULL;
  bool is_final = ReachedFinal();
  if (!is_final) {
    pk_hashlist_elem_t *elem = cur_toks_.head;
    while (elem) {
      pk_decoder_token_t *token = (pk_decoder_token_t *)elem->value;
      if (best_tok == NULL || best_tok->cost > token->cost) {
        best_tok = token;
      }
      elem = elem->next;
    }
  } else {
    double infinity = std::numeric_limits<double>::infinity();
    double best_cost = infinity;

    pk_hashlist_elem_t *elem = cur_toks_.head;
    while (elem) {
      int state = elem->key;
      pk_decoder_token_t *token = (pk_decoder_token_t *)elem->value;
      double this_cost = token->cost + pk_fst_final(fst_, state);
      if (this_cost != infinity && this_cost < best_cost) {
        best_cost = this_cost;
        best_tok = token;
      }
      elem = elem->next;
    }
  }
  if (best_tok == NULL) return false;  // No output.

  token_vector_t arcs_reverse;
  token_vector_init(&arcs_reverse, &alloc_);
  for (pk_decoder_token_t *tok = best_tok; tok != NULL; tok = tok->previous) {
    token_vector_push_back(&arcs_reverse, tok);
  }
  pk_decoder_token_t *back = token_vector_back(&arcs_reverse);
  assert(back->next_state == pk_fst_startstate(fst_));
  
  // that was a "fake" token... gives no info.
  token_vector_pop_back(&arcs_reverse);  

  int32_vector_t alignment;
  int32_vector_t words;
  int32_vector_init(&alignment, &alloc_);
  int32_vector_init(&words, &alloc_);
  float weight = 0;
  for (int idx = arcs_reverse.size - 1; idx >= 0; idx--) {
    if (arcs_reverse.data[idx]->input_label != 0)  {
      int32_vector_push_back(&alignment, arcs_reverse.data[idx]->input_label);
    }
    if (arcs_reverse.data[idx]->output_label != 0) {
      int32_vector_push_back(&words, arcs_reverse.data[idx]->output_label);
    }
    weight += arcs_reverse.data[idx]->arc_weight;
  }

  best_path->alignment = (int32_t *)pk_alloc(
      &alloc_,
      alignment.size * sizeof(int32_t));
  for (int idx = 0; idx < alignment.size; ++idx) {
    best_path->alignment[idx] = alignment.data[idx];
  }
  best_path->alignment_size = alignment.size;
  best_path->words = (int32_t *)pk_alloc(
      &alloc_,
      words.size * sizeof(int32_t));
  for (int idx = 0; idx < words.size; ++idx) {
    best_path->words[idx] = words.data[idx];
  }
  best_path->size = words.size;
  best_path->weight = weight;
  best_path->alloc = &alloc_;

  if (is_final && use_final_probs) {
    best_path->weight += pk_fst_final(fst_, best_tok->next_state);
  }

  token_vector_destroy(&arcs_reverse);
  int32_vector_destroy(&alignment);
  int32_vector_destroy(&words);
  return true;
}


void PkSimpleDecoder::ProcessEmitting(DecodableInterface *decodable) {
  int32 frame = num_frames_decoded_;
  // Processes emitting arcs for one frame.  Propagates from
  // prev_toks_ to cur_toks_.
  double cutoff = std::numeric_limits<BaseFloat>::infinity();
  pk_hashlist_elem_t *elem = prev_toks_.head;
  while (elem) {
    int state = elem->key;
    pk_decoder_token_t *tok = (pk_decoder_token_t *)elem->value;
    KALDI_ASSERT(state == tok->next_state);
    pk_fst_iter_t arc_iter;
    pk_fst_iterate_arc(fst_, state, &arc_iter);
    const pk_fst_arc_t *arc = NULL;
    while ((arc = pk_fst_iter_next(&arc_iter)) != NULL) {
      if (arc->input_label != 0) {  // propagate..
        BaseFloat acoustic_cost = -decodable->LogLikelihood(
            frame,
            arc->input_label);
        double total_cost = tok->cost + arc->weight + acoustic_cost;
        
        if (total_cost > cutoff) continue;
        if (total_cost + beam_ < cutoff) {
          cutoff = total_cost + beam_;
        }
        pk_decoder_token_t *new_tok = pk_decoder_token_new(
            &alloc_,
            arc,
            acoustic_cost,
            tok); 
        pk_hashlist_elem_t *new_elem = pk_hashlist_find(
            &cur_toks_,
            (pk_hashlist_key_t)arc->next_state);
        if (new_elem == NULL) {
          pk_hashlist_insert(
              &cur_toks_,
              (pk_hashlist_key_t)arc->next_state,
              (pk_hashlist_value_t)new_tok);
        } else {
          pk_decoder_token_t *token = (pk_decoder_token_t *)new_elem->value;
          if (token->cost > new_tok->cost) {
            pk_decoder_token_delete(&alloc_, token);
            new_elem->value = (pk_hashlist_value_t)new_tok;
          } else {
            pk_decoder_token_delete(&alloc_, new_tok);
          }
        }
      }
    }
    elem = elem->next;
  }
  num_frames_decoded_++;
}

void PkSimpleDecoder::ProcessNonemitting() {
  // Processes nonemitting arcs for one frame.  Propagates within
  // cur_toks_.
  int_vector_t queue;
  int_vector_init(&queue, &alloc_);
  double infinity = std::numeric_limits<double>::infinity();
  double best_cost = infinity;
  pk_hashlist_elem_t *elem = cur_toks_.head;
  while (elem) {
    int state = elem->key;
    pk_decoder_token_t *token = (pk_decoder_token_t *)elem->value;

    int_vector_push_back(&queue, state);
    best_cost = std::min(best_cost, (double)token->cost);
    elem = elem->next;
  }
  double cutoff = best_cost + beam_;
  
  while (!int_vector_empty(&queue)) {
    int state = int_vector_back(&queue);
    int_vector_pop_back(&queue);
    elem = pk_hashlist_find(&cur_toks_, (pk_hashlist_key_t)state);
    KALDI_ASSERT(elem != NULL);
    pk_decoder_token_t *tok = (pk_decoder_token_t *)elem->value;
    KALDI_ASSERT(tok != NULL && state == tok->next_state);
    pk_fst_iter_t arc_iter;
    pk_fst_iterate_arc(fst_, state, &arc_iter);
    const pk_fst_arc_t *arc = NULL;
    while ((arc = pk_fst_iter_next(&arc_iter)) != NULL) {
      if (arc->input_label == 0) {  // propagate nonemitting only...
        const BaseFloat acoustic_cost = 0.0;
        pk_decoder_token_t *new_tok = pk_decoder_token_new(
            &alloc_,
            arc,
            acoustic_cost,
            tok);
        if (new_tok->cost > cutoff) {
          pk_decoder_token_delete(&alloc_, new_tok);
        } else {
          elem = pk_hashlist_find(&cur_toks_, arc->next_state);
          if (elem == NULL) {
            pk_hashlist_insert(
                &cur_toks_,
                arc->next_state,
                (pk_hashlist_value_t)new_tok);
            int_vector_push_back(&queue, arc->next_state);
          } else {
            pk_decoder_token_t *token = (pk_decoder_token_t *)elem->value;
            if (token->cost > new_tok->cost) {
              pk_decoder_token_delete(&alloc_, token);
              elem->value = (pk_hashlist_value_t)new_tok;
              int_vector_push_back(&queue, arc->next_state);
            } else {
              pk_decoder_token_delete(&alloc_, new_tok);
            }
          }
        }
      }
    }
  }

  int_vector_destroy(&queue);
}

// static
void PkSimpleDecoder::ClearToks(pk_hashlist_t *tok) {
  pk_hashlist_elem_t *elem = tok->head;
  while (elem) {
    pk_decoder_token_t *token = (pk_decoder_token_t *)elem->value;
    pk_decoder_token_delete(&alloc_, token);
    elem = elem->next;
  }
  pk_hashlist_clear(tok);
}

// static
void PkSimpleDecoder::PruneToks(BaseFloat beam, pk_hashlist_t *toks) {
  if (toks->size == 0) {
    KALDI_VLOG(2) <<  "No tokens to prune.\n";
    return;
  }
  double best_cost = std::numeric_limits<double>::infinity();
  pk_hashlist_elem_t *elem = toks->head;
  while (elem) {
    pk_decoder_token_t *token = (pk_decoder_token_t *)elem->value;
    best_cost = std::min(best_cost, (double)token->cost);
    elem = elem->next;
  }

  int_vector_t retained;
  int_vector_init(&retained, &alloc_);
  double cutoff = best_cost + beam;
  elem = toks->head;
  while (elem) {
    int state = (int)elem->key;
    pk_decoder_token_t *token = (pk_decoder_token_t *)elem->value;

    if (token->cost < cutoff) {
      int_vector_push_back(&retained, state);
    } else {
      pk_decoder_token_delete(&alloc_, token);
    }
    elem = elem->next;
  }
  pk_hashlist_clear(&tmp_toks_);
  for (int i = 0; i < retained.size; i++) {
    int state = retained.data[i];
    elem = pk_hashlist_find(toks, (pk_hashlist_key_t)state);
    KALDI_ASSERT(elem);
    pk_hashlist_insert(&tmp_toks_, (pk_hashlist_key_t)state, elem->value);
  }
  KALDI_VLOG(2) <<  "Pruned to " << retained.size << " toks.\n";
  pk_hashlist_swap(&tmp_toks_, toks);

  int_vector_destroy(&retained);
}

} // end namespace kaldi.
