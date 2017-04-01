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

#include <stdio.h>
#include <math.h>
#include "decoder.h"
#include "fst.h"
#include "list.h"
#include "util.h"

PKLIST_DEFINE(pk_decoder_token_t *, token_list)
PKLIST_DEFINE(int, int_list)
PKLIST_DEFINE(int32_t, int32_list)

void pk_decoder_result_init(pk_decoder_result_t *self) {
  self->alignment = NULL;
  self->alignment_size = 0;
  self->words = NULL;
  self->size = 0;
  self->weight = 0.0f;
}

void pk_decoder_result_destroy(pk_decoder_result_t *best_path) {
  best_path->alignment_size = 0;
  pk_free(best_path->alignment);
  best_path->alignment = NULL;

  best_path->size = 0;
  pk_free(best_path->words);
  best_path->words = NULL;

  best_path->weight = 0;
}

pk_decoder_token_t *pk_decoder_token_new(
    const pk_fst_arc_t *arc,
    double acoustic_cost,
    pk_decoder_token_t *previous) {
  pk_decoder_token_t *self = (pk_decoder_token_t *)pk_alloc(
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

void pk_decoder_token_delete(pk_decoder_token_t *self) {
  pk_decoder_token_t *token = self;
  while (--(token->ref_count) == 0) {
    pk_decoder_token_t *previous = token->previous;
    token->previous = NULL;
    token->ref_count = 0;
    token->input_label = 0;
    token->output_label = 0;
    token->arc_weight = 0;
    token->next_state = 0;
    token->cost = 0;
    pk_free(token);
    if (previous == NULL) break;
    token = previous;
  }
}

// Clear and free all tokens in hashlist
static void clear_toks(pk_hashlist_t *tok) {
  pk_hashlist_elem_t *elem = tok->head;
  while (elem) {
    pk_decoder_token_t *token = (pk_decoder_token_t *)elem->value;
    pk_decoder_token_delete(token);
    elem = elem->next;
  }
  pk_hashlist_clear(tok);
}

// Processes nonemitting arcs for one frame.  Propagates within cur_toks_.
static void process_nonemitting(pk_decoder_t *self) {
  int_list_t queue;
  int_list_init(&queue);
  double best_cost = INFINITY;
  pk_hashlist_elem_t *elem = self->cur_toks.head;
  while (elem) {
    int state = elem->key;
    pk_decoder_token_t *token = (pk_decoder_token_t *)elem->value;

    int_list_push_back(&queue, state);
    best_cost = PK_MIN(best_cost, (double)token->cost);
    elem = elem->next;
  }
  double cutoff = best_cost + self->beam;
  
  while (!int_list_empty(&queue)) {
    int state = int_list_back(&queue);
    int_list_pop_back(&queue);
    elem = pk_hashlist_find(&self->cur_toks, (pk_hashlist_key_t)state);
    assert(elem != NULL);
    pk_decoder_token_t *tok = (pk_decoder_token_t *)elem->value;
    assert(tok != NULL && state == tok->next_state);
    pk_fst_iter_t arc_iter;
    pk_fst_iterate_arc(self->fst, state, &arc_iter);
    const pk_fst_arc_t *arc = NULL;
    while ((arc = pk_fst_iter_next(&arc_iter)) != NULL) {
      if (arc->input_label == 0) {  // propagate nonemitting only...
        const float acoustic_cost = 0.0;
        pk_decoder_token_t *new_tok = pk_decoder_token_new(
            arc,
            acoustic_cost,
            tok);
        if (new_tok->cost > cutoff) {
          pk_decoder_token_delete(new_tok);
        } else {
          elem = pk_hashlist_find(&self->cur_toks, arc->next_state);
          if (elem == NULL) {
            pk_hashlist_insert(
                &self->cur_toks,
                arc->next_state,
                (pk_hashlist_value_t)new_tok);
            int_list_push_back(&queue, arc->next_state);
          } else {
            pk_decoder_token_t *token = (pk_decoder_token_t *)elem->value;
            if (token->cost > new_tok->cost) {
              pk_decoder_token_delete(token);
              elem->value = (pk_hashlist_value_t)new_tok;
              int_list_push_back(&queue, arc->next_state);
            } else {
              pk_decoder_token_delete(new_tok);
            }
          }
        }
      }
    }
  }

  int_list_destroy(&queue);
}

// Initializes the decoding
static void init_decoding(pk_decoder_t *self) {
  // Clean up from last time
  clear_toks(&self->cur_toks);
  clear_toks(&self->prev_toks);

  // Initialize decoding:
  int start_state = pk_fst_startstate(self->fst);
  assert(start_state >= 0);
  pk_fst_arc_t dummy_arc;
  dummy_arc.input_label = 0;
  dummy_arc.output_label = 0;
  dummy_arc.next_state = start_state;
  dummy_arc.weight = 0;

  pk_decoder_token_t *start_token = pk_decoder_token_new(&dummy_arc, 0, NULL);
  pk_hashlist_insert(&self->cur_toks, start_state, start_token);
  self->num_frames_decoded = 0;
  process_nonemitting(self);
}


// Processes emitting arcs for one frame.  Propagates from prev_toks_ to
// cur_toks_.
static void process_emitting(
    pk_decoder_t *self,
    pk_decodable_t *decodable) {
  int32_t frame = self->num_frames_decoded;

  double cutoff = INFINITY;
  pk_hashlist_elem_t *elem = self->prev_toks.head;
  while (elem) {
    int state = elem->key;
    pk_decoder_token_t *tok = (pk_decoder_token_t *)elem->value;
    assert(state == tok->next_state);
    pk_fst_iter_t arc_iter;
    pk_fst_iterate_arc(self->fst, state, &arc_iter);
    const pk_fst_arc_t *arc = NULL;
    while ((arc = pk_fst_iter_next(&arc_iter)) != NULL) {
      if (arc->input_label != 0) {  // propagate..
        float acoustic_cost = -pk_decodable_loglikelihood(
            decodable,
            frame,
            arc->input_label);
        double total_cost = tok->cost + arc->weight + acoustic_cost;
        
        if (total_cost > cutoff) continue;
        if (total_cost + self->beam < cutoff) {
          cutoff = total_cost + self->beam;
        }
        pk_decoder_token_t *new_tok = pk_decoder_token_new(
            arc,
            acoustic_cost,
            tok); 
        pk_hashlist_elem_t *new_elem = pk_hashlist_find(
            &self->cur_toks,
            (pk_hashlist_key_t)arc->next_state);
        if (new_elem == NULL) {
          pk_hashlist_insert(
              &self->cur_toks,
              (pk_hashlist_key_t)arc->next_state,
              (pk_hashlist_value_t)new_tok);
        } else {
          pk_decoder_token_t *token = (pk_decoder_token_t *)new_elem->value;
          if (token->cost > new_tok->cost) {
            pk_decoder_token_delete(token);
            new_elem->value = (pk_hashlist_value_t)new_tok;
          } else {
            pk_decoder_token_delete(new_tok);
          }
        }
      }
    }
    elem = elem->next;
  }
  self->num_frames_decoded++;
}

static void prune_toks(
    pk_decoder_t *self,
    float beam,
    pk_hashlist_t *toks) {
  if (toks->size == 0) {
    return;
  }
  double best_cost = INFINITY;
  pk_hashlist_elem_t *elem = toks->head;
  while (elem) {
    pk_decoder_token_t *token = (pk_decoder_token_t *)elem->value;
    best_cost = PK_MIN(best_cost, (double)token->cost);
    elem = elem->next;
  }

  int_list_t retained;
  int_list_init(&retained);
  double cutoff = best_cost + beam;
  elem = toks->head;
  while (elem) {
    int state = (int)elem->key;
    pk_decoder_token_t *token = (pk_decoder_token_t *)elem->value;

    if (token->cost < cutoff) {
      int_list_push_back(&retained, state);
    } else {
      pk_decoder_token_delete(token);
    }
    elem = elem->next;
  }
  pk_hashlist_clear(&self->tmp_toks);
  for (int i = 0; i < retained.size; i++) {
    int state = retained.data[i];
    elem = pk_hashlist_find(toks, (pk_hashlist_key_t)state);
    assert(elem);
    pk_hashlist_insert(&self->tmp_toks, (pk_hashlist_key_t)state, elem->value);
  }
  pk_hashlist_swap(&self->tmp_toks, toks);

  int_list_destroy(&retained);
}

void pk_decoder_init(pk_decoder_t *self, const pk_fst_t *fst, float beam) {
  self->fst = fst;
  self->beam = beam;

  pk_hashlist_init(&self->cur_toks);
  pk_hashlist_init(&self->prev_toks);
  pk_hashlist_init(&self->tmp_toks);
}

void pk_decoder_destroy(pk_decoder_t *self) {
  clear_toks(&self->cur_toks);
  clear_toks(&self->prev_toks);
  pk_hashlist_destroy(&self->cur_toks);
  pk_hashlist_destroy(&self->prev_toks);
  pk_hashlist_destroy(&self->tmp_toks);
}

bool pk_decoder_decode(pk_decoder_t *self, pk_decodable_t *decodable) {
  init_decoding(self);
  while(!pk_decodable_islastframe(decodable, self->num_frames_decoded - 1)) {
    clear_toks(&self->prev_toks);
    pk_hashlist_swap(&self->cur_toks, &self->prev_toks);
    process_emitting(self, decodable);
    process_nonemitting(self);
    prune_toks(self, self->beam, &self->cur_toks);
  }
  return self->cur_toks.size > 0;
}

bool pk_decoder_reachedfinal(pk_decoder_t *self) {
  pk_hashlist_elem_t *elem = self->cur_toks.head;
  while (elem) {
    int state = elem->key;
    pk_decoder_token_t *token = (pk_decoder_token_t *)elem->value;
    if (token->cost != INFINITY &&
        pk_fst_final(self->fst, state) != INFINITY) {
      return true;
    }
    elem = elem->next;
  }
  return false;
}

// Outputs best_path corresponding to the single best path through the lattice.
bool pk_decoder_bestpath(
    pk_decoder_t *self,
    pk_decoder_result_t *best_path,
    bool use_final_probs)  {
  pk_decoder_token_t *best_tok = NULL;
  bool is_final = pk_decoder_reachedfinal(self);
  if (!is_final) {
    pk_hashlist_elem_t *elem = self->cur_toks.head;
    while (elem) {
      pk_decoder_token_t *token = (pk_decoder_token_t *)elem->value;
      if (best_tok == NULL || best_tok->cost > token->cost) {
        best_tok = token;
      }
      elem = elem->next;
    }
  } else {
    double infinity = INFINITY;
    double best_cost = infinity;

    pk_hashlist_elem_t *elem = self->cur_toks.head;
    while (elem) {
      int state = elem->key;
      pk_decoder_token_t *token = (pk_decoder_token_t *)elem->value;
      double this_cost = token->cost + pk_fst_final(self->fst, state);
      if (this_cost != infinity && this_cost < best_cost) {
        best_cost = this_cost;
        best_tok = token;
      }
      elem = elem->next;
    }
  }
  if (best_tok == NULL) return false;  // No output.

  token_list_t arcs_reverse;
  token_list_init(&arcs_reverse);
  for (pk_decoder_token_t *tok = best_tok; tok != NULL; tok = tok->previous) {
    token_list_push_back(&arcs_reverse, tok);
  }
  pk_decoder_token_t *back = token_list_back(&arcs_reverse);
  assert(back->next_state == pk_fst_startstate(self->fst));
  
  // that was a "fake" token... gives no info.
  token_list_pop_back(&arcs_reverse);  

  int32_list_t alignment;
  int32_list_t words;
  int32_list_init(&alignment);
  int32_list_init(&words);
  float weight = 0;
  for (int idx = arcs_reverse.size - 1; idx >= 0; idx--) {
    if (arcs_reverse.data[idx]->input_label != 0)  {
      int32_list_push_back(&alignment, arcs_reverse.data[idx]->input_label);
    }
    if (arcs_reverse.data[idx]->output_label != 0) {
      int32_list_push_back(&words, arcs_reverse.data[idx]->output_label);
    }
    weight += arcs_reverse.data[idx]->arc_weight;
  }

  best_path->alignment = (int32_t *)pk_alloc(
      alignment.size * sizeof(int32_t));
  for (int idx = 0; idx < alignment.size; ++idx) {
    best_path->alignment[idx] = alignment.data[idx];
  }
  best_path->alignment_size = alignment.size;
  best_path->words = (int32_t *)pk_alloc(
      words.size * sizeof(int32_t));
  for (int idx = 0; idx < words.size; ++idx) {
    best_path->words[idx] = words.data[idx];
  }
  best_path->size = words.size;
  best_path->weight = weight;

  if (is_final && use_final_probs) {
    best_path->weight += pk_fst_final(self->fst, best_tok->next_state);
  }

  token_list_destroy(&arcs_reverse);
  int32_list_destroy(&alignment);
  int32_list_destroy(&words);
  return true;
}

