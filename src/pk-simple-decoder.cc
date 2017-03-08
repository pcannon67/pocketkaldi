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

#include "pk-simple-decoder.h"
#include "fstext/remove-eps-local.h"
#include <algorithm>
#include <stdio.h>

namespace kaldi {

PkSimpleDecoder::PkSimpleDecoder(
    const fst::Fst<fst::StdArc> &fst,
    BaseFloat beam): fst_(fst), beam_(beam) {
  pk_alloc_init(&alloc_);
  pk_hashlist_init(&cur_toks_, &alloc_);
  pk_hashlist_init(&prev_toks_, &alloc_);
  pk_hashlist_init(&tmp_toks_, &alloc_);
}

PkSimpleDecoder::~PkSimpleDecoder() {
  puts("~PkSimpleDecoder");
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
  StateId start_state = fst_.Start();
  KALDI_ASSERT(start_state != fst::kNoStateId);
  StdArc dummy_arc(0, 0, StdWeight::One(), start_state);
  Token *start_token = new Token(dummy_arc, 0.0, NULL);
  pk_hashlist_insert(&cur_toks_, start_state, start_token);
  num_frames_decoded_ = 0;
  ProcessNonemitting();
}

void PkSimpleDecoder::AdvanceDecoding(DecodableInterface *decodable,
                                      int32 max_num_frames) {
  KALDI_ASSERT(num_frames_decoded_ >= 0 &&
               "You must call InitDecoding() before AdvanceDecoding()");
  int32 num_frames_ready = decodable->NumFramesReady();
  // num_frames_ready must be >= num_frames_decoded, or else
  // the number of frames ready must have decreased (which doesn't
  // make sense) or the decodable object changed between calls
  // (which isn't allowed).
  KALDI_ASSERT(num_frames_ready >= num_frames_decoded_);
  int32 target_frames_decoded = num_frames_ready;
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
    StateId state = (StateId)elem->key;
    Token *token = (Token *)elem->value;
    if (token->cost_ != std::numeric_limits<BaseFloat>::infinity() &&
        fst_.Final(state) != StdWeight::Zero()) {
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
    StateId state = (StateId)elem->key;
    Token *token = (Token *)elem->value;

    // Note: Plus is taking the minimum cost, since we're in the tropical
    // semiring.
    best_cost = std::min(best_cost, token->cost_);
    best_cost_with_final = std::min(
        best_cost_with_final,
        token->cost_ + fst_.Final(state).Value());
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
bool PkSimpleDecoder::GetBestPath(Lattice *fst_out, bool use_final_probs) const {
  fst_out->DeleteStates();
  Token *best_tok = NULL;
  bool is_final = ReachedFinal();
  if (!is_final) {
    pk_hashlist_elem_t *elem = cur_toks_.head;
    while (elem) {
      Token *token = (Token *)elem->value;
      if (best_tok == NULL || *best_tok < *token) {
        best_tok = token;
      }
      elem = elem->next;
    }
  } else {
    double infinity = std::numeric_limits<double>::infinity();
    double best_cost = infinity;

    pk_hashlist_elem_t *elem = cur_toks_.head;
    while (elem) {
      StateId state = (StateId)elem->key;
      Token *token = (Token *)elem->value;
      double this_cost = token->cost_ + fst_.Final(state).Value();
      if (this_cost != infinity && this_cost < best_cost) {
        best_cost = this_cost;
        best_tok = token;
      }
      elem = elem->next;
    }
  }
  if (best_tok == NULL) return false;  // No output.

  std::vector<LatticeArc> arcs_reverse;  // arcs in reverse order.
  for (Token *tok = best_tok; tok != NULL; tok = tok->prev_)
    arcs_reverse.push_back(tok->arc_);
  KALDI_ASSERT(arcs_reverse.back().nextstate == fst_.Start());
  arcs_reverse.pop_back();  // that was a "fake" token... gives no info.

  StateId cur_state = fst_out->AddState();
  fst_out->SetStart(cur_state);
  for (ssize_t i = static_cast<ssize_t>(arcs_reverse.size())-1; i >= 0; i--) {
    LatticeArc arc = arcs_reverse[i];
    arc.nextstate = fst_out->AddState();
    fst_out->AddArc(cur_state, arc);
    cur_state = arc.nextstate;
  }
  if (is_final && use_final_probs)
    fst_out->SetFinal(cur_state,
                      LatticeWeight(fst_.Final(best_tok->arc_.nextstate).Value(),
                                    0.0));
  else
    fst_out->SetFinal(cur_state, LatticeWeight::One());
  fst::RemoveEpsLocal(fst_out);
  return true;
}


void PkSimpleDecoder::ProcessEmitting(DecodableInterface *decodable) {
  int32 frame = num_frames_decoded_;
  // Processes emitting arcs for one frame.  Propagates from
  // prev_toks_ to cur_toks_.
  double cutoff = std::numeric_limits<BaseFloat>::infinity();
  pk_hashlist_elem_t *elem = prev_toks_.head;
  while (elem) {
    StateId state = (StateId)elem->key;
    Token *tok = (Token *)elem->value;
    if (state != tok->arc_.nextstate) printf("%d, %d\n", state, tok->arc_.nextstate);
    KALDI_ASSERT(state == tok->arc_.nextstate);
    for (fst::ArcIterator<fst::Fst<StdArc> > aiter(fst_, state);
         !aiter.Done();
         aiter.Next()) {
      const StdArc &arc = aiter.Value();
      if (arc.ilabel != 0) {  // propagate..
        BaseFloat acoustic_cost = -decodable->LogLikelihood(frame, arc.ilabel);
        double total_cost = tok->cost_ + arc.weight.Value() + acoustic_cost;
        
        if (total_cost > cutoff) continue;
        if (total_cost + beam_  < cutoff)
          cutoff = total_cost + beam_;
        Token *new_tok = new Token(arc, acoustic_cost, tok);
        pk_hashlist_elem_t *new_elem = pk_hashlist_find(
            &cur_toks_,
            (pk_hashlist_key_t)arc.nextstate);
        if (new_elem == NULL) {
          pk_hashlist_insert(
              &cur_toks_,
              (pk_hashlist_key_t)arc.nextstate,
              (pk_hashlist_value_t)new_tok);
        } else {
          Token *token = (Token *)new_elem->value;
          if (*token < *new_tok) {
            Token::TokenDelete(token);
            new_elem->value = (pk_hashlist_value_t)new_tok;
          } else {
            Token::TokenDelete(new_tok);
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
  std::vector<StateId> queue_;
  double infinity = std::numeric_limits<double>::infinity();
  double best_cost = infinity;
  pk_hashlist_elem_t *elem = cur_toks_.head;
  while (elem) {
    StateId state = (StateId)elem->key;
    Token *token = (Token *)elem->value;

    queue_.push_back(state);
    best_cost = std::min(best_cost, token->cost_);
    elem = elem->next;
  }
  double cutoff = best_cost + beam_;
  
  while (!queue_.empty()) {
    StateId state = queue_.back();
    queue_.pop_back();
    elem = pk_hashlist_find(&cur_toks_, (pk_hashlist_key_t)state);
    KALDI_ASSERT(elem != NULL);
    Token *tok = (Token *)elem->value;
    KALDI_ASSERT(tok != NULL && state == tok->arc_.nextstate);
    for (fst::ArcIterator<fst::Fst<StdArc> > aiter(fst_, state);
         !aiter.Done();
         aiter.Next()) {
      const StdArc &arc = aiter.Value();
      if (arc.ilabel == 0) {  // propagate nonemitting only...
        const BaseFloat acoustic_cost = 0.0;
        Token *new_tok = new Token(arc, acoustic_cost, tok);
        if (new_tok->cost_ > cutoff) {
          Token::TokenDelete(new_tok);
        } else {
          elem = pk_hashlist_find(&cur_toks_, (pk_hashlist_key_t)arc.nextstate);
          if (elem == NULL) {
            pk_hashlist_insert(
                &cur_toks_,
                (pk_hashlist_key_t)arc.nextstate,
                (pk_hashlist_value_t)new_tok);
            queue_.push_back(arc.nextstate);
          } else {
            Token *token = (Token *)elem->value;
            if (*token < *new_tok) {
              Token::TokenDelete(token);
              elem->value = (pk_hashlist_value_t)new_tok;
              queue_.push_back(arc.nextstate);
            } else {
              Token::TokenDelete(new_tok);
            }
          }
        }
      }
    }
  }
}

// static
void PkSimpleDecoder::ClearToks(pk_hashlist_t *tok) {
  pk_hashlist_elem_t *elem = tok->head;
  while (elem) {
    Token *token = (Token *)elem->value;
    Token::TokenDelete(token);
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
    Token *token = (Token *)elem->value;
    best_cost = std::min(best_cost, token->cost_);
    elem = elem->next;
  }

  std::vector<StateId> retained;
  double cutoff = best_cost + beam;
  elem = toks->head;
  while (elem) {
    StateId state = (StateId)elem->key;
    Token *token = (Token *)elem->value;

    if (token->cost_ < cutoff) {
      retained.push_back(state);
    } else {
      Token::TokenDelete(token);
    }
    elem = elem->next;
  }
  pk_hashlist_clear(&tmp_toks_);
  for (size_t i = 0; i < retained.size(); i++) {
    StateId state = retained[i];
    elem = pk_hashlist_find(toks, (pk_hashlist_key_t)state);
    KALDI_ASSERT(elem);
    pk_hashlist_insert(&tmp_toks_, (pk_hashlist_key_t)state, elem->value);
  }
  KALDI_VLOG(2) <<  "Pruned to " << (retained.size()) << " toks.\n";
  pk_hashlist_swap(&tmp_toks_, toks);
}

} // end namespace kaldi.
