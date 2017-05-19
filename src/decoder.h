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
#include <vector>
#include <unordered_map>
#include <unordered_map>
#include "transition.h"
#include "decodable.h"
#include "hashtable.h"
#include "vector.h"
#include "pool.h"
#include "fst.h"
#include "am.h"

namespace pocketkaldi {

// Decoder is the core class in pocketkaldi. It decodes the Decodable object
// using viterbi algorithm and stores the best result into Decoder::Hypothesis
// class
class Decoder {
 public:
  static constexpr int kBeamSize = 30000;
  static constexpr float kBeamDelta = 0.5;
  static constexpr int kOLabelBeginIdx = -1;
  static constexpr int kNotExist = -1;
  static constexpr int kCutoffSamples = 200;
  static constexpr int kCutoffRandSeed = 0x322;

  // Stores the decoding result
  class Hypothesis;

  // Initialize the decoder with the FST graph fst. It just borrows the pointer
  // of fst and not own it.
  Decoder(const Fst *fst);
  ~Decoder();

  // Decodes the Decodable object, the best path could be obtain by
  // Decoder::BestPath()
  bool Decode(pk_decodable_t *decodable);

  // Return true if reached the final state
  bool ReachedFinal();

  // Get best hypothesis from lattice.
  Hypothesis BestPath();

 private:
  // Token represents a state in the viterbi lattice. olabel_idx is the index
  // of its corresponded outpu label link-list in the list impl->olabels
  class Token;

  // OLabel is the struct records the list of output labels of a tok in beam. 
  // previous pointes to the previous OLabel like a link list.
  class OLabel;

  // Get the weight cutoff from prev_toks_. We won't go throuth all the toks in
  // beam here to calculate cutoff, Since it takes a long time. Instead, we
  // randomly sample N costs from the beam and GUESS the cutoff value.
  // \param adaptive_beam the gap of cost between best token and the token of
  // beam_size
  // \param best_tokidx index of best token in prev_toks_
  // \return the cutoff of cost to beam size in prev_toks_
  double GetCutoff(float *adaptive_beam, Token **best_tok);

  // Initialize decoding and put the root state into beam
  void InitDecoding();

  // Insert tok into self->toks_. The tok is created according to arc. And it
  // will either insert a new token or update existing token in the beam.
  // \param olabel_idx the olabel-index of previous token in viterbi
  // \param cost the cost of new token
  // \return true if successfully inserted. Otherwise, when the cost of
  // existing tok is less than new one, return false
  bool InsertTok(const Fst::Arc *arc, OLabel *olabel, float cost);

  // Processes nonemitting arcs for one frame. Propagates within cur_toks_.
  void ProcessNonemitting(double cutoff);

  // Process the emitting (non-epsilon) arcs of each states in the beam
  // \return cutoff of next weight
  float ProcessEmitting(pk_decodable_t *decodable);

  // Only used in get_cutoff()
  std::vector<float> costs_;

  // FST graph used for decoding
  const Fst *fst_;

  // Frames decoded
  int num_frames_decoded_;

  // Tokens used in decoding. toks_ is the current beam, all generated tokens
  // will be placed here. And after frame advanced, the toks_ will be swapped
  // with prev_toks_
  Pool<Token> toks_pool_;
  std::vector<Token *> toks_;
  std::vector<Token *> prev_toks_;

  // Stores the map between state-id and the index of corresponded token in
  // toks_
  HashTable<int32_t, int32_t> state_idx_;

  // Storea all output-label nodes
  GCPool<OLabel> olabels_pool_;

  // Beam threshold
  float beam_;
};

// Stores the decoding result
class Decoder::Hypothesis {
 public:
  Hypothesis(const std::vector<int> &words, float weight);

  // Word-ids in the decoding result
  const std::vector<int> &words() const { return words_; }

  // Weight for this utterance
  float weight() const { return weight_; }

 private:
  std::vector<int> words_;  
  float weight_;
};

class Decoder::Token {
 public:
  Token(int state, float cost, OLabel *olabel);

  // The state in FST
  int state() const { return state_; }

  // Current cost
  float cost() const { return cost_; }

  // Head of output label chain
  OLabel *olabel() const { return olabel_; }

 private:
  int state_;
  float cost_;
  OLabel *olabel_;
};

class Decoder::OLabel : public Collectable {
 public:
  OLabel(OLabel *previous, int olabel);

  // Index of previous olabel
  OLabel *previous() const { return previous_; }

  // Output label of current node
  int olabel() const { return olabel_; }

 private:
  OLabel *previous_;
  int olabel_;
};

}  // namespace pocketkaldi

#endif  // POCKETKALDI_DECODER_H_
