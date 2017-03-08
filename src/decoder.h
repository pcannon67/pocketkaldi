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

#ifndef KALDI_DECODER_SIMPLE_DECODER_H_
#define KALDI_DECODER_SIMPLE_DECODER_H_

#include <stdbool.h>
#include <stdint.h>
#include "hashlist.h"
#include "fst.h"
#include "util/stl-utils.h"
#include "lat/kaldi-lattice.h"
#include "itf/decodable-itf.h"

typedef struct {
  int32_t *alignment;
  int32_t alignment_size;
  int32_t *words;
  int32_t size;
  float weight;
} pk_decoder_result_t;

// Destroy the pk_decoder_result_t.
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

// Allocate and initialize a token for decoder, return the pointer to it
pk_decoder_token_t *pk_decoder_token_new(
    const pk_fst_arc_t *arc,
    double acoustic_cost,
    pk_decoder_token_t *previous);

// Delete the token when there is no reference to it
void pk_decoder_token_delete(pk_decoder_token_t *self);


namespace kaldi {

/** Simplest possible decoder, included largely for didactic purposes and as a
    means to debug more highly optimized decoders.  See \ref decoders_simple
    for more information.
 */
class PkSimpleDecoder {
 public:
  PkSimpleDecoder(const pk_fst_t *fst, BaseFloat beam);
  ~PkSimpleDecoder();

  /// Decode this utterance.
  /// Returns true if any tokens reached the end of the file (regardless of
  /// whether they are in a final state); query ReachedFinal() after Decode()
  /// to see whether we reached a final state.
  bool Decode(DecodableInterface *decodable);

  bool ReachedFinal() const;

  // GetBestPath gets the decoding traceback. If "use_final_probs" is true
  // AND we reached a final state, it limits itself to final states;
  // otherwise it gets the most likely token not taking into account final-probs.
  // fst_out will be empty (Start() == kNoStateId) if nothing was available due to
  // search error.
  // If Decode() returned true, it is safe to assume GetBestPath will return true.
  // It returns true if the output lattice was nonempty (i.e. had states in it);
  // using the return value is deprecated.
  bool GetBestPath(pk_decoder_result_t *best_path, bool use_final_probs = true) ;
  
  /// *** The next functions are from the "new interface". ***
  
  /// FinalRelativeCost() serves the same function as ReachedFinal(), but gives
  /// more information.  It returns the difference between the best (final-cost plus
  /// cost) of any token on the final frame, and the best cost of any token
  /// on the final frame.  If it is infinity it means no final-states were present
  /// on the final frame.  It will usually be nonnegative.
  BaseFloat FinalRelativeCost() const;

  /// InitDecoding initializes the decoding, and should only be used if you
  /// intend to call AdvanceDecoding().  If you call Decode(), you don't need
  /// to call this.  You can call InitDecoding if you have already decoded an
  /// utterance and want to start with a new utterance. 
  void InitDecoding();  

  /// This will decode until there are no more frames ready in the decodable
  /// object, but if max_num_frames is >= 0 it will decode no more than
  /// that many frames.  If it returns false, then no tokens are alive,
  /// which is a kind of error state.
  void AdvanceDecoding(DecodableInterface *decodable,
                         int32 max_num_frames = -1);
  
  /// Returns the number of frames already decoded.  
  int32 NumFramesDecoded() const { return num_frames_decoded_; }

 private:
  // ProcessEmitting decodes the frame num_frames_decoded_ of the
  // decodable object, then increments num_frames_decoded_.
  void ProcessEmitting(DecodableInterface *decodable);

  void ProcessNonemitting();
  
  pk_hashlist_t cur_toks_;
  pk_hashlist_t prev_toks_;
  pk_hashlist_t tmp_toks_;

  const pk_fst_t *fst_;
  BaseFloat beam_;
  // Keep track of the number of frames decoded in the current file.
  int32 num_frames_decoded_;
  
  void ClearToks(pk_hashlist_t *tok);

  void PruneToks(BaseFloat beam, pk_hashlist_t *toks);
  
  KALDI_DISALLOW_COPY_AND_ASSIGN(PkSimpleDecoder);
};


} // end namespace kaldi.


#endif
