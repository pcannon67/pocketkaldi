// Created at 2017-03-06

#ifndef POCKETKALDI_ONLINE_CMVN_H_
#define POCKETKALDI_ONLINE_CMVN_H_

#include <stdint.h>
#include "matrix.h"
#include "util.h"

#define PK_ONLINECMVN_WINDOW 600
#define PK_ONLINECMVN_GLOBALFRAMES 200

namespace pocketkaldi {

// Computing CMVN of an utterance. It will stores the intermediate state.
// Currently only mean computation is supported
class CMVN {
 public:
  // Initialize online cmvn. global_stats is the accumulated CMVN stats from
  // training data. raw_feats is raw feature matrix. (global_stats and raw_feats
  // are borrowed)
  CMVN(const pk_vector_t *global_stats, const pk_matrix_t *raw_feats);
  ~CMVN();

  // Gets the fests of frame after appling online CMVN.
  void GetFrame(int frame, pk_vector_t *feats);

 private:
  const pk_matrix_t *raw_feats_;
  Vector<float> global_stats_;
  Vector<float> cached_stats_;
  int cached_frame_;

  // Apply CMVN to frame
  void Apply(const VectorBase<float> &stats, VectorBase<float> *feats);

  // Computes the raw CMVN stats for frame. Store it into `stats`
  void ComputeStats(int frame, Vector<float> *stats);

  // Smooth the CMVN stats "stats", by possibly adding some stats from
  // global_stats
  void SmoothStats(Vector<float> *stats);
};

}  // namepsace pocketkaldi

#endif  // POCKETKALDI_ONLINE_CMVN_H_

