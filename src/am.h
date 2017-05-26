// Created at 2017-03-22

#ifndef POCKETKALDI_AM_H_
#define POCKETKALDI_AM_H_

#include <vector>
#include "nnet.h"
#include "util.h"
#include "configuration.h"

#define PK_AM_SECTION "AM~0"

using pocketkaldi::Nnet;

namespace pocketkaldi {

// Acoustic model in ASR, inculde
//   - Neural network model
//   - Prior for each CD-state
//   - Map from transition-id to pdf-id
class AcousticModel {
 public:
  AcousticModel();
  ~AcousticModel();

  // Read AcousticModel from configuration file
  Status Read(const Configuration &conf);

  // Convert transition-id to pdf-id
  int TransitionIdToPdfId(int transition_id) const {
    return tid2pdf_(transition_id);
  }

  // Compute the log-likelihood of the feature matrix
  void Compute(const pk_matrix_t *frames, pk_matrix_t *loglikelihood);
 
  // Number of PDFs in this AM
  int num_pdfs() const { return num_pdfs_; }

 private:
  Nnet nnet_;
  Vector<float> log_prior_;
  int left_context_;
  int right_context_;
  int num_pdfs_;
  Vector<int32_t> tid2pdf_;

  // Splice the feature matrix 'feats' with the contexts specified by
  // 'left_context' and 'right_context'. Store the result into 'spliced_feats'
  void SpliceFeats(const pk_matrix_t *feats, pk_matrix_t *spliced_feats);
};

}  // namespace pocketkaldi


#endif
