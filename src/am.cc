// Created at 2017-03-22

#include "am.h"

#include "status.h"
#include "matrix.h"
#include "math.h"

namespace pocketkaldi {

AcousticModel::AcousticModel() :
    left_context_(0),
    right_context_(0),
    num_pdfs_(0) {
}

AcousticModel::~AcousticModel() {
  left_context_ = 0;
  right_context_ = 0;
  num_pdfs_ = 0;
}

Status AcousticModel::Read(const Configuration &conf) {
  Status status;

  // Read nnet
  std::string nnet_filename;
  status = conf.GetPath("nnet", &nnet_filename);
  if (!status.ok()) return status;
  pk_status_t c_status;
  pk_status_init(&c_status);

  util::ReadableFile fd;
  PK_CHECK_STATUS(fd.Open(nnet_filename));
  PK_CHECK_STATUS(nnet_.Read(&fd));
  fd.Close();

  // Read prior
  std::string prior_filename;
  PK_CHECK_STATUS(conf.GetPath("prior", &prior_filename));
  PK_CHECK_STATUS(fd.Open(prior_filename));
  PK_CHECK_STATUS(log_prior_.Read(&fd));
  log_prior_.ApplyLog();
  fd.Close();

  // Read left and right context
  status = conf.GetInteger("left_context", &left_context_);
  status = conf.GetInteger("right_context", &right_context_);
  if (!status.ok()) return status;

  // Read tid2pdf_
  std::string tid2pdf_filename;
  status = conf.GetInteger("num_pdfs", &num_pdfs_);
  if (!status.ok()) return status;
  status = conf.GetPath("tid2pdf", &tid2pdf_filename);
  if (!status.ok()) return status;
  status = fd.Open(tid2pdf_filename);
  if (!status.ok()) return status;
  status = tid2pdf_.Read(&fd);
  if (!status.ok()) return status;

  return Status::OK();
}

void AcousticModel::SpliceFeats(
    const pk_matrix_t *feats,
    pk_matrix_t *spliced_feats) {
  for (int frame_idx = 0; frame_idx < feats->ncol; ++frame_idx) {
    pk_vector_t spliced_x = pk_matrix_getcol(spliced_feats, frame_idx);
    int offset = 0;
    for (int f = -left_context_; f <= right_context_; ++f) {
      int from_frame_idx = frame_idx + f;

      // Padding boundary feats
      if (from_frame_idx < 0) from_frame_idx = 0;
      if (from_frame_idx >= feats->ncol) from_frame_idx = feats->ncol - 1;

      pk_vector_t frame = pk_matrix_getcol(feats, from_frame_idx);
      for (int d = 0; d < frame.dim; ++d) {
        assert(offset + d < spliced_x.dim);
        spliced_x.data[offset + d] = frame.data[d];
      }

      // Increase offset
      offset += frame.dim;
    }
  }
}

void AcousticModel::Compute(
    const pk_matrix_t *frames,
    pk_matrix_t *loglikelihood) {
  // Prepare spliced feature matrix
  int feats_dim = frames->nrow;
  int spliced_dim = (left_context_ + right_context_ + 1) * feats_dim;
  pk_matrix_t nn_input;
  pk_matrix_init(&nn_input, spliced_dim, frames->ncol);
  
  // Splice the feats
  SpliceFeats(frames, &nn_input);

  // Propogate through the neural network
  nnet_.Propagate(&nn_input, loglikelihood);

  // Compute log-likelihood
  for (int col_idx = 0; col_idx < loglikelihood->ncol; ++col_idx) {
    pk_vector_t c_col = pk_matrix_getcol(loglikelihood, col_idx);
    SubVector<float> col(c_col.data, c_col.dim);
    col.ApplyFloor(1.0e-20);
    col.ApplyLog();
    col.AddVec(-1.0f, log_prior_);
  }

  pk_matrix_destroy(&nn_input);
}

}  // namespace pocketkaldi
