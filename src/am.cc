// Created at 2017-03-22

#include "am.h"

#include "status.h"
#include "matrix.h"
#include "math.h"

using pocketkaldi::Status;

void pk_am_init(pk_am_t *self) {
  pk_vector_init(&self->log_prior, 0, NAN);
  self->left_context = 0;
  self->right_context = 0;
}

void pk_am_read(pk_am_t *self, pk_readable_t *fd, pk_status_t *status) {
  Status vn_status;
  pk_bytebuffer_t bytebuffer;
  pk_bytebuffer_init(&bytebuffer);

  int section_size = pk_readable_readsectionhead(
      fd,
      PK_AM_SECTION,
      status);
  if (!status->ok) goto pk_am_read_failed;
  pk_bytebuffer_reset(&bytebuffer, section_size);

  if (section_size != 8) {
    PK_STATUS_CORRUPTED(
        status,
        "pk_am_read: section_size == %d expected, but %d found (%s)",
        8,
        section_size,
        fd->filename);
    goto pk_am_read_failed;
  }

  // Read context
  pk_readable_readbuffer(fd, &bytebuffer, status);
  if (!status->ok) goto pk_am_read_failed;
  self->left_context = pk_bytebuffer_readint32(&bytebuffer);
  self->right_context = pk_bytebuffer_readint32(&bytebuffer);

  // Read nnet
  vn_status = self->nnet.Read(fd);
  if (!vn_status.ok()) goto pk_am_read_failed;

  // Read prior and apply log
  pk_vector_read(&self->log_prior, fd, status);
  if (!status->ok) goto pk_am_read_failed;
  pk_vector_log(&self->log_prior);

  if (false) {
pk_am_read_failed:
    pk_am_destroy(self);
  }

  pk_bytebuffer_destroy(&bytebuffer);
}

void pk_am_destroy(pk_am_t *self) {
  pk_vector_destroy(&self->log_prior);
  self->left_context = 0;
  self->right_context = 0;
}

// Splice the feature matrix 'feats' with the contexts specified by
// 'left_context' and 'right_context'. Store the result into 'spliced_feats'
static void splice_feats(
    const pk_matrix_t *feats,
    int left_context,
    int right_context,
    pk_matrix_t *spliced_feats) {
  for (int frame_idx = 0; frame_idx < feats->ncol; ++frame_idx) {
    pk_vector_t spliced_x = pk_matrix_getcol(spliced_feats, frame_idx);
    int offset = 0;
    for (int f = -left_context; f <= right_context; ++f) {
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


void pk_am_compute(
    const pk_am_t *self,
    const pk_matrix_t *frames,
    pk_matrix_t *loglikelihood) {
  // Prepare spliced feature matrix
  int feats_dim = frames->nrow;
  int spliced_dim = (self->left_context + self->right_context + 1) * feats_dim;
  pk_matrix_t nn_input;
  pk_matrix_init(&nn_input, spliced_dim, frames->ncol);
  
  // Splice the feats
  splice_feats(frames, self->left_context, self->right_context, &nn_input);

  // Propogate through the neural network
  self->nnet.Propagate(&nn_input, loglikelihood);

  // Compute log-likelihood
  for (int col_idx = 0; col_idx < loglikelihood->ncol; ++col_idx) {
    pk_vector_t col = pk_matrix_getcol(loglikelihood, col_idx);
    pk_vector_floor(&col, 1.0e-20);
    pk_vector_log(&col);
    pk_vector_addscale(&col, -1.0f, &self->log_prior);
  }

  pk_matrix_destroy(&nn_input);
}
