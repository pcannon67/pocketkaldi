// Created at 2017-03-22

#include "am.h"
#include "matrix.h"
#include "math.h"

void pk_am_init(pk_am_t *self) {
  pk_nnet_init(&self->nnet);
  pk_vector_init(&self->log_prior, 0, NAN);
  self->left_context = 0;
  self->right_context = 0;
}

void pk_am_read(pk_am_t *self, pk_readable_t *fd, pk_status_t *status) {
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
  pk_nnet_read(&self->nnet, fd, status);
  if (!status->ok) goto pk_am_read_failed;

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
  pk_nnet_destroy(&self->nnet);
  pk_vector_destroy(&self->log_prior);
  self->left_context = 0;
  self->right_context = 0;
}

void pk_am_compute(
    const pk_am_t *self,
    const pk_matrix_t *frames,
    int frame_idx,
    pk_vector_t *loglikelihood) {
  assert(frame_idx < frames->ncol && "frame_idx out-of-boundary");

  // Prepare spliced input vector
  int dim = frames->nrow;
  int total_dim = (self->left_context + self->right_context + 1) * dim;
  pk_vector_t x;
  pk_vector_init(&x, total_dim, NAN);
  
  int offset = 0;
  for (int f = -self->left_context; f <= self->right_context; ++f) {
    int from_frame_idx = frame_idx + f;

    // Padding boundary frames
    if (from_frame_idx < 0) from_frame_idx = 0;
    if (from_frame_idx >= frames->ncol) from_frame_idx = frames->ncol - 1;

    pk_vector_t frame = pk_matrix_getcol(frames, from_frame_idx);
    for (int d = 0; d < frame.dim; ++d) {
      x.data[offset + d] = frame.data[d];
    }

    // Increase offset
    offset += frame.dim;
  }

  // Propogate through the neural network
  pk_nnet_propagate(&self->nnet, &x, loglikelihood);

  // Compute log-likelihood
  pk_vector_floor(loglikelihood, 1.0e-20);
  pk_vector_log(loglikelihood);
  pk_vector_addscale(loglikelihood, -1.0f, &self->log_prior);

  pk_vector_destroy(&x);
}
