// Created at 2017-03-13

#include "nnet.h"

#include <assert.h>
#include <math.h>

namespace pocketkaldi {

LinearLayer::LinearLayer(const pk_matrix_t *W, const pk_vector_t *b) {
  assert(b->dim == W->ncol && "linear layer: dimension mismatch in W and b");
  pk_matrix_init(&W_, 0, 0);
  pk_vector_init(&b_, 0, NAN);
  pk_matrix_copy(&W_, W);
  pk_vector_copy(&b_, b);
}

LinearLayer::~LinearLayer() {
  pk_matrix_destroy(&W_);
  pk_vector_destroy(&b_);
}

void LinearLayer::Propagate(const pk_matrix_t *in, pk_matrix_t *out) const {
  pk_matrix_resize(out, W_.ncol, in->ncol);
  pk_matrix_matmat(&W_, in, out);
  for (int col_idx = 0; col_idx < out->ncol; ++col_idx) {
    pk_vector_t col = pk_matrix_getcol(out, col_idx);
    pk_vector_add(&col, &b_);
  }
}

void SoftmaxLayer::Propagate(const pk_matrix_t *in, pk_matrix_t *out) const {
  pk_matrix_resize(out, in->nrow, in->ncol);
  for (int col_idx = 0; col_idx < in->ncol; ++col_idx) {
    pk_vector_t col_in = pk_matrix_getcol(in, col_idx);
    pk_vector_t col_out = pk_matrix_getcol(out, col_idx);

    float sum = 0.0f;
    for (int i = 0; i < col_in.dim; ++i) {
      float exp_d = expf(col_in.data[i]);
      col_out.data[i] = exp_d;
      sum += exp_d;
    }
    for (int i = 0; i < col_out.dim; ++i) {
      col_out.data[i] /= sum;
    }
  }
}

void ReLULayer::Propagate(const pk_matrix_t *in, pk_matrix_t *out) const {
  pk_matrix_resize(out, in->nrow, in->ncol);
  for (int col_idx = 0; col_idx < in->ncol; ++col_idx) {
    pk_vector_t col_in = pk_matrix_getcol(in, col_idx);
    pk_vector_t col_out = pk_matrix_getcol(out, col_idx);

    for (int i = 0; i < col_in.dim; ++i) {
      col_out.data[i] = col_in.data[i] > 0.0f ? col_in.data[i] : 0.0f;
    }
  }
}

void NormalizeLayer::Propagate(const pk_matrix_t *in, pk_matrix_t *out) const {
  float D = in->nrow;
  pk_matrix_resize(out, in->nrow, in->ncol);
  for (int col_idx = 0; col_idx < in->ncol; ++col_idx) {
    pk_vector_t col_in = pk_matrix_getcol(in, col_idx);
    pk_vector_t col_out = pk_matrix_getcol(out, col_idx);

    double squared_sum = 0.0;
    for (int i = 0; i < col_in.dim; ++i) {
      squared_sum += col_in.data[i] * col_in.data[i];
    }
    float scale = (float)sqrt(D / squared_sum);
    for (int i = 0; i < col_in.dim; ++i) {
      col_out.data[i] = col_in.data[i] * scale;
    }
  }
}

Nnet::Nnet() {
}

Status Nnet::ReadLayer(pk_readable_t *fd) {
  pk_status_t status;
  pk_status_init(&status);

  float scale;
  int layer_type;
  int expected_size;

  pk_bytebuffer_t bytebuffer;
  pk_matrix_t W;
  pk_vector_t b;
  pk_bytebuffer_init(&bytebuffer);
  pk_matrix_init(&W, 0, 0);
  pk_vector_init(&b, 0, NAN);

  int section_size = pk_readable_readsectionhead(
      fd,
      PK_NNET_LAYER_SECTION,
      &status);
  if (!status.ok) goto read_layer_failed;
  pk_bytebuffer_reset(&bytebuffer, section_size);

  pk_readable_readbuffer(fd, &bytebuffer, &status);
  if (!status.ok) goto read_layer_failed;
  layer_type = pk_bytebuffer_readint32(&bytebuffer);

  // Check size of this section
  expected_size = 4;
  if (expected_size != section_size) {
    PK_STATUS_CORRUPTED(
        &status,
        "read_layer: section_size == %d expected, but %d found (%s)",
        expected_size,
        section_size,
        fd->filename);
    goto read_layer_failed;
  }

  // Read additional parameters and initialize layer
  scale = 0.0f;
  LinearLayer *layer;
  switch (layer_type) {
  case Layer::kLinear:
    pk_matrix_read(&W, fd, &status);
    pk_vector_read(&b, fd, &status);

    layers_.emplace_back(new LinearLayer(&W, &b));
    break;
  case Layer::kReLU:
    layers_.emplace_back(new ReLULayer());
    break;
  case Layer::kNormalize:
    layers_.emplace_back(new NormalizeLayer());
    break;
  case Layer::kSoftmax:
    layers_.emplace_back(new SoftmaxLayer());
    break;
  default:
    PK_STATUS_CORRUPTED(
        &status,
        "read_layer: unexpected layer_type %d (%s)",
        layer_type,
        fd->filename);
    goto read_layer_failed;
  }

  if (false) {
read_layer_failed:
    return Status::IOError(status.message);
  }

  pk_vector_destroy(&b);
  pk_matrix_destroy(&W);
  pk_bytebuffer_destroy(&bytebuffer);

  return Status::OK();
}

Status Nnet::Read(pk_readable_t *fd) {
  Status vn_status;
  pk_status_t status;
  int num_layers;
  int section_size;
  
  pk_bytebuffer_t bytebuffer;
  pk_bytebuffer_init(&bytebuffer);
  pk_status_init(&status);

  section_size = pk_readable_readsectionhead(
      fd,
      PK_NNET_SECTION,
      &status);
  if (!status.ok) return Status::IOError(status.message);
  pk_bytebuffer_reset(&bytebuffer, section_size);

  pk_readable_readbuffer(fd, &bytebuffer, &status);
  if (!status.ok) return Status::IOError(status.message);
  num_layers = pk_bytebuffer_readint32(&bytebuffer);

  // Read each layers
  for (int layer_idx = 0; layer_idx < num_layers; ++layer_idx) {
    vn_status = ReadLayer(fd);
    if (!vn_status.ok()) goto pk_nnet_read_failed;
  }

  if (false) {
pk_nnet_read_failed:
    return Status::IOError(status.message);
  }
  
  pk_bytebuffer_destroy(&bytebuffer);
  return Status::OK();
}

void Nnet::Propagate(const pk_matrix_t *in, pk_matrix_t *out) const {
  pk_matrix_t *x = NULL,
              *y = NULL,
              *t = NULL,
              forward_data[2];
  pk_matrix_init(&forward_data[0], 0, 0);
  pk_matrix_init(&forward_data[1], 0, 0);
  x = &forward_data[0];
  y = &forward_data[1];

  pk_matrix_copy(x, in);
  for (const std::unique_ptr<Layer> &layer : layers_) {
    layer->Propagate(x, y);

    // Swap x and y
    t = x;
    x = y;
    y = t;
  }

  pk_matrix_copy(out, x);
  pk_matrix_destroy(&forward_data[0]);
  pk_matrix_destroy(&forward_data[1]);
}

}  // namespace pocketkaldi
