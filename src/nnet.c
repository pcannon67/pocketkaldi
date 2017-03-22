// Created at 2017-03-13

#include "nnet.h"

#include <assert.h>
#include <math.h>

void pk_nnet_layer_linear_destroy(pk_nnet_layer_t *self) {
  pk_nnet_layer_linear_t *self_linear = (pk_nnet_layer_linear_t *)self;

  pk_matrix_destroy(&self_linear->W);
  pk_vector_destroy(&self_linear->b);
  self_linear->base.propagate = NULL;
  self_linear->base.destroy = NULL;
}

void pk_nnet_layer_linear_propagate(
    const pk_nnet_layer_t *self,
    const pk_vector_t *in,
    pk_vector_t *out) {
  pk_nnet_layer_linear_t *self_linear = (pk_nnet_layer_linear_t *)self;
  pk_vector_dotmat(in, &self_linear->W, out);
  pk_vector_add(out, &self_linear->b);
}

pk_nnet_layer_t *pk_nnet_layer_linear_init(
    pk_nnet_layer_linear_t *self,
    const pk_matrix_t *W,
    const pk_vector_t *b) {
  assert(b->dim == W->ncol && "linear layer: dimension mismatch in W and b");

  pk_matrix_init(&self->W, 0, 0);
  pk_vector_init(&self->b, 0, NAN);
  pk_matrix_copy(&self->W, W);
  pk_vector_copy(&self->b, b);

  self->base.destroy = &pk_nnet_layer_linear_destroy;
  self->base.propagate = &pk_nnet_layer_linear_propagate;

  return (pk_nnet_layer_t *)self;
}

void pk_nnet_layer_softmax_propagate(
    const pk_nnet_layer_t *self,
    const pk_vector_t *in,
    pk_vector_t *out) {
  PK_UNUSED(self);

  pk_vector_resize(out, in->dim);
  float sum = 0.0f;
  for (int i = 0; i < in->dim; ++i) {
    float exp_d = expf(in->data[i]);
    out->data[i] = exp_d;
    sum += exp_d;
  }
  for (int i = 0; i < out->dim; ++i) {
    out->data[i] /= sum;
  }
}

pk_nnet_layer_t *pk_nnet_layer_softmax_init(pk_nnet_layer_softmax_t *self) {
  self->base.destroy = NULL;
  self->base.propagate = &pk_nnet_layer_softmax_propagate;

  return (pk_nnet_layer_t *)self;
}

void pk_nnet_layer_relu_propagate(
    const pk_nnet_layer_t *self,
    const pk_vector_t *in,
    pk_vector_t *out) {
  PK_UNUSED(self);

  pk_vector_resize(out, in->dim);
  float sum = 0.0f;
  for (int i = 0; i < in->dim; ++i) {
    out->data[i] = in->data[i] > 0.0f ? in->data[i] : 0.0f;
  }
}

pk_nnet_layer_t *pk_nnet_layer_relu_init(pk_nnet_layer_relu_t *self) {
  self->base.destroy = NULL;
  self->base.propagate = &pk_nnet_layer_relu_propagate;

  return (pk_nnet_layer_t *)self;
}

void pk_nnet_layer_normalize_propagate(
    const pk_nnet_layer_t *self,
    const pk_vector_t *in,
    pk_vector_t *out) {
  PK_UNUSED(self);
  float D = in->dim;

  pk_vector_resize(out, in->dim);
  double squared_sum = 0.0;
  for (int i = 0; i < in->dim; ++i) {
    squared_sum += in->data[i] * in->data[i];
  }
  float scale = (float)sqrt(D / squared_sum);
  for (int i = 0; i < in->dim; ++i) {
    out->data[i] = in->data[i] * scale;
  }
}

pk_nnet_layer_t *pk_nnet_layer_normalize_init(pk_nnet_layer_normalize_t *self) {
  self->base.destroy = NULL;
  self->base.propagate = &pk_nnet_layer_normalize_propagate;

  return (pk_nnet_layer_t *)self;
}

void pk_nnet_layer_add_propagate(
    const pk_nnet_layer_t *self,
    const pk_vector_t *in,
    pk_vector_t *out) {
  pk_nnet_layer_add_t *add_layer = (pk_nnet_layer_add_t *)self;
  assert(in->dim == add_layer->b.dim && "add_layer: vector size mismatch");

  pk_vector_resize(out, in->dim);
  for (int i = 0; i < in-> dim; ++i) {
    out->data[i] = in->data[i] + add_layer->scale * add_layer->b.data[i];
  }
}

void pk_nnet_layer_add_destroy(pk_nnet_layer_t *self) {
  pk_nnet_layer_add_t *add_layer = (pk_nnet_layer_add_t *)self;

  pk_vector_destroy(&add_layer->b);
  add_layer->base.propagate = NULL;
  add_layer->base.destroy = NULL;
  add_layer->scale = 0.0f;
}

pk_nnet_layer_t *pk_nnet_layer_add_init(
    pk_nnet_layer_add_t *self,
    float scale,
    const pk_vector_t *b) {
  self->base.destroy = &pk_nnet_layer_add_destroy;
  self->base.propagate = &pk_nnet_layer_add_propagate;

  // Copy and initialize b
  pk_vector_init(&self->b, 0, NAN);
  pk_vector_copy(&self->b, b);

  self->scale = scale;
}

void pk_nnet_init(pk_nnet_t *self) {
  self->layers = NULL;
  self->num_layers = 0;
}

// Reads a layer from fd, returns a pointer to that layer. And to free this 
// pointer, it should call pointer->destroy(pointer) first when pointer != NULL.
// Then call free to free it. When failed, status->ok == false
static pk_nnet_layer_t *read_layer(pk_readable_t *fd, pk_status_t *status) {
  pk_nnet_layer_t *layer = NULL;
  pk_bytebuffer_t bytebuffer;
  pk_matrix_t W;
  pk_vector_t b;
  pk_bytebuffer_init(&bytebuffer);
  pk_matrix_init(&W, 0, 0);
  pk_vector_init(&b, 0, NAN);

  int section_size = pk_readable_readsectionhead(
      fd,
      PK_NNET_LAYER_SECTION,
      status);
  if (!status->ok) goto read_layer_failed;
  pk_bytebuffer_reset(&bytebuffer, section_size);

  pk_readable_readbuffer(fd, &bytebuffer, status);
  if (!status->ok) goto read_layer_failed;
  int layer_type = pk_bytebuffer_readint32(&bytebuffer);

  // Check size of this section
  int expected_size = 4;
  if (layer_type == PK_NNET_ADD_LAYER) expected_size = 8;
  if (expected_size != section_size) {
    PK_STATUS_CORRUPTED(
        status,
        "read_layer: section_size == %d expected, but %d found (%s)",
        expected_size,
        section_size,
        fd->filename);
    goto read_layer_failed;
  }

  // Read additional parameters and initialize layer
  pk_nnet_layer_linear_t *linear_layer = NULL;
  pk_nnet_layer_relu_t *relu_layer = NULL;
  pk_nnet_layer_normalize_t *norm_layer = NULL;
  pk_nnet_layer_softmax_t *softmax_layer = NULL;
  pk_nnet_layer_add_t *add_layer = NULL;
  float scale = 0.0f;
  switch (layer_type) {
  case PK_NNET_LINEAR_LAYER:
    pk_matrix_read(&W, fd, status);
    if (!status->ok) goto read_layer_failed;

    pk_vector_read(&b, fd, status);
    if (!status->ok) goto read_layer_failed;

    linear_layer = (pk_nnet_layer_linear_t *)malloc(
        sizeof(pk_nnet_layer_linear_t));
    layer = pk_nnet_layer_linear_init(linear_layer, &W, &b);
    break;
  case PK_NNET_RELU_LAYER:
    relu_layer = (pk_nnet_layer_relu_t *)malloc(
        sizeof(pk_nnet_layer_relu_t));
    layer = pk_nnet_layer_relu_init(relu_layer);
    break;
  case PK_NNET_NORMALIZE_LAYER:
    norm_layer = (pk_nnet_layer_normalize_t *)malloc(
        sizeof(pk_nnet_layer_normalize_t));
    layer = pk_nnet_layer_normalize_init(norm_layer);
    break;
  case PK_NNET_SOFTMAX_LAYER:
    softmax_layer = (pk_nnet_layer_softmax_t *)malloc(
        sizeof(pk_nnet_layer_softmax_t));
    layer = pk_nnet_layer_softmax_init(softmax_layer);
    break;
  case PK_NNET_ADD_LAYER:
    scale = pk_bytebuffer_readfloat(&bytebuffer);
    pk_vector_read(&b, fd, status);
    if (!status->ok) goto read_layer_failed;

    add_layer = (pk_nnet_layer_add_t *)malloc(
        sizeof(pk_nnet_layer_add_t));
    layer = pk_nnet_layer_add_init(add_layer, scale, &b);
    break;
  default:
    PK_STATUS_CORRUPTED(
        status,
        "read_layer: unexpected layer_type %d (%s)",
        read_layer,
        fd->filename);
    goto read_layer_failed;
  }

  if (false) {
read_layer_failed:
    if (layer != NULL) {
      if (layer->destroy) layer->destroy(layer);
      layer->destroy = NULL;
      layer->propagate = NULL;
      free(layer);
    }
    layer = NULL;
  }

  pk_vector_destroy(&b);
  pk_matrix_destroy(&W);
  pk_bytebuffer_destroy(&bytebuffer);

  return layer;
}

void pk_nnet_read(pk_nnet_t *self, pk_readable_t *fd, pk_status_t *status) {
  pk_bytebuffer_t bytebuffer;
  pk_bytebuffer_init(&bytebuffer);

  int section_size = pk_readable_readsectionhead(
      fd,
      PK_NNET_SECTION,
      status);
  if (!status->ok) goto pk_nnet_read_failed;
  pk_bytebuffer_reset(&bytebuffer, section_size);

  pk_readable_readbuffer(fd, &bytebuffer, status);
  if (!status->ok) goto pk_nnet_read_failed;
  int num_layers = pk_bytebuffer_readint32(&bytebuffer);

  // Initialize the layer pointers
  self->layers = (pk_nnet_layer_t **)malloc(
      sizeof(pk_nnet_layer_t) * num_layers);
  for (int layer_idx = 0; layer_idx < num_layers; ++layer_idx) {
    self->layers[layer_idx] = NULL;
  }
  self->num_layers = num_layers;

  // Read each layers
  for (int layer_idx = 0; layer_idx < num_layers; ++layer_idx) {
    pk_nnet_layer_t *layer = read_layer(fd, status);
    if (!status->ok) goto pk_nnet_read_failed;
    self->layers[layer_idx] = layer;
  }

  if (false) {
pk_nnet_read_failed:
    pk_nnet_destroy(self);
  }
  
  pk_bytebuffer_destroy(&bytebuffer);
}

void pk_nnet_propagate(
    const pk_nnet_t *self,
    const pk_vector_t *in,
    pk_vector_t *out) {
  pk_vector_t *x = NULL,
              *y = NULL,
              *t = NULL,
              forward_data[2];
  pk_vector_init(&forward_data[0], 0, NAN);
  pk_vector_init(&forward_data[1], 0, NAN);
  x = &forward_data[0];
  y = &forward_data[1];

  pk_vector_copy(x, in);
  for (int layer_idx = 0; layer_idx < self->num_layers; ++layer_idx) {
    self->layers[layer_idx]->propagate(self->layers[layer_idx], x, y);

    // Swap x and y
    t = x;
    x = y;
    y = t;
  }

  pk_vector_copy(out, x);
  pk_vector_destroy(&forward_data[0]);
  pk_vector_destroy(&forward_data[1]);
}

void pk_nnet_destroy(pk_nnet_t *self) {
  for (int layer_idx = 0; layer_idx < self->num_layers; ++layer_idx) {
    pk_nnet_layer_t *layer = self->layers[layer_idx];
    if (layer == NULL) continue;

    if (layer->destroy) layer->destroy(layer);
    layer->destroy = NULL;
    layer->propagate = NULL;
    free(layer);
    self->layers[layer_idx] = NULL;
  }
  self->num_layers = 0;
}
