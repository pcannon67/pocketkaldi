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
  assert(self);  // self is an unused parameter

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
  assert(self);  // self is an unused parameter

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
  assert(self);  // self is an unused parameter

  pk_vector_resize(out, in->dim);
  double squared_sum = 0.0;
  for (int i = 0; i < in->dim; ++i) {
    squared_sum += in->data[i] * in->data[i];
  }
  float invstd = (float)(1.0 / sqrt(squared_sum));
  for (int i = 0; i < in-> dim; ++i) {
    out->data[i] = in->data[i] * invstd;
  }
}

pk_nnet_layer_t *pk_nnet_layer_normalize_init(pk_nnet_layer_normalize_t *self) {
  self->base.destroy = NULL;
  self->base.propagate = &pk_nnet_layer_normalize_propagate;

  return (pk_nnet_layer_t *)self;
}
