// Created at 2017-03-13

#ifndef POCKETKALDI_NNET_H_
#define POCKETKALDI_NNET_H_

#include "matrix.h"

typedef struct pk_nnet_layer_t {
  void (*propagate)(
      const struct pk_nnet_layer_t *self,
      const pk_vector_t *in,
      pk_vector_t *out);
  void (*destroy)(struct pk_nnet_layer_t *self);
} pk_nnet_layer_t;

// Linear layer: x^T dot W + b
typedef struct pk_nnet_layer_linear_t {
  pk_nnet_layer_t base;
  pk_matrix_t W;
  pk_vector_t b;
} pk_nnet_layer_linear_t;

// Softmax layer
typedef struct pk_nnet_layer_softmax_t {
  pk_nnet_layer_t base;
} pk_nnet_layer_softmax_t;

// ReLU layer
typedef struct pk_nnet_layer_relu_t {
  pk_nnet_layer_t base;
} pk_nnet_layer_relu_t;

// Normalize layer: ensure root-mean-square equals 1.0
typedef struct pk_nnet_layer_normalize_t {
  pk_nnet_layer_t base;
} pk_nnet_layer_normalize_t;

// Initialize the linear layer with parameter W and b. It just copies the values
// from W and b. Then return a pointer of pk_nnet_layer_t, which points to the
// same address as self.
POCKETKALDI_EXPORT
pk_nnet_layer_t *pk_nnet_layer_linear_init(
    pk_nnet_layer_linear_t *self,
    const pk_matrix_t *W,
    const pk_vector_t *b);

// Destroy the linear layer. It just destroys W and b, and doesn't free the
// pointer self
POCKETKALDI_EXPORT
void pk_nnet_layer_linear_destroy(pk_nnet_layer_t *self);

// Propagate through the linear layer
POCKETKALDI_EXPORT
void pk_nnet_layer_linear_propagate(
    const pk_nnet_layer_t *self,
    const pk_vector_t *in,
    pk_vector_t *out);

// Initialize the softmax layer
POCKETKALDI_EXPORT
pk_nnet_layer_t *pk_nnet_layer_softmax_init(pk_nnet_layer_softmax_t *self);

// Propagate through the softmax layer
POCKETKALDI_EXPORT
void pk_nnet_layer_softmax_propagate(
    const pk_nnet_layer_t *,
    const pk_vector_t *in,
    pk_vector_t *out);

// Initialize the ReLU layer
POCKETKALDI_EXPORT
pk_nnet_layer_t *pk_nnet_layer_relu_init(pk_nnet_layer_relu_t *self);

// Propagate through the ReLU layer
POCKETKALDI_EXPORT
void pk_nnet_layer_relu_propagate(
    const pk_nnet_layer_t *,
    const pk_vector_t *in,
    pk_vector_t *out);

// Initialize the Normalize layer
POCKETKALDI_EXPORT
pk_nnet_layer_t *pk_nnet_layer_normalize_init(pk_nnet_layer_normalize_t *self);

// Propagate through the Normalize layer
POCKETKALDI_EXPORT
void pk_nnet_layer_normalize_propagate(
    const pk_nnet_layer_t *,
    const pk_vector_t *in,
    pk_vector_t *out);

#endif
