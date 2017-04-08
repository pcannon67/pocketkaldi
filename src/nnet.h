// Created at 2017-03-13

#ifndef POCKETKALDI_NNET_H_
#define POCKETKALDI_NNET_H_

#include "matrix.h"
#include "util.h"

#define PK_NNET_SECTION "NNT0"
#define PK_NNET_LAYER_SECTION "LAY0"

// Layer types
#define PK_NNET_LINEAR_LAYER 0
#define PK_NNET_RELU_LAYER 1
#define PK_NNET_NORMALIZE_LAYER 2
#define PK_NNET_SOFTMAX_LAYER 3
#define PK_NNET_ADD_LAYER 4
#define PK_NNET_MUL_LAYER 5

typedef struct pk_nnet_layer_t {
  void (*propagate)(
      const struct pk_nnet_layer_t *self,
      const pk_matrix_t *in,
      pk_matrix_t *out);
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

// Normalize layer: ensure y * y^T == D
typedef struct pk_nnet_layer_normalize_t {
  pk_nnet_layer_t base;
} pk_nnet_layer_normalize_t;

typedef struct pk_nnet_t {
  pk_nnet_layer_t **layers;
  int num_layers;
} pk_nnet_t;

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
    const pk_matrix_t *in,
    pk_matrix_t *out);

// Initialize the softmax layer
POCKETKALDI_EXPORT
pk_nnet_layer_t *pk_nnet_layer_softmax_init(pk_nnet_layer_softmax_t *self);

// Propagate through the softmax layer
POCKETKALDI_EXPORT
void pk_nnet_layer_softmax_propagate(
    const pk_nnet_layer_t *,
    const pk_matrix_t *in,
    pk_matrix_t *out);

// Initialize the ReLU layer
POCKETKALDI_EXPORT
pk_nnet_layer_t *pk_nnet_layer_relu_init(pk_nnet_layer_relu_t *self);

// Propagate through the ReLU layer
POCKETKALDI_EXPORT
void pk_nnet_layer_relu_propagate(
    const pk_nnet_layer_t *,
    const pk_matrix_t *in,
    pk_matrix_t *out);

// Initialize the Normalize layer
POCKETKALDI_EXPORT
pk_nnet_layer_t *pk_nnet_layer_normalize_init(pk_nnet_layer_normalize_t *self);

// Propagate through the Normalize layer
POCKETKALDI_EXPORT
void pk_nnet_layer_normalize_propagate(
    const pk_nnet_layer_t *,
    const pk_matrix_t *in,
    pk_matrix_t *out);

// Initialize the neural network
POCKETKALDI_EXPORT
void pk_nnet_init(pk_nnet_t *self);

// Read the neural network from fd. If failed, status->ok == false
POCKETKALDI_EXPORT
void pk_nnet_read(pk_nnet_t *self, pk_readable_t *fd, pk_status_t *status);

// Propagate through the netral network
POCKETKALDI_EXPORT
void pk_nnet_propagate(
    const pk_nnet_t *self,
    const pk_matrix_t *in,
    pk_matrix_t *out);

// Destroy the nnet
POCKETKALDI_EXPORT
void pk_nnet_destroy(pk_nnet_t *self);

#endif
