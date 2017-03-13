// Created at 2017-03-13

#include "nnet.h"

#include <assert.h>
#include <math.h>
#include "matrix.h"

bool CheckEq(float a, float b) {
  return abs(a - b) < 1e-6;
}

void TestLinearLayer() {
  pk_matrix_t W;
  pk_vector_t b;

  // Matrix W
  //   0.1 0.4 0.2 0.4
  //   0.8 0.2 0.1 0.3
  //   0.9 0.7 0.1 0.2
  pk_matrix_init(&W, 3, 4);
  float W_data[] = {
    0.1, 0.8, 0.9,
    0.4, 0.2, 0.7,
    0.2, 0.1, 0.1,
    0.4, 0.3, 0.2
  };
  for (int i = 0; i < W.ncol * W.nrow; ++i) {
    W.data[i] = W_data[i];
  }

  // Vector b
  //   0.1 -0.1 0.2 -0.2
  pk_vector_init(&b, 4, NAN);
  float b_data[] = {0.1, -0.1, 0.2, -0.2};
  for (int d = 0; d < b.dim; ++d) {
    b.data[d] = b_data[d];
  }

  // Create the linear layer
  pk_nnet_layer_linear_t linear;
  pk_nnet_layer_t *layer = pk_nnet_layer_linear_init(&linear, &W, &b);

  // Propagation
  // Vector x
  //   0.3 -0.1 0.9
  pk_vector_t x;
  pk_vector_t y;
  pk_vector_init(&x, 3, NAN);
  pk_vector_init(&y, 0, NAN);
  float x_data[] = {0.3, -0.1, 0.9};
  for (int d = 0; d < x.dim; ++d) {
    x.data[d] = x_data[d];
  }
  layer->propagate(layer, &x, &y);

  // Check results
  assert(y.dim == 4);
  assert(CheckEq(y.data[0], 0.86f));
  assert(CheckEq(y.data[1], 0.63f));
  assert(CheckEq(y.data[2], 0.34f));
  assert(CheckEq(y.data[3], 0.07f));

  // Destory
  layer->destroy(layer);
  pk_vector_destroy(&x);
  pk_vector_destroy(&y);
  pk_vector_destroy(&b);
  pk_matrix_destroy(&W);
}

int main() {
  TestLinearLayer();
  return 0;
}
