// Created at 2017-03-13

#include "nnet.h"

#include <assert.h>
#include <math.h>
#include "matrix.h"

bool CheckEq(float a, float b) {
  return fabs(a - b) < 1e-6;
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
  pk_matrix_t x;
  pk_matrix_t y;
  pk_matrix_init(&x, 3, 1);
  pk_matrix_init(&y, 0, 0);
  float x_data[] = {0.3, -0.1, 0.9};
  for (int d = 0; d < x.nrow * x.ncol; ++d) {
    x.data[d] = x_data[d];
  }
  layer->propagate(layer, &x, &y);

  // Check results
  assert(y.nrow == 4 && y.ncol == 1);
  assert(CheckEq(y.data[0], 0.86f));
  assert(CheckEq(y.data[1], 0.63f));
  assert(CheckEq(y.data[2], 0.34f));
  assert(CheckEq(y.data[3], 0.07f));

  // Destory
  layer->destroy(layer);
  pk_matrix_destroy(&x);
  pk_matrix_destroy(&y);
  pk_vector_destroy(&b);
  pk_matrix_destroy(&W);
}


void TestSoftmaxLayer() {
  // Create the Softmax layer
  pk_nnet_layer_softmax_t softmax;
  pk_nnet_layer_t *layer = pk_nnet_layer_softmax_init(&softmax);

  // Propagation
  pk_matrix_t x;
  pk_matrix_t y;
  pk_matrix_init(&x, 4, 1);
  pk_matrix_init(&y, 0, 0);
  float x_data[] = {0.3, -0.1, 0.9, 0.2};
  for (int d = 0; d < x.nrow * x.ncol; ++d) {
    x.data[d] = x_data[d];
  }
  layer->propagate(layer, &x, &y);

  // Check results
  assert(y.nrow == 4);
  assert(y.ncol == 1);
  assert(CheckEq(y.data[0], 0.2274135f));
  assert(CheckEq(y.data[1], 0.15243983f));
  assert(CheckEq(y.data[2], 0.41437442f));
  assert(CheckEq(y.data[3], 0.20577225f));

  // Destory
  if (layer->destroy) layer->destroy(layer);
  pk_matrix_destroy(&x);
  pk_matrix_destroy(&y);
}

void TestReLULayer() {
  // Create the ReLU layer
  pk_nnet_layer_relu_t relu;
  pk_nnet_layer_t *layer = pk_nnet_layer_relu_init(&relu);

  // Propagation
  pk_matrix_t x;
  pk_matrix_t y;
  pk_matrix_init(&x, 4, 1);
  pk_matrix_init(&y, 0, 0);
  float x_data[] = {0.3, -0.1, 0.9, 0.2};
  for (int d = 0; d < x.nrow * x.ncol; ++d) {
    x.data[d] = x_data[d];
  }
  layer->propagate(layer, &x, &y);

  // Check results
  assert(y.nrow == 4);
  assert(y.ncol == 1);
  assert(CheckEq(y.data[0], 0.3f));
  assert(CheckEq(y.data[1], 0.0f));
  assert(CheckEq(y.data[2], 0.9f));
  assert(CheckEq(y.data[3], 0.2f));

  // Destory
  if (layer->destroy) layer->destroy(layer);
  pk_matrix_destroy(&x);
  pk_matrix_destroy(&y);
}

void TestNormalizeLayer() {
  // Create the normalize layer
  pk_nnet_layer_normalize_t normalize;
  pk_nnet_layer_t *layer = pk_nnet_layer_normalize_init(&normalize);

  // Propagation
  pk_matrix_t x;
  pk_matrix_t y;
  pk_matrix_init(&x, 4, 1);
  pk_matrix_init(&y, 0, 0);
  float x_data[] = {0.3, -0.1, 0.9, 0.2};
  for (int d = 0; d < x.nrow * x.ncol; ++d) {
    x.data[d] = x_data[d];
  }
  layer->propagate(layer, &x, &y);

  // Check results
  double sum = 0.0;
  for (int d = 0; d < 4; ++d) {
    sum += y.data[d] * y.data[d];
  }
  assert(fabs(sum - 4.0) < 0.0001);

  // Destory
  if (layer->destroy) layer->destroy(layer);
  pk_matrix_destroy(&x);
  pk_matrix_destroy(&y);
}

int main() {
  TestLinearLayer();
  TestSoftmaxLayer();
  TestReLULayer();
  TestNormalizeLayer();
  return 0;
}
