// Created at 2017-03-13

#include "nnet.h"

#include <assert.h>
#include <math.h>
#include "matrix.h"

using pocketkaldi::Layer;
using pocketkaldi::LinearLayer;
using pocketkaldi::SoftmaxLayer;
using pocketkaldi::ReLULayer;
using pocketkaldi::NormalizeLayer;
using pocketkaldi::Matrix;
using pocketkaldi::SubMatrix;
using pocketkaldi::Vector;
using pocketkaldi::SubVector;

bool CheckEq(float a, float b) {
  return fabs(a - b) < 1e-6;
}

void TestLinearLayer() {
  // Matrix W
  float W_data[] = {
    0.1, 0.8, 0.9,
    0.4, 0.2, 0.7,
    0.2, 0.1, 0.1,
    0.4, 0.3, 0.2
  };
  SubMatrix<float> W(W_data, 4, 3, 3);

  // Vector b
  float b_data[] = {0.1, -0.1, 0.2, -0.2};
  SubVector<float> b(b_data, 4);

  // Create the linear layer
  LinearLayer linear(W, b);

  // Propagation
  // Vector x
  //   0.3 -0.1 0.9
  float x_data[] = {0.3, -0.1, 0.9};
  SubMatrix<float> x(x_data, 1, 3, 3);
  Matrix<float> y;
  linear.Propagate(x, &y);

  // Check results
  assert(y.NumCols() == 4 && y.NumRows() == 1);
  assert(CheckEq(y(0, 0), 0.86f));
  assert(CheckEq(y(0, 1), 0.63f));
  assert(CheckEq(y(0, 2), 0.34f));
  assert(CheckEq(y(0, 3), 0.07f));
}


void TestSoftmaxLayer() {
  // Create the Softmax layer
  SoftmaxLayer softmax;

  float x_data[] = {0.3, -0.1, 0.9, 0.2};
  SubMatrix<float> x(x_data, 1, 4, 4);
  Matrix<float> y;
  softmax.Propagate(x, &y);

  // Check results
  assert(y.NumCols() == 4);
  assert(y.NumRows() == 1);
  assert(CheckEq(y(0, 0), 0.2274135f));
  assert(CheckEq(y(0, 1), 0.15243983f));
  assert(CheckEq(y(0, 2), 0.41437442f));
  assert(CheckEq(y(0, 3), 0.20577225f));
}

void TestReLULayer() {
  // Create the ReLU layer
  ReLULayer relu;

  // Propagation
  float x_data[] = {0.3, -0.1, 0.9, 0.2};
  SubMatrix<float> x(x_data, 1, 4, 4);
  Matrix<float> y;
  relu.Propagate(x, &y);

  // Check results
  assert(y.NumCols() == 4);
  assert(y.NumRows() == 1);
  assert(CheckEq(y(0, 0), 0.3f));
  assert(CheckEq(y(0, 1), 0.0f));
  assert(CheckEq(y(0, 2), 0.9f));
  assert(CheckEq(y(0, 3), 0.2f));
}

void TestNormalizeLayer() {
  // Create the normalize layer
  NormalizeLayer normalize;

  // Propagation
  float x_data[] = {0.3, -0.1, 0.9, 0.2};
  SubMatrix<float> x(x_data, 1, 4, 4);
  Matrix<float> y;
  normalize.Propagate(x, &y);

  // Check results
  double sum = 0.0;
  for (int d = 0; d < 4; ++d) {
    sum += y(0, d) * y(0, d);
  }
  assert(fabs(sum - 4.0) < 0.0001);
}

int main() {
  TestLinearLayer();
  TestSoftmaxLayer();
  TestReLULayer();
  TestNormalizeLayer();
  return 0;
}
